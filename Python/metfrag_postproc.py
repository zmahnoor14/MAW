#!/usr/bin/env python
# coding: utf-8

#!/usr/bin/env python
#make executable in bash chmod +x PyRun

# Libraries

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem import PandasTools

import sys

import pandas as pd

import pubchempy as pcp
import numpy as np


import os
import glob
import re


def metfrag_postproc(input_dir, input_tablecsv, sl= True):
    def isNaN(string):
        return string != string
    
    """metfrag_postproc function re-checks the Suspect list, if present 
    or given as a parameter, whether the candidates have a high 
    similarity with compounds in Suspect List. It also calculates the 
    Maximum Common Substructure (MCSS). This function adds top candidates
    from PubChem and KEGG as these two databases are used with MetFrag

    Parameters:
    input_dir (str): This is the input directory where all the .mzML 
    files and their respective result directories are stored. For this 
    function this directory must contain a csv file that has a column 
    named "SMILES".
    
    input_tablecsv (str): This is the table in csv format (defined in R), 
    which stores a csv table containing columns "mzml_files", which 
    contains liat of all input files with their relative paths, second
    column is "ResultFileName" which is a list of the corresponding
    result relative directories to each mzml files. Lastly, "file_id", 
    contains a file directory. This table will be used to read the 
    MetFrag csv files

    Returns:
    csv: a result file with additional columns such as those for suspect
    list if one is used. It also adds columns on MCSS., named as 
    "input_dir/ResultFileName/insilico/MetFragResults.csv". It 
    contains columns for KEGG and PubChem
    
    
    Usage:
    metfrag_postproc(input_dir = "/user/project/", 
    input_table = "/user/project/suspectlist.csv", sl = True, slistcsv)


    """
    
    # Describe the heavy atoms to be considered for MCSS
    heavy_atoms = ['C', 'N', 'P', 'O', 'S']
    
    input_table = pd.read_csv(input_tablecsv)
    
    for m, row in input_table.iterrows():
        
        # read SIRIUS results:
        
        #siriusResults = pd.read_csv(input_dir + (input_table['ResultFileNames'][m] + '/insilico/SiriusResults.csv'))
    
        # Result directory
        result = input_dir + (input_table['ResultFileNames'][m] + 
                                 '/insilico/MetFrag').replace("./", "")

        # list of all the csv files in the result directory result_dir/inislico/MetFrag/
        files_met = (glob.glob(result+'/*.csv'))

        # read the csv file that contains all the features from the input .mzml file
        file1  = pd.read_csv(input_dir + (input_table['ResultFileNames'][m] + '/insilico/MS1DATA.csv').replace("./", ""))
    
        # for each feature in the MS1DATA.csv file
        for i, row in file1.iterrows():
        
            # take id as a pattern to differentiate between different ids
            pattern = file1.loc[i, "id_X"]
        
            #check which of the csv result files have the same pattern in their names
            results = [i for i in files_met if pattern in i]
        
            # find which of the files with that id have KEGG in their names,
            KEGG = [i for i in results if "KEGG" in i]
        
            # if kegg present in the name
            if KEGG:
            
                # read the KEGG csv file for that feature
                KEGG_file = pd.read_csv((KEGG)[0])
            
                # if the KEGG file isn't empty
                if len(KEGG_file)>0:
                
                    # extract only the columns with >0.75 score
                    KEGG_file = KEGG_file.drop(KEGG_file[KEGG_file.Score < 0.98].index)
                    
                    #s_best_kg = []
                    #for kg, rows in KEGG_file.iterrows():
                        #kg_smiles = Chem.MolToSmiles(Chem.MolFromInchi(KEGG_file["InChI"][kg]))
                        #SSmsk = [Chem.MolFromSmiles(kg_smiles), Chem.MolFromSmiles(siriusResults["SMILES"][0])]
                        #SSfpsk = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in SSmsk]
                        #SStn2k = DataStructs.FingerprintSimilarity(SSfpsk[0],SSfpsk[1])
                        #s_best_kg.append(SStn2k)
                    #index_kg = np.argmax(s_best_kg)
                        
                    # add the relevavnt information to the original MS1DATA csv
                    file1.loc[i, 'KG_ID'] = KEGG_file.loc[0, 'Identifier']
                    file1.loc[i, 'KG_Name'] = KEGG_file.loc[0, 'CompoundName']
                    file1.loc[i, 'KG_Formula'] = KEGG_file.loc[0, 'MolecularFormula']
                    file1.loc[i, 'KG_expPeaks'] = KEGG_file.loc[0, 'NoExplPeaks']
                    file1.loc[i, 'KG_SMILES'] = Chem.MolToSmiles(Chem.MolFromInchi(KEGG_file["InChI"][0]))
                    file1.loc[i, 'KG_Score'] = KEGG_file.loc[0, 'Score']
                    if sl:
                        file1.loc[i, 'KGSL_Score'] = KEGG_file.loc[0, 'SuspectListScore']
                    file1.loc[i, 'KG_file'] = KEGG[0]
                
                    #create empty list of KEGG top smiles
                    Kegg_smiles = []
                
                    # extract only the InChI of the top 5
                    for j in KEGG_file["InChI"][0:5].tolist():
                        # convert the InChI to SMILES
                        mol = Chem.MolToSmiles(Chem.MolFromInchi(j))
                        mol2 = Chem.MolFromSmiles(mol)
                        Kegg_smiles.append(mol2)
                    # if there are more than 1 top smiles
                    if len(Kegg_smiles) > 1:
                        #calculate the MCSS
                        res = rdFMCS.FindMCS(Kegg_smiles)
                        sm_res = Chem.MolToSmiles(Chem.MolFromSmarts(res.smartsString))
                        # if there are atleast 3 heavy atoms in the MCSS, then add it to the result file
                        elem = [ele for ele in heavy_atoms if(ele in sm_res)]
                        if elem and len(sm_res)>=3:
                            file1.loc[i, 'KG_MCSSstring'] = res.smartsString
                            file1.loc[i, 'KG_MCSS_SMILES'] = Chem.MolToSmiles(Chem.MolFromSmarts(res.smartsString))
                            
            #start here for PubChem; find which of the files with that id have PubChem in their names,
            PubChem = [i for i in results if "PubChem" in i]
            
            if PubChem:

                PubChem_file = pd.read_csv(PubChem[0])
                
                # if more candidates
                if len(PubChem_file)>0:
                    
                    # take the ones with more than 0.80 score
                    PubChem_file = PubChem_file.drop(PubChem_file[PubChem_file.Score < 0.80].index)
                    #s_best_pb = []
                    #for pb, rows in PubChem_file.iterrows():
                        #SSmsp = [Chem.MolFromSmiles(PubChem_file["SMILES"][pb]), Chem.MolFromSmiles(siriusResults["SMILES"][0])]
                        #SSfpsp = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in SSmsp]
                        #SStn2p = DataStructs.FingerprintSimilarity(SSfpsp[0],SSfpsp[1])
                        #s_best_pb.append(SStn2p)
                    #index_pb = np.argmax(s_best_pb)
                    # add the relavnt information to the original MS1DATA csv
                    file1.loc[i, 'PC_ID'] = PubChem_file.loc[0, 'Identifier']
                    file1.loc[i, 'PC_Name'] = PubChem_file.loc[0, 'IUPACName']
                    file1.loc[i, 'PC_Formula'] = PubChem_file.loc[0, 'MolecularFormula']
                    file1.loc[i, 'PC_expPeaks'] = PubChem_file.loc[0, 'NoExplPeaks']
                    file1.loc[i, 'PC_SMILES'] = PubChem_file["SMILES"][0]
                    file1.loc[i, 'PC_Score'] = PubChem_file["Score"][0]
                    if sl:
                        file1.loc[i, 'PCSL_Score'] = PubChem_file.loc[0, 'SuspectListScore']
                    file1.loc[i, 'PC_file'] = PubChem[0]
                    
                    # empty object
                    Pubchem_smiles = []
                    
                    # extract only the SMILES of the top 5
                    for j in PubChem_file["SMILES"][0:5].tolist():
                        
                        # Concert smiles to mol
                        sm2 = Chem.MolFromSmiles(j)
                        # store mol in Pubchem_smiles
                        Pubchem_smiles.append(sm2)
                    
                    if len(Pubchem_smiles) > 1:
                        # calculate MCSS
                        res2 = rdFMCS.FindMCS(Pubchem_smiles)
                        sm_res = Chem.MolToSmiles(Chem.MolFromSmarts(res2.smartsString))
                        # If atleast 3 heavy atoms present
                        elem = [ele for ele in heavy_atoms if(ele in sm_res)]
                        if elem and len(sm_res)>=3:
                            file1.loc[i, 'PC_MCSSstring']= res2.smartsString
                            file1.loc[i, 'PC_MCSS_SMILES'] = Chem.MolToSmiles(Chem.MolFromSmarts(res2.smartsString))
        file1.to_csv(input_dir + (input_table['ResultFileNames'][m] + '/insilico/MetFragResults.csv').replace("./", ""))
    
metfrag_postproc(sys.argv[1], sys.argv[2], sys.argv[3])

