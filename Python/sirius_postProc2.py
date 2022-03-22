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

def sirius_postProc2(input_dir, input_tablecsv, slistcsv ,sl = True):
    
    def isNaN(string):
        return string != string
    """sirius_postProc2 is the second part of the function 
    sirius_postProc defined in R part of the workflow. This function
    re-checks the Suspect list, if present or given as a parameter, 
    whether the candidates have a high similarity with compounds in
    Suspect List. It also calculates the Maximum Common Substructure
    (MCSS)

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
    SIRIUS json files
    
    sl (bool): True if a suspct list is to be used
    
    slistcsv (list): This is the csv file that contains a column of 
    "SMILES". Additionally this file can contain other information 
    about the compounds, but for this function, column of "SMILES", 
    named as "SMILES" is necessary.

    Returns:
    csv: a result file with additional columns such as those for suspect
    list if one is used. It also adds columns on MCSS., named as 
    "input_dir/ResultFileName/insilico/SiriusResults.csv"
    
    
    Usage:
    sirius_postProc2(input_dir = "/user/project/", 
    input_table = "/user/project/suspectlist.csv", sl = True, slistcsv)


    """
    
    # Describe the heavy atoms to be considered for MCSS
    heavy_atoms = ['C', 'N', 'P', 'O', 'S']
    
    input_table = pd.read_csv(input_tablecsv)
    
    for m, row in input_table.iterrows():
        
        # Read the file result_dir/insilico/MS1DATAsirius.csv. 
        # This file has been produced in R workflow and contains 
        # SIRIUS results.

        file1 = pd.read_csv(input_dir + (input_table['ResultFileNames'][m] + '/insilico/MS1DATAsirius.csv').replace("./", ""))
        
        for i, row in file1.iterrows():
            
            # if the entry has SMILES extracted for MCSS calculation
            if not isNaN(file1['SMILESforMCSS'][i]):
                
                # split the SMILES using |
                top_smiles = file1['SMILESforMCSS'][i].split("|")
                
                # if there are more than 1 smiles in the top smiles, 
                if len(top_smiles) > 1:
                    
                    # define empty mol list to add mol objects from the top Smiles
                    mol = []
                    
                    # is sl = True
                    if sl:
                        
                        # Add columns 
                        file1['Top_can_SL'] = np.nan # top candidate among the top 5 candidates, according to similarity with a compound in suspect list
                        file1['tanimotoSLvsCAN'] = np.nan # tanimoto score
                        file1['SL_comp'] = np.nan # Smiles of the suspect listr compund with  high similairity with the one of the top 5 candidates
                        
                        # read the suspect list
                        slist = pd.read_csv(slistcsv)
                        
                        for j in top_smiles:
                            
                            for k, row in slist.iterrows():
                                
                                # calculate the tanimoto between SMILES in top smiles and the suspect list
                                SSms = [Chem.MolFromSmiles(j), Chem.MolFromSmiles(slist['SMILES'][k])]
                                SSfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in SSms]
                                SStn = DataStructs.FingerprintSimilarity(SSfps[0],SSfps[1])
                                
                                # if there is high similarity, keep that entity in the following columns to consider later
                                if SStn >= 0.8:
                                    file1.loc[i,'Top_can_SL'] = j
                                    file1.loc[i,'tanimotoSLvsCAN'] = SStn
                                    file1.loc[i,'SL_comp'] = slist['SMILES'][k]
                            # calculate the mol object from each smiles
                            mol_object = Chem.MolFromSmiles(j)
                            # add all these mol objects to mol
                            mol.append(mol_object)
                    # list of mol used to calaculate the MCSS
                    res = rdFMCS.FindMCS(mol)
                    sm_res = Chem.MolToSmiles(Chem.MolFromSmarts(res.smartsString))
                    
                    # Check if the MCSS has one of the heavy atoms and whether they are
                    # more than 3
                    elem = [ele for ele in heavy_atoms if(ele in sm_res)]
                    if elem and len(sm_res)>=3:
                        file1.loc[i, 'MCSSstring'] = res.smartsString
                        file1.loc[i, 'MCSS_SMILES'] = Chem.MolToSmiles(Chem.MolFromSmarts(res.smartsString))

        file1.to_csv(input_dir + (input_table['ResultFileNames'][m] + '/insilico/SiriusResults.csv').replace("./", ""))
        
        
        
sirius_postProc2(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

