#!/usr/bin/env python
# coding: utf-8

# ## SIRIUS_Metfrag_SList

# In[1]:


import pandas as pd

import pubchempy as pcp
import numpy as np
def isNaN(string):
    return string != string
import os
import glob
import re
from pybatchclassyfire import *
import csv 
import time
import json
from pandas import json_normalize


#import openpyxl
import statistics
import sys


# In[2]:


from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem import PandasTools


# In[3]:


# make sure your Smiles entries in the suspect list csv are in a column named "SMILES"
def slist_metfrag(input_dir, slist_csv):
    """slist_metfrag is used to create a txt file that contains a list of 
    InChIKeys. This list is later used by MetFrag to use these compounds 
    as a Suspect List.

    Parameters:
    input_dir (str): This is the input directory where all the .mzML 
    files and their respective result directories are stored. For this 
    function this directory must contain a csv file that has a column 
    named "SMILES".
    
    slist_csv (str): This is the csv file that contains a column of 
    "SMILES". Additionally this file can contain other information 
    about the compounds, but for this function, column of "SMILES", 
    named as "SMILES" is necessary.

    Returns:
    list: list of InChIKeys
    txt: a txt file of list of InChIKeys, is stored in input_dir
    
    Usage:
    slist_metfrag(input_dir = "/user/project/", slist_csv = 
    "suspectlist.csv")
    
    """
    sl = pd.read_csv(slist_csv)
    sl_mtfrag= []
    for i, rows in sl.iterrows():
        mols = Chem.MolFromSmiles(sl['SMILES'][i])
        sl.loc[i, 'InChIKey'] = Chem.inchi.MolToInchiKey(mols)
        sl_mtfrag.append(sl['InChIKey'][i])
    return(sl_mtfrag)
    with open((input_dir + 'SL_metfrag.txt'), 'w') as filehandle:
        for listitem in sl_mtfrag:
            filehandle.write('%s\n' % listitem)


# In[4]:


#print(slist_metfrag.__doc__)


# In[ ]:





# In[5]:


def slist_sirius(input_dir, slist_csv, substring = None):
    
    """slist_sirius is used to create a tsv file that contains a list of 
    SMILES. The function also runs the sirius command custom db to create
    fingerprints for each SMILES in a folder that we by default name as
    SL_Frag/. This fingerprints folder is later used by SIRIUS to use 
    these compounds as a another small list of compounds to match against
    the input spectra fingerprints.
    Since SIRIUS doesn't take disconnected structure, Multiply charged, 
    Incorrect syntax, wild card(*) in smiles; this function removes all
    such SMILES from the Suspect List.

    Parameters:
    input_dir (str): This is the input directory where all the .mzML 
    files and their respective result directories are stored. For this 
    function this directory must contain a csv file that has a column 
    named "SMILES".
    
    slist_csv (str): This is the csv file that contains a column of 
    "SMILES". Additionally this file can contain other information 
    about the compounds, but for this function, column of "SMILES", 
    named as "SMILES" is necessary.
    
    substring (list): provide a list of strings of SMILES that 
    shouldn't be considered, provide a list even if there is one string
    that shouldnt be considered. e.g: "[Fe+2]". 

    Returns:
    tsv: a tsv file of list of SMILES, named as SL_Sirius.tsv, is stored 
    in input_dir
    directory: a directory with compound fragmentations will be created 
    in a folder named SL_Frag/ within the same input_dir
    
    
    Usage:
    slist_sirius("/user/project/", "suspectlist.csv", 
    substring = None)

    """
    
    sl = pd.read_csv(slist_csv)
    
    # define function to neutralize the charged SMILES
    def neutralize_atoms(mol):
        
        pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
        at_matches = mol.GetSubstructMatches(pattern)
        at_matches_list = [y[0] for y in at_matches]
        if len(at_matches_list) > 0:
            for at_idx in at_matches_list:
                atom = mol.GetAtomWithIdx(at_idx)
                chg = atom.GetFormalCharge()
                hcount = atom.GetTotalNumHs()
                atom.SetFormalCharge(0)
                atom.SetNumExplicitHs(hcount - chg)
                atom.UpdatePropertyCache()
        return mol
    

    for i, row in sl.iterrows():
        # remove SMILES with wild card
        if "*" in sl["SMILES"][i]:
            sl = sl.drop(labels = i, axis = 0) 
    for i, row in sl.iterrows():
        # remove SMILES with any string present in the substring
        if substring:
            if bool([ele for ele in substring if(ele in sl["SMILES"][i])]):
                sl = sl.drop(labels = i, axis = 0)
    for i, row in sl.iterrows():
        if "." in sl["SMILES"][i]:
            sl.loc[i, "SMILES"] = sl["SMILES"][i].split('.')[0]
    # Neutralize the charged SMILES
    for i, row in sl.iterrows():
        if "+" in sl["SMILES"][i] or "-" in sl["SMILES"][i]:
            mol = Chem.MolFromSmiles(sl["SMILES"][i])
            neutralize_atoms(mol)
            sl.loc[i, "SMILES"] = Chem.MolToSmiles(mol)
            
            # Remove multiple charged SMILES
            if "+" in sl["SMILES"][i] or "-" in sl["SMILES"][i]:
                pos = sl["SMILES"][i].count('+')
                neg = sl["SMILES"][i].count('-')
                charge = pos + neg 
                if charge > 1:
                    sl = sl.drop(labels = i, axis = 0) 
                    
    slsirius = pd.DataFrame({'smiles':sl["SMILES"]})
    slsirius.to_csv(input_dir+ "SL_Sirius.tsv", sep = "\t", header = False, index = False)
    os.system("sirius --input " + input_dir + "SL_Sirius.tsv custom-db --name=SL_Frag --output "+ input_dir)
# Usage:
# slist_sirius("/Users/user/project/SuspectList/", "SUSPECT_LIST.csv", substring = None)


# In[6]:


#print(slist_sirius.__doc__)


# ### SIRIUS post processing

# ### SIRIUS Result Post Processing

# In[7]:


def sirius_postProc2(input_dir, input_table, slistcsv ,sl = True):
    
    
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
    
    input_table (str): This is the table in csv format (defined in R), 
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
        file1.to_csv(input_table['ResultFileNames'][m] + '/insilico/SiriusResults.csv')


# In[8]:


#print(sirius_postProc2.__doc__)


# ### MetFrag Result Post Processing

# In[39]:


def metfrag_postproc(input_dir, input_table, slistcsv ,sl = True):
    
    
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
    
    input_table (str): This is the table in csv format (defined in R), 
    which stores a csv table containing columns "mzml_files", which 
    contains liat of all input files with their relative paths, second
    column is "ResultFileName" which is a list of the corresponding
    result relative directories to each mzml files. Lastly, "file_id", 
    contains a file directory. This table will be used to read the 
    MetFrag csv files
    
    sl (bool): True if a suspct list is to be used
    
    slistcsv (list): This is the csv file that contains a column of 
    "SMILES". Additionally this file can contain other information 
    about the compounds, but for this function, column of "SMILES", 
    named as "SMILES" is necessary.

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
    
    
    for m, row in input_table.iterrows():
    
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
                    KEGG_file = KEGG_file.drop(KEGG_file[KEGG_file.Score < 0.75].index)
                
                    # add the relavnt information to the original MS1DATA csv
                    #file1.loc[i, 'KG_ID'] = KEGG_file.loc[0, 'Identifier']
                    file1.loc[i, 'KG_Name'] = KEGG_file.loc[0, 'CompoundName']
                    file1.loc[i, 'KG_Formula'] = KEGG_file.loc[0, 'MolecularFormula']
                    file1.loc[i, 'KG_expPeaks'] = KEGG_file.loc[0, 'NoExplPeaks']
                    file1.loc[i, 'KG_SMILES'] = Chem.MolToSmiles(Chem.MolFromInchi(KEGG_file["InChI"][0]))
                    file1.loc[i, 'KG_file'] = KEGG
                
                    #create empty list of KEGG top smiles
                    Kegg_smiles = []
                
                    # extract only the InChI of the top 5
                    for j in KEGG_file["InChI"][0:5].tolist():
                        # convert the InChI to SMILES
                        mol = Chem.MolToSmiles(Chem.MolFromInchi(j))
                        if sl:
                            # read the suspect list
                            slist = pd.read_csv(slistcsv)
                            
                            # Add columns 
                            file1['KG_Top_can_SL'] = np.nan # top candidate among the top 5 candidates, according to similarity with a compound in suspect list
                            file1['KG_tanimotoSLvsCAN'] = np.nan # tanimoto score
                            file1['KG_SL_comp'] = np.nan # Smiles of the suspect listr compund with  high similairity with the one of the top 5 candidates
                            
                            # for each smiles in suspect list
                            for k, row in slist.iterrows():
                                # Calculate the tanimoto score
                                SSms = [Chem.MolFromSmiles(mol), Chem.MolFromSmiles(slist['SMILES'][k])]
                                SSfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in SSms]
                                SStn = DataStructs.FingerprintSimilarity(SSfps[0],SSfps[1])
                                if SStn >= 0.8:
                                    file1.loc[i, 'KG_Top_can_SL'] = j
                                    file1.loc[i, 'KG_tanimotoSLvsCAN'] = SStn
                                    file1.loc[i, 'KG_SL_comp'] = slist['SMILES'][k]
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
                    
                    # take the ones with more than 0.75 score
                    PubChem_file = PubChem_file.drop(PubChem_file[PubChem_file.Score < 0.75].index)
                    
                    # add the relavnt information to the original MS1DATA csv
                    file1.loc[i, 'PC_ID'] = PubChem_file.loc[0, 'Identifier']
                    file1.loc[i, 'PC_Name'] = PubChem_file.loc[0, 'IUPACName']
                    file1.loc[i, 'PC_Formula'] = PubChem_file.loc[0, 'MolecularFormula']
                    file1.loc[i, 'PC_expPeaks'] = PubChem_file.loc[0, 'NoExplPeaks']
                    file1.loc[i, 'PC_SMILES'] = PubChem_file["SMILES"][0]
                    file1.loc[i, 'PC_file'] = PubChem
                    
                    # empty object
                    Pubchem_smiles = []
                    
                    # extract only the SMILES of the top 5
                    for j in PubChem_file["SMILES"][0:5].tolist():
                        
                        # if sl = True
                        if sl:
                            
                            # read the suspect list
                            slist = pd.read_csv(slistcsv)
                            
                            # Add columns 
                            file1['PC_Top_can_SL'] = np.nan # top candidate among the top 5 candidates, according to similarity with a compound in suspect list
                            file1['PC_tanimotoSLvsCAN'] = np.nan # tanimoto score
                            file1['PC_SL_comp'] = np.nan # Smiles of the suspect listr compund with  high similairity with the one of the top 5 candidates
                            # calculate tanimoto
                            for n, row in slist.iterrows():
                                
                                SSms = [Chem.MolFromSmiles(j), Chem.MolFromSmiles(slist['SMILES'][n])]
                                SSfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in SSms]
                                SStn2 = DataStructs.FingerprintSimilarity(SSfps[0],SSfps[1])
                                
                                if SStn2 >= 0.8:
                                    file1.loc[i, 'PC_Top_can_SL'] = j
                                    file1.loc[i, 'PC_tanimotoSLvsCAN'] = SStn2
                                    file1.loc[i, 'PC_SL_comp'] = slist['SMILES'][n]
                                    
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


# In[10]:


#print(metfrag_postproc.__doc__)


# ### COMBINE IN SILICO -All files with SIRIUS results separate and with MetFragresults separate

# In[11]:


def combine_insilico(input_dir, input_table, Source = "SIRIUS"):
    
    """combine_insilico function combines the Sirius results from all
    result directories for each input mzml file. It does same for 
    Metfrag.

    Parameters:
    input_dir (str): This is the input directory where all the .mzML 
    files and their respective result directories are stored.
    
    input_table (str): This is the table in csv format (defined in R), 
    which stores a csv table containing columns "mzml_files", which 
    contains liat of all input files with their relative paths, second
    column is "ResultFileName" which is a list of the corresponding
    result relative directories to each mzml files. Lastly, "file_id", 
    contains a file directory. This table will be used to read the 
    Sirius and MetFrag result csv files
    
    Source (str): either "SIRIUS" or "MetFrag"

    Returns:
    
    dataframe: of combined SIRIUS/MetFrag results
    
    csv: stores the dataframe in a csv, named as 
    "input_dir/ResultFileName/MetabolomicsResults/SIRIUS_combined.csv" 
    OR/AND 
    "input_dir/ResultFileName/MetabolomicsResults/MetFrag_combined.csv"
    
    
    Usage:
    combine_insilico(input_dir = "/user/project/", 
    input_table = "/user/project/suspectlist.csv", Source = "SIRIUS")


    """
    # create a new directory to store all results /MetabolomicsResults/
    path = os.path.join(input_dir, "MetabolomicsResults")
    if not os.path.isdir(path):
        os.mkdir(path)    
    # if Sirius results are to be combined
    if Source == "SIRIUS":
        
        # store all files paths here
        all_files = []
        for n, row in input_table.iterrows():
            all_files.append(input_dir + input_table['ResultFileNames'][n].replace("./", "") + '/insilico/SiriusResults.csv')
        
        # store all dataframes of the results here
        li = []
    
        for filename in all_files:
            df = pd.read_csv(filename, index_col=None, header=0)
            df["ResultFileNames"] = filename
            li.append(df)
            
        # join all resulst dataframe
        frame = pd.concat(li, axis=0, ignore_index=True)
        

        for i, row in frame.iterrows():
            if frame["FormulaRank"][i] == 1.0:
                sep = 'json/'
                strpd = frame["dir"][i].split(sep, 1)[0] +"json/canopus_summary.tsv"
                if os.path.isfile(strpd):

                    canopus = pd.read_csv(strpd, sep='\t')
                    if len(canopus) > 0:
                        frame.loc[i, 'most_specific_class'] = canopus["most specific class"][0]
                        frame.loc[i, 'level _5'] = canopus["level 5"][0]
                        frame.loc[i, 'subclass'] = canopus["subclass"][0]
                        frame.loc[i, 'class'] = canopus["class"][0]
                        frame.loc[i, 'superclass'] = canopus["superclass"][0]
                        frame.loc[i, 'all_classifications'] = canopus["all classifications"][0]
                        frame.loc[i, 'Classification_Source'] = 'CANOPUS'
        frame.to_csv(input_dir + '/MetabolomicsResults/SIRIUS_combined.csv')
        return(frame)
    
    # if MetFrag results are to be combined
    elif Source == "MetFrag":
        
        # store all files paths here
        all_files = []
        for m, row in input_table.iterrows():
            all_files.append(input_dir + input_table['ResultFileNames'][m].replace("./", "") + '/insilico/MetFragResults.csv')
        li = []

        for filename in all_files:
            df = pd.read_csv(filename, index_col=None, header=0)
            df["result_dir"] = filename
            li.append(df)

        frame = pd.concat(li, axis=0, ignore_index=True)
        frame.to_csv(input_dir+'MetabolomicsResults/MetFrag_combined.csv')
        return(frame)


# In[12]:


#print(combine_insilico.__doc__)


# # Spectral DB Dereplication Results Post Processing

# ### GNPS, MassBank and HMDB Results post processing

# In[3]:


def spec_postproc(input_dir, Source = "all"):
    
    """spec_postproc function processes the resulst from dereplication 
    using different spectral DBs. 

    Parameters:
    input_dir (str): This is the input directory where all the .mzML 
    files and their respective result directories are stored.
    
    Source (str): either "mbank" or "hmdb" or "gnps", or "all"

    Returns:
    
    dataframe: of the paths of the processed DB results
    
    
    Usage:
    spec_postproc(input_dir = "/user/project/", Source = "all")

    """
    
    # empty lists of csv files paths for each database
    GNPScsvfiles = []
    HMDBcsvfiles = []
    MassBankcsvfiles = []
    
    #list all files and directories
    for entry in os.listdir(input_dir):
        if os.path.isdir(os.path.join(input_dir, entry)):
            
            # enter the directory with /spectral_dereplication/ results
            sub_dir = input_dir + entry + '/spectral_dereplication'
            if os.path.exists(sub_dir):
                files = (glob.glob(sub_dir+'/*.csv'))

                for f in files:
                    if 'gnps.' in f: 
                        GNPScsvfiles.append(f)
                    if 'hmdb.' in f: 
                        HMDBcsvfiles.append(f)
                    if 'mbank.' in f: 
                        MassBankcsvfiles.append(f)
                            
    
    if Source == "hmdb" or "all":
        #download SDF structures
        #os.system("wget https://hmdb.ca/system/downloads/current/structures.zip")
        #os.system("unzip "+ input_dir + "structures.zip")
            
        # Load the sdf
        dframe = PandasTools.LoadSDF((input_dir+"structures.sdf"),
                                     idName='HMDB_ID',smilesName='SMILES',
                                     molColName='Molecule', includeFingerprints=False)
        
        #### read sdf file from HMDB to collect names and smiles ####
    
        #HMDB CSV Result file pre_processing
        
        #open another csv path holding empty list, which will be filled 
        #with post processed csv results
        HMDBcsvfiles2 = []
        
        for k in HMDBcsvfiles:
            
            # read the csv files
            hmdb_df = pd.read_csv(k)
            
            # merge on basis of id, frame and hmdb result files
            SmilesHM = pd.merge(hmdb_df, dframe, left_on=hmdb_df.HMDBcompoundID, right_on=dframe.DATABASE_ID)
            
            
            for i, row in hmdb_df.iterrows():
                
                for j, row in SmilesHM.iterrows():
                    
                    # where index for both match, add the name and SMILES
                    if hmdb_df['id_X'][i]== SmilesHM['id_X'][j]:
                        hmdb_df.loc[i, 'HMDBSMILES'] = SmilesHM['SMILES'][j]#add SMILES
                        hmdb_df.loc[i, 'HMDBcompound_name'] = SmilesHM["GENERIC_NAME"][j]#add name
                        hmdb_df.loc[i, 'HMDBformula'] = SmilesHM["FORMULA"][j]#add formula
                
            csvname = (os.path.splitext(k)[0])+"proc"+".csv" # name for writing it in a new file
            hmdb_df.to_csv(csvname) #write
            HMDBcsvfiles2.append(csvname)# add to a list
        
    #MassBank CSV Result file pre_processing
    
    if Source == "mbank" or "all":
        
        #open another csv path holding empty list, which will be filled 
        #with post processed csv results
        MassBankcsvfiles2 = []
        
        for l in MassBankcsvfiles:
            
            # read mbank csv file
            mbank_df = pd.read_csv(l)
            
            for i, row in mbank_df.iterrows():
                
                inchiK = str(mbank_df["MBinchiKEY"][i])
                
                #extract inchikeys
                y = pcp.get_compounds(inchiK, 'inchikey')#compound based on inchikey
                
                for compound in y:
                    
                    #add smiles
                    smles = compound.isomeric_smiles   
                    mbank_df.loc[i, 'MBSMILES'] = smles
                    
            csvname = (os.path.splitext(l)[0])+"proc"+".csv"
            mbank_df.to_csv(csvname)
            MassBankcsvfiles2.append(csvname)
    
    # GNPS CSV Result file pre_processing
    if Source == "gnps" or "all":
        
        #currently only these subsets are removed from the names from GNPS
        matches = ["M+","[M", "M-", "2M", "M*" "20.0", "50.0", "30.0", "40.0", "60.0", "70.0", "eV", "Massbank"
               , "Spectral", "Match", "to", "from", "NIST14", "MoNA", '[IIN-based:',  '[IIN-based', 'on:', 'CCMSLIB00003136269]']
        
        #open another csv path holding empty list, which will be filled 
        #with post processed csv results
        GNPScsvfiles2 = []
        
        for l in GNPScsvfiles:
            
            gnps_df = pd.read_csv(l)
    
            for i, row in gnps_df.iterrows():
            
                # if compound name is present
                if not isNaN(gnps_df['GNPScompound_name'][i]):
                    
                    # split if there is a gap in the names
                    string_chng = (gnps_df['GNPScompound_name'][i].split(" "))
                    
                    # create an empty list
                    newstr = []
                    
                    # for each part of the string in the names
                    chng = []
                    
                    for j in range(len(string_chng)):
                        
                        # check if the substrings re present in the matches and no - is present
                        if not any(x in string_chng[j] for x in matches) and not '-' == string_chng[j]:
                            
                            # IF | and ! not in the substring
                            if '|' not in string_chng[j] or '!' not in string_chng[j]:
                                newstr.append(string_chng[j])
                                
                            # if | present in the substring   
                            elif '|' in string_chng[j]:
                                
                                #split the string
                                jlen = string_chng[j].split("|")
                                #how many substrings are left now
                                lst = len(jlen)-1
                                #append this to chng
                                chng.append(jlen[lst])
                                break
                                
                    # now append chng to newstr            
                    chng.append(' '.join(newstr))
                    #save this as the correct name
                    gnps_df.loc[i, "corr_names"] = chng[0]
                    
                    if not isNaN(gnps_df['GNPSSMILES'][i]):
                        if chng == '':
                            break
                        elif gnps_df['GNPSSMILES'][i].isalpha():
                            s = pcp.get_compounds(chng[0], 'name')
                            if s:
                                for comp in s:
                                    gnps_df["GNPSSMILES"][i] = comp.isomeric_smiles
                            else:
                                gnps_df["GNPSSMILES"][i] = ''
                else:
                    gnps_df["GNPSSMILES"][i] = ''
                    
            for i, row in gnps_df.iterrows():
                if not isNaN(gnps_df['GNPSSMILES'][i]):
                    try:
                        sx = pcp.get_compounds(gnps_df['GNPSSMILES'][i], 'smiles')
                        if sx:
                            sx = str(sx)
                            comp = pcp.Compound.from_cid([int(x) for x in re.findall(r'\b\d+\b', sx)])
                            gnps_df.loc[i, 'GNPSformula'] = comp.molecular_formula
                    except:
                        gnps_df.loc[i, 'GNPSformula'] = ''
                        
            csvname = (os.path.splitext(l)[0])+"_with_cor_names"+".csv"
            gnps_df.to_csv(csvname)
            GNPScsvfiles2.append(csvname)
    

    if Source == "all":
        
        dict1 = {'GNPSr': GNPScsvfiles2, 'HMDBr': HMDBcsvfiles2, 'MBr': MassBankcsvfiles2} 
        df = pd.DataFrame(dict1)

        return(df)
    


# In[14]:


#print(spec_postproc.__doc__)


# ### Combine_all Spectral DBs for one file

# In[4]:


def combine_specdb(input_dir):
    
    """combine_specdb function combines all results from different
    spectral dbs. Can only be used if more than one db is used 

    Parameters:
    input_dir (str): This is the input directory where all the .mzML 
    files and their respective result directories are stored.

    Returns:
    dataframe: of the paths of the merged results
    
    
    Usage:
    combine_specdb(input_dir)

    """
    
    # empty lists of csv files paths for each database
    GNPScsvfiles2 = []
    HMDBcsvfiles2 = []
    MassBankcsvfiles2 = []
    
    #list all files and directories
    for entry in os.listdir(input_dir):
        if os.path.isdir(os.path.join(input_dir, entry)):
            
            # enter the directory with /spectral_dereplication/ results
            sub_dir = input_dir + entry + '/spectral_dereplication'
            if os.path.exists(sub_dir):
                files = (glob.glob(sub_dir+'/*.csv'))

                for f in files:
                    if 'gnps_with_' in f: 
                        GNPScsvfiles2.append(f)
                    if 'hmdbproc.' in f: 
                        HMDBcsvfiles2.append(f)
                    if 'mbankproc.' in f: 
                        MassBankcsvfiles2.append(f)
    
    dict1 = {'GNPSr': GNPScsvfiles2, 'HMDBr': HMDBcsvfiles2, 'MBr': MassBankcsvfiles2} 
    df = pd.DataFrame(dict1)
    
    Merged_Result_df = []
    
    
    for i, row in df.iterrows():
        
        
        CSVfileG = pd.read_csv(df["GNPSr"][i])
        CSVfileH = pd.read_csv(df["HMDBr"][i])
        CSVfileM = pd.read_csv(df["MBr"][i])
        
        # merge on the basis of Idx
        MergedRE = CSVfileG.merge(CSVfileH,on='id_X').merge(CSVfileM,on='id_X')

        csvname = (df["GNPSr"][i]).replace("gnps_with_cor_names", "mergedR")
        MergedRE.to_csv(csvname)
        Merged_Result_df.append(csvname)
        
    return(Merged_Result_df)


# In[16]:


#print(combine_specdb.__doc__)


# ### Combine all files for spectral db dereplication

# In[17]:


def combine_allspec(input_dir):
    
    """combine_allspec function combines all results from different
    spectral dbs. Can only be used if more than one db is used 

    Parameters:
    input_dir (str): This is the input directory where all the .mzML 
    files and their respective result directories are stored.
    df (dataframe): dataframe from combine_specdb
    
    Returns:
    dataframe: of the paths of the merged results from all files
    
    Usage:
    combine_allspec(input_dir = "usr/project/", comb_df)

    """
    
    Mergedcsvfiles = []
    
    #list all files and directories
    for entry in os.listdir(input_dir):
        if os.path.isdir(os.path.join(input_dir, entry)):
            
            # enter the directory with /spectral_dereplication/ results
            sub_dir = input_dir + entry + '/spectral_dereplication'
            if os.path.exists(sub_dir):
                files = (glob.glob(sub_dir+'/*.csv'))

                for f in files:
                    if 'mergedR.csv' in f: 
                        Mergedcsvfiles.append(f)
    
    combined_csv = pd.concat([pd.read_csv(l) for l in Mergedcsvfiles], ignore_index=True)
    
    for i, row in combined_csv.iterrows():
        if combined_csv['GNPSSMILES'][i] == ' ' or isNaN(combined_csv['GNPSSMILES'][i]):
            combined_csv['GNPSSMILES'][i] = ''
            
    #for i, row in combined_csv.iterrows():
        #if not isNaN(combined_csv['MBinchiKEY'][i]):
            #try:
                #y = pcp.get_compounds(combined_csv['MBinchiKEY'][i], 'inchikey')
                #if len(y)>1:
                    #combined_csv['MBSMILES'][i] = y[0].isomeric_smiles
            #except:
                #pass
            
    combined_csv.to_csv(input_dir + 'MetabolomicsResults/SD_post_processed_combined_results.csv')
    return(combined_csv)


# In[18]:


#print(combine_allspec.__doc__)


# ### Scoring Scheme for Spectral DB Dereplication

# In[ ]:





# In[19]:


def scoring_spec(input_dir, combined):
    
    """scoring_spec extracts the candidates with high scores from
    the results from combine_allspec function 

    Parameters:
    input_dir (str): This is the input directory where all the .mzML 
    files and their respective result directories are stored.
    combined (dataframe): dataframe from combine_allspec
    
    Returns:
    dataframe: of the all features and their results
    csv: CSV reuslt file named MetabolomicsResults/combinedSpecDB.csv
    which contains all the features and their Spec DB annotations
    
    Usage:
    scoring_spec(input_dir = "usr/project/", combined)

    """
    
    # the scoring highly depends on the following information:
    # similarity scores should be higher than 0.75
    # intScore >=0.50
    # mzScore >= 0.50
    # ratio of the matchingpeaks by the totalpeaks in the query >= 0.50
    
    combined = pd.read_csv(combined)
    
    def HMDB_Scoring(db, i):
        if db['HMDBmax_similarity'][i] >= 0.75 and db['HMDBintScore'][i] >= 0.50 and db['HMDBmzScore'][i] >= 0.50 and db['HQMatchingPeaks'][i]/db['hQueryTotalPeaks'][i] >= 0.50:
            return True
        else:
            return False
    
    
    def GNPS_Scoring(db, i):
        if db['GNPSmax_similarity'][i] >= 0.75 and db['GNPSintScore'][i] >= 0.50 and db['GNPSmzScore'][i] >= 0.50 and db['GQMatchingPeaks'][i]/db['gQueryTotalPeaks'][i] >= 0.50:
            return True
        else:
            return False
    
    
    def MB_Scoring(db, i):
        if db['MBmax_similarity'][i] >= 0.75 and db['MBintScore'][i] >= 0.50 and db['MBmzScore'][i] >= 0.50 and db['MQMatchingPeaks'][i]/db['mQueryTotalPeaks'][i] >= 0.50:
            return True
        else:
            return False
    
    
    for i, row in combined.iterrows():
        
        # if all DBs show good candidates accorindg to the scoring
        if HMDB_Scoring(combined, i) and GNPS_Scoring(combined, i) and MB_Scoring(combined, i) and not isNaN(combined['GNPSSMILES'][i]) and not isNaN(combined['MBSMILES'][i]) and not isNaN(combined['HMDBSMILES'][i]):
            
            # calulate the tanimoto similarity between the candidates from three DBs
            
            # hmdb and gnps
            HGms = [Chem.MolFromSmiles(combined['HMDBSMILES'][i]), Chem.MolFromSmiles(combined['GNPSSMILES'][i])]
            HGfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in HGms]
            HGtn = DataStructs.FingerprintSimilarity(HGfps[0],HGfps[1])
            
            # gnps and mbank
            GMms = [Chem.MolFromSmiles(combined['GNPSSMILES'][i]), Chem.MolFromSmiles(combined['MBSMILES'][i])]
            GMfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in GMms]
            GMtn = DataStructs.FingerprintSimilarity(GMfps[0],GMfps[1])
            
            # mbank and hmdb
            HMms = [Chem.MolFromSmiles(combined['HMDBSMILES'][i]), Chem.MolFromSmiles(combined['MBSMILES'][i])]
            HMfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in HMms]
            HMtn = DataStructs.FingerprintSimilarity(HMfps[0],HMfps[1])
            
            # add the following columns
            combined.loc[i, 'annotation'] = 'HMDB, GNPS, MassBank'
            combined.loc[i, 'tanimotoHG'] = HGtn
            combined.loc[i, 'tanimotoGM'] = GMtn
            combined.loc[i, 'tanimotoHM'] = HMtn
            combined.loc[i, 'occurence'] = 3
        
        # if HMDB and GNPS show good candidates accorindg to the scoring
        if HMDB_Scoring(combined, i) and GNPS_Scoring(combined, i) and not MB_Scoring(combined, i) and not isNaN(combined['GNPSSMILES'][i]) and not isNaN(combined['HMDBSMILES'][i]):
            HGms = [Chem.MolFromSmiles(combined['HMDBSMILES'][i]), Chem.MolFromSmiles(combined['GNPSSMILES'][i])]
            HGfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in HGms]
            HGtn = DataStructs.FingerprintSimilarity(HGfps[0],HGfps[1])
        
            combined.loc[i, 'annotation'] = 'HMDB, GNPS'
            combined.loc[i, 'tanimotoHG'] = HGtn
            combined.loc[i, 'tanimotoGM'] = np.nan
            combined.loc[i, 'tanimotoHM'] = np.nan
            combined.loc[i, 'occurence'] = 2
        
        # if MassBank and GNPS show good candidates accorindg to the scoring
        if not HMDB_Scoring(combined, i) and GNPS_Scoring(combined, i) and MB_Scoring(combined, i) and not isNaN(combined['MBSMILES'][i]) and not isNaN(combined['GNPSSMILES'][i]):
            GMms = [Chem.MolFromSmiles(combined['GNPSSMILES'][i]), Chem.MolFromSmiles(combined['MBSMILES'][i])]
            GMfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in GMms]
            GMtn = DataStructs.FingerprintSimilarity(GMfps[0],GMfps[1])
        
            combined.loc[i, 'annotation'] = 'GNPS, MassBank'
            combined.loc[i, 'tanimotoHG'] = np.nan
            combined.loc[i, 'tanimotoGM'] = GMtn
            combined.loc[i, 'tanimotoHM'] = np.nan
            combined.loc[i, 'occurence'] = 2
        
        # if MassBank and HMDB show good candidates accorindg to the scoring
        if HMDB_Scoring(combined, i) and not GNPS_Scoring(combined, i) and MB_Scoring(combined, i) and not isNaN(combined['MBSMILES'][i]) and not isNaN(combined['HMDBSMILES'][i]):
            HMms = [Chem.MolFromSmiles(combined['HMDBSMILES'][i]), Chem.MolFromSmiles(combined['MBSMILES'][i])]
            HMfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in HMms]
            HMtn = DataStructs.FingerprintSimilarity(HMfps[0],HMfps[1])
        
            combined.loc[i, 'annotation'] = 'HMDB, MassBank'
            combined.loc[i, 'tanimotoHG'] = np.nan
            combined.loc[i, 'tanimotoGM'] = np.nan
            combined.loc[i, 'tanimotoHM'] = HMtn
            combined.loc[i, 'occurence'] = 2
        
        # only HMDB
        if HMDB_Scoring(combined, i) and not GNPS_Scoring(combined, i) and not MB_Scoring(combined, i):
        
            combined.loc[i, 'annotation'] = 'HMDB'
            combined.loc[i, 'tanimotoHG'] = np.nan
            combined.loc[i, 'tanimotoGM'] = np.nan
            combined.loc[i, 'tanimotoHM'] = np.nan
            combined.loc[i, 'occurence'] = 1
        
        # only GNPS
        if not HMDB_Scoring(combined, i) and GNPS_Scoring(combined, i) and not MB_Scoring(combined, i):
        
            combined.loc[i, 'annotation'] = 'GNPS'
            combined.loc[i, 'tanimotoHG'] = np.nan
            combined.loc[i, 'tanimotoGM'] = np.nan
            combined.loc[i, 'tanimotoHM'] = np.nan
            combined.loc[i, 'occurence'] = 1
        
        # only MassBank
        if not HMDB_Scoring(combined, i) and not GNPS_Scoring(combined, i) and MB_Scoring(combined, i):
        
            combined.loc[i, 'annotation'] = 'MassBank'
            combined.loc[i, 'tanimotoHG'] = np.nan
            combined.loc[i, 'tanimotoGM'] = np.nan
            combined.loc[i, 'tanimotoHM'] = np.nan
            combined.loc[i, 'occurence'] = 1
        
        # none
        if not HMDB_Scoring(combined, i) and not GNPS_Scoring(combined, i) and not MB_Scoring(combined, i):
            combined.loc[i, 'annotation'] = 'none'
            combined.loc[i, 'tanimotoHG'] = np.nan
            combined.loc[i, 'tanimotoGM'] = np.nan
            combined.loc[i, 'tanimotoHM'] = np.nan
            combined.loc[i, 'occurence'] = 0
    
    combined.to_csv(input_dir + "MetabolomicsResults/combinedSpecDB.csv")
    return(combined)


# In[20]:


#print(scoring_spec.__doc__)


# ### Suspect List Screening

# In[21]:


def suspectListScreening(input_dir, slistcsv, SpectralDB_Results):
    
    """suspectListScreening runs tanoimoto similarity score to between
    compounds from the results from spectral DBs and suspect list

    Parameters:
    input_dir (str): This is the input directory where all the .mzML 
    files and their respective result directories are stored.
    slistcsv (str): path to suspect list
    SpectralDB_Results (dataframe): dataframe from scoring_spec
    
    Returns:
    dataframe: all features and specDB reults and suspect list screening 
    results
    csv: CSV reuslt file named MetabolomicsResults/SpecDBvsSL.csv
    which contains all the features and their Spec DB annotations
    and suspect list occurences if any
    
    Usage:
    suspectListScreening(input_dir = "usr/project/",
    slistcsv = "usr/project/suspect_list.csv", 
    SpectralDB_Results)

    """
    
    SpectralDB_Results = pd.read_csv(SpectralDB_Results)
    
    # add columns to the resulst from scoring_spec
    # these columns are for high similiarity canidtes between the databases and suspect list
    SpectralDB_Results['HLsmiles'] = np.nan
    SpectralDB_Results['HLname'] = np.nan
    SpectralDB_Results['GLsmiles'] = np.nan
    SpectralDB_Results['GLname'] = np.nan
    SpectralDB_Results['MLsmiles'] = np.nan
    SpectralDB_Results['MLname'] = np.nan
    
    Suspect_list = pd.read_csv(slistcsv)
    
    for i, row in SpectralDB_Results.iterrows():
        if not isNaN(SpectralDB_Results['HMDBSMILES'][i]) and SpectralDB_Results['HMDBSMILES'][i] != " ":
            for j, row in Suspect_list.iterrows():
                LHms2 = [Chem.MolFromSmiles(SpectralDB_Results['HMDBSMILES'][i]), Chem.MolFromSmiles(Suspect_list['SMILES'][j])]
                LHfps2 = [AllChem.GetMorganFingerprintAsBitVect(x2,2, nBits=2048) for x2 in LHms2]
                LHtn2 = DataStructs.FingerprintSimilarity(LHfps2[0],LHfps2[1])
                if LHtn2 >= 0.9:
                    SpectralDB_Results.loc[i, 'HLsmiles'] = Suspect_list['SMILES'][j]
                    SpectralDB_Results.loc[i, 'HLname'] = Suspect_list['Name'][j]
                    
                    
        if not isNaN(SpectralDB_Results['GNPSSMILES'][i]) and SpectralDB_Results['GNPSSMILES'][i] != " ":
            for k, row in Suspect_list.iterrows():
                LGms2 = [Chem.MolFromSmiles(SpectralDB_Results['GNPSSMILES'][i]), Chem.MolFromSmiles(Suspect_list['SMILES'][k])]
                LGfps2 = [AllChem.GetMorganFingerprintAsBitVect(x2,2, nBits=2048) for x2 in LGms2]
                LGtn2 = DataStructs.FingerprintSimilarity(LGfps2[0],LGfps2[1])
                if LGtn2 >= 0.9:
                    SpectralDB_Results.loc[i, 'GLsmiles'] = Suspect_list['SMILES'][k]
                    SpectralDB_Results.loc[i, 'GLname'] = Suspect_list['Name'][k]
                    
                    
        if not isNaN(SpectralDB_Results['MBSMILES'][i]) and SpectralDB_Results['MBSMILES'][i] != " ":
            for l, row in Suspect_list.iterrows():
                LMms2 = [Chem.MolFromSmiles(SpectralDB_Results['MBSMILES'][i]), Chem.MolFromSmiles(Suspect_list['SMILES'][l])]
                LMfps2 = [AllChem.GetMorganFingerprintAsBitVect(x2,2, nBits=2048) for x2 in LMms2]
                LMtn2 = DataStructs.FingerprintSimilarity(LMfps2[0],LMfps2[1])
                if LMtn2 >= 0.9:
                    SpectralDB_Results.loc[i, 'MLsmiles'] = Suspect_list['SMILES'][l]
                    SpectralDB_Results.loc[i, 'MLname'] = Suspect_list['Name'][l]
    
    # add annotations and occurences
    for i, row in SpectralDB_Results.iterrows():
        if not isNaN(SpectralDB_Results['HLname'][i]) or not isNaN(SpectralDB_Results['GLname'][i]) or not isNaN(SpectralDB_Results['MLname'][i]):
            SpectralDB_Results['occurence'][i] = SpectralDB_Results['occurence'][i] + 1
            SpectralDB_Results['annotation'][i] = SpectralDB_Results['annotation'][i] + ', Suspect_List'
            
    SpectralDB_Results.to_csv(input_dir + "MetabolomicsResults/SpecDBvsSL.csv")
    return(SpectralDB_Results)


# In[22]:


#print(suspectListScreening.__doc__)


# In[ ]:





# # Final Candidate List Curation

# ## MetFrag Curation

# In[ ]:





# In[23]:


def metfrag_curation(input_dir, metfragcsv, sl = True):
    
    
    """metfrag_curation checks which database produced results. If both 
    did, it checks whether it was the same compound as candidate, if not,
    add PubChem or any of the two databases with similarity to Suspect
    list

    Parameters:
    input_dir (str): This is the input directory where all the .mzML 
    files and their respective result directories are stored.
    metfragcsv (str): path to combined metfrag results:
    MetabolomicsResults/MetFrag_combined.csv
    
    Returns:
    dataframe: dataframe of curated metfrag results
    csv: MetabolomicsResults/metfrag_curated.csv
    
    Usage:
    metfrag_curation(input_dir = "usr/project/", 
    metfragcsv = "usr/project/MetabolomicsResults/MetFrag_combined.csv")

    """
    
    metfrag = pd.read_csv(metfragcsv)
    for i, row in metfrag.iterrows():
    
        # If only Pubchem
        if not isNaN(metfrag['PC_SMILES'][i]) and isNaN(metfrag['KG_SMILES'][i]):
            metfrag.loc[i, 'Annotation_M'] = 'PubChem'
            #metfrag.loc[i, 'SMILES_final'] = metfrag['PC_SMILES'][i]
    
        # If only KEGG
        elif not isNaN(metfrag['KG_SMILES'][i]) and isNaN(metfrag['PC_SMILES'][i]):
            metfrag.loc[i, 'Annotation_M'] = 'KEGG'
            #metfrag.loc[i, 'SMILES_final'] = metfrag['KG_SMILES'][i]
    
        # If both, calculate the similarity
        elif not isNaN(metfrag['PC_SMILES'][i]) and not isNaN(metfrag['KG_SMILES'][i]):
        
            PKms = [Chem.MolFromSmiles(metfrag['KG_SMILES'][i]), Chem.MolFromSmiles(metfrag['PC_SMILES'][i])]
            PKfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in PKms]
            PKtn = DataStructs.FingerprintSimilarity(PKfps[0],PKfps[1])
        
            # if both are similar, add both
            if PKtn == 1:
                metfrag.loc[i, 'Annotation_M'] = 'KEGG, PubChem'
                #metfrag.loc[i, 'SMILES_final'] = metfrag['PC_SMILES'][i]
        
            #if not similar:
            else:
                # if there is NO entry from suspect list, then add PuBchem
                if isNaN(metfrag['KG_SL_comp'][i]) and isNaN(metfrag['PC_SL_comp'][i]):
                    metfrag.loc[i, 'Annotation_M'] = 'PubChem'
                    #metfrag.loc[i, 'SMILES_final'] = metfrag['PC_SMILES'][i]
                else:
                    if sl:
                        
                        #if there is an entry from suspect list WITH kegg
                        if not isNaN(metfrag['KG_SL_comp'][i]) and isNaN(metfrag['PC_SL_comp'][i]):
                            metfrag.loc[i, 'Annotation_M'] = 'KEGG, SuspectList'
                            #metfrag.loc[i, 'SMILES_final'] = metfrag['KG_SMILES'][i]
                
                        #if there is an entry from suspect list WITH kegg
                        elif not isNaN(metfrag['PC_SL_comp'][i]) and isNaN(metfrag['KG_SL_comp'][i]):
                            metfrag.loc[i, 'Annotation_M'] = 'PubChem, SuspectList'
                            #metfrag.loc[i, 'SMILES_final'] = metfrag['PC_SMILES'][i]
                    
                                
    metfrag.to_csv(input_dir + "MetabolomicsResults/metfrag_curated.csv")  
    return(metfrag)
    


# In[24]:


#print(metfrag_curation.__doc__)


# ## SIRIUS Results Curation

# In[25]:


def sirius_curation(input_dir, siriuscsv, sl = True):
    
    """sirius_curation checks if candidate selected has a good score for 
    explained intensity. It also checks if there was any similarity to
    a compound from Suspect list

    Parameters:
    input_dir (str): This is the input directory where all the .mzML 
    files and their respective result directories are stored.
    siriuscsv (str): path to combined metfrag results:
    MetabolomicsResults/Sirius_combined.csv
    
    Returns:
    dataframe: dataframe of curated sirius results
    csv: MetabolomicsResults/sirius_curated.csv
    
    Usage:
    sirius_curation(input_dir = "usr/project/", 
    siriuscsv = "usr/project/MetabolomicsResults/Sirius_combined.csv")

    """
    
    sirius = pd.read_csv(siriuscsv)
    for i, row in sirius.iterrows():
    
        # If the explained intensity is greater than 0.70 and there is no suspect list entry
        if sirius['exp_int'][i] >= 0.70 and isNaN(sirius['SL_comp'][i]):
            sirius.loc[i, 'Annotation_S'] = 'SIRIUS'
            #sirius.loc[i, 'SMILES_final'] = sirius['SMILES'][i]
        else:
            if sl:
                
                #If the explained intensity is greater than 0.70 and there is an entry from suspect list
                if sirius['exp_int'][i] >= 0.70 and not isNaN(sirius['SL_comp'][i]):
                    sirius.loc[i, 'Annotation_S'] = 'SIRIUS, SuspectList'
                    #sirius.loc[i, 'SMILES_final'] = sirius['SMILES'][i]
    
                # if the intensity is less thna 0.70 but it still is similar to an entry in Suspect list,
                elif sirius['exp_int'][i] < 0.70 and not isNaN(sirius['SL_comp'][i]):
                    sirius.loc[i, 'Annotation_S'] = 'SIRIUS, SuspectList'
                    #sirius.loc[i, 'SMILES_final'] = sirius['SMILES'][i]
        
    sirius.to_csv(input_dir + "MetabolomicsResults/sirius_curated.csv")
    return(sirius)


# In[26]:


#print(sirius_curation.__doc__)


# ## combine curated S and M results

# In[27]:


def combineSM(input_dir, metfrag, sirius):
    
    """combineSM prioritizes Sirius and Suspect list over PubChem and
    KEGG

    Parameters:
    input_dir (str): This is the input directory where all the .mzML 
    files and their respective result directories are stored.
    sirius (dataframe): result of sirius_curation
    metfrag (dataframe): result of metfrag_curation
    
    Returns:
    dataframe: dataframe of combined curated sirius and metfrag results
    csv: "MetabolomicsResults/combinedSM.csv"
    
    Usage:
    combineSM(input_dir = "usr/project/", metfrag, sirius)

    """

    S_M_CSV = pd.concat([sirius, metfrag], axis = 1, levels = ["id_X"])
    for i, rows in S_M_CSV.iterrows():
    
        # if results has Sirius Structure annotation, and the explained inetnsity is >= 0.70, keep the annotation as is.
        if S_M_CSV["Result"][i] == "SIRIUS_STR" and S_M_CSV['exp_int'][i] >= 0.70:
            S_M_CSV.loc[i, 'Annotation_C'] = S_M_CSV['Annotation_S'][i]
        
            # to add to that annotation
            if not isNaN(S_M_CSV["Annotation_M"][i]):
            
                # if annotation has PubChem, by default add SIRIUS
                if "PubChem" in S_M_CSV["Annotation_M"][i]:
                    S_M_CSV.loc[i, 'Annotation_C'] = S_M_CSV['Annotation_S'][i]
    
                    # And calculate the similarity between pubchem and SIrius results
                    try:
                        PSms = [Chem.MolFromSmiles(S_M_CSV['SMILES'][i]), Chem.MolFromSmiles(S_M_CSV['PC_SMILES'][i])]
                        PSfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in PSms]
                        PStn = DataStructs.FingerprintSimilarity(PSfps[0],PSfps[1])
                        # if similar strcutres, then add Pubchme and sirius
                        if PStn == 1:
                            print(i)
                            S_M_CSV.loc[i, 'Annotation_C'] = S_M_CSV['Annotation_S'][i] + ', PubChem'
                            #S_M_CSV.loc[i, 'SMILES_final'] = S_M_CSV['SMILES'][i]
                        # if not then just keep sirius
                        else:
                            S_M_CSV.loc[i, 'Annotation_C'] = S_M_CSV['Annotation_S'][i]
                            #S_M_CSV.loc[i, 'SMILES_final'] = S_M_CSV['SMILES'][i]
                    except:
                        pass
                #apply similar approach to KEGG, and by default add SIRIUS
                elif not isNaN(S_M_CSV["KG_SMILES"][i]):
                    S_M_CSV.loc[i, 'Annotation_C'] = S_M_CSV['Annotation_S'][i]

                    try:
                        SKms = [Chem.MolFromSmiles(S_M_CSV['SMILES'][i]), Chem.MolFromSmiles(S_M_CSV['KG_SMILES'][i])]
                        SKfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in SKms]
                        SKtn = DataStructs.FingerprintSimilarity(SKfps[0],SKfps[1])
                        if SKtn == 1:
                            print(i)
                            S_M_CSV.loc[i, 'Annotation_C'] = S_M_CSV['Annotation_S'][i] +', KEGG'

                        else:
                            S_M_CSV.loc[i, 'Annotation_C'] = S_M_CSV['Annotation_S'][i]

                    except:
                        pass
    
        # if there is no annotation from SIRIUS, then use Metfrag results
        elif S_M_CSV["Result"][i] != "SIRIUS_STR" and S_M_CSV['exp_int'][i] < 0.70:
            if not isNaN(S_M_CSV['Annotation_M'][i]):
                if 'PubChem' in S_M_CSV['Annotation_M'][i]:
                    S_M_CSV.loc[i, 'Annotation_C'] = S_M_CSV['Annotation_M'][i]
                    #S_M_CSV.loc[i, 'SMILES'] = S_M_CSV['PC_SMILES'][i]
                else:
                    S_M_CSV.loc[i, 'Annotation_C'] = S_M_CSV['Annotation_M'][i]
                    #S_M_CSV.loc[i, 'SMILES'] = S_M_CSV['KG_SMILES'][i]
                    
    S_M_CSV.to_csv(input_dir + "MetabolomicsResults/combinedSM.csv")
    return(S_M_CSV)


# In[28]:


#print(combineSM.__doc__)


# ## Spec DB Curation

# In[29]:


def specDB_Curation(input_dir, combinedx, sl = True):
    
    """specDB_Curation prioritizes in the following manner: gnps>
    mbank>suspectlist>hmdb

    Parameters:
    input_dir (str): This is the input directory where all the .mzML 
    files and their respective result directories are stored.
    
    combined: dataframe from either suspectListScreening function if
    sl = True OR from scoring_spec if sl = False
    
    Returns:
    dataframe: with curated Spectral DB results
    csv: "MetabolomicsResults/curatedSDB.csv"
    
    Usage:
    specDB_Curation(input_dir = "usr/project/",combinedx, sl = True)

    """

    # define scores again to remove the entries with lower scores
    def HMDB_Scoring(db, i):
        if db['HMDBmax_similarity'][i] >= 0.75 and db['HMDBintScore'][i] >= 0.50 and db['HMDBmzScore'][i] >= 0.50 and db['HQMatchingPeaks'][i]/db['hQueryTotalPeaks'][i] >= 0.50:
            return True
        else:
            return False
    def GNPS_Scoring(db, i):
        if db['GNPSmax_similarity'][i] >= 0.75 and db['GNPSintScore'][i] >= 0.50 and db['GNPSmzScore'][i] >= 0.50 and db['GQMatchingPeaks'][i]/db['gQueryTotalPeaks'][i] >= 0.50:
            return True
        else:
            return False
    def MB_Scoring(db, i):
        if db['MBmax_similarity'][i] >= 0.75 and db['MBintScore'][i] >= 0.50 and db['MBmzScore'][i] >= 0.50 and db['MQMatchingPeaks'][i]/db['mQueryTotalPeaks'][i] >= 0.50:
            return True
        else:
            return False
    
    
    combined = pd.read_csv(combinedx)
    
    
    # remove the similarity scores from low scoring candidates
    for i, row in combined.iterrows():
        if HMDB_Scoring(combined, i) or GNPS_Scoring(combined, i) or MB_Scoring(combined, i):
            print(i)
        else:
            combined['GNPSspectrumID'][i] = np.nan
            combined['MBspectrumID'][i] = np.nan
            combined['HMDBcompoundID'][i] = np.nan
    
    # if sl = True
    if sl:
    
        for i, row in combined.iterrows():
    
        ##### When there is an annotaion from all DBs #####
    
        #all entries with a high scoring annotation in all DBs,
            if not isNaN(combined['GNPSspectrumID'][i]) and not isNaN(combined['MBspectrumID'][i]) and not isNaN(combined['HMDBcompoundID'][i]):
        
                # entries with same candidate from all Spectral DBs
                if combined['tanimotoHG'][i] == 1.0 and combined['tanimotoGM'][i] == 1.0 and combined['tanimotoHM'][i] == 1.0:
                    combined.loc[i, 'Annotation'] = 'GNPS, HMDB, MassBank'
                    #entries with same candidate in suspect list, as in all Spectral DBs
                    if combined['GLname'][i] == combined['HLname'][i]== combined['MLname'][i]:
                        combined.loc[i, 'Annotation'] = 'GNPS, HMDB, MassBank, SuspectList'
                
                # same candidate from GNPS and HMDB        
                if combined['tanimotoHG'][i] == 1.0 and combined['tanimotoGM'][i] != 1.0 and combined['tanimotoHM'][i] != 1.0:
                    combined.loc[i, 'Annotation'] = 'GNPS, HMDB'
                    # if its present in Suspect List
                    if combined['GLname'][i] == combined['HLname'][i]:
                        combined.loc[i, 'Annotation'] = 'GNPS, HMDB, SuspectList'
        
                # same candidate from GNPS and MassBank        
                if combined['tanimotoHG'][i] != 1.0 and combined['tanimotoGM'][i] == 1.0 and combined['tanimotoHM'][i] != 1.0:
                    combined.loc[i, 'Annotation'] = 'GNPS, MassBank'
                    # if its present in Suspect List
                    if combined['GLname'][i] == combined['MLname'][i]:
                        combined.loc[i, 'Annotation'] = 'GNPS, MassBank, SuspectList'
                
                # same candidate from MassBank and HMDB        
                if combined['tanimotoHG'][i] != 1.0 and combined['tanimotoGM'][i] != 1.0 and combined['tanimotoHM'][i] == 1.0:
                    combined.loc[i, 'Annotation'] = 'GNPS, HMDB'
                    # if its present in Suspect List
                    if combined['GLname'][i] == combined['HLname'][i]:
                        combined.loc[i, 'Annotation'] = 'HMDB, MassBank, SuspectList'
                
                # only one database must be selected based on SuspectList annotation
                if combined['tanimotoHG'][i] != 1.0 and combined['tanimotoGM'][i] != 1.0 and combined['tanimotoHM'][i] != 1.0:
            
                    # only GNPS has SuspectList annotation
                    if not isNaN(combined['GLname'][i]):

                        combined.loc[i, 'Annotation'] = 'GNPS, SuspectList'
            
            
                    # only MassBank has SuspectList annotation
                    elif not isNaN(combined['MLname'][i]):
                        combined.loc[i, 'Annotation'] = 'MassBank, SuspectList'
            
            
                    # only HMDB has SuspectList annotation
                    elif not isNaN(combined['HLname'][i]):
                        combined.loc[i, 'Annotation'] = 'HMDB, SuspectList'
            
        
                    # all different annotations, take GNPS
                    else:
                        if not isNaN(combined['GNPSSMILES'][i]):
                            combined.loc[i, 'Annotation'] = 'GNPS'
                        else:
                            combined.loc[i, 'Annotation'] = 'MassBank'
            
        
            ##### When there is an annotation from two DBs #####
    
    
            # only GNPS and HMDB
            if not isNaN(combined['GNPSspectrumID'][i]) and isNaN(combined['MBspectrumID'][i]) and not isNaN(combined['HMDBcompoundID'][i]):

                if not isNaN(combined['GLname'][i]):
                    combined.loc[i, 'Annotation'] = 'GNPS, SuspectList'
                elif not isNaN(combined['HLname'][i]):
                    combined.loc[i, 'Annotation'] = 'HMDB, SuspectList'
                elif not isNaN(combined['GNPSSMILES'][i]):
                    combined.loc[i, 'Annotation'] = 'GNPS'
                else:
                    combined.loc[i, 'Annotation'] = 'HMDB'
        
    
            # only GNPS and MassBank
            if not isNaN(combined['GNPSspectrumID'][i]) and not isNaN(combined['MBspectrumID'][i]) and isNaN(combined['HMDBcompoundID'][i]):

                if not isNaN(combined['GLname'][i]):
                    combined.loc[i, 'Annotation'] = 'GNPS, SuspectList'
                elif not isNaN(combined['MLname'][i]):
                    combined.loc[i, 'Annotation'] = 'MassBank, SuspectList'
                elif not isNaN(combined['GNPSSMILES'][i]):
                    combined.loc[i, 'Annotation'] = 'GNPS'
                else:
                    combined.loc[i, 'Annotation'] = 'MassBank'
    
            # only MassBank and HMDB
    
            if isNaN(combined['GNPSspectrumID'][i]) and not isNaN(combined['MBspectrumID'][i]) and not isNaN(combined['HMDBcompoundID'][i]):

                if not isNaN(combined['MLname'][i]):
                    combined.loc[i, 'Annotation'] = 'MassBank, SuspectList'
                elif not isNaN(combined['HLname'][i]):
                    combined.loc[i, 'Annotation'] = 'HMDB, SuspectList'
                else:
                    combined.loc[i, 'Annotation'] = 'MassBank'
    
    
    
            ##### When there is an annotation from one DBs #####
    
    
            # only GNPS
            if not isNaN(combined['GNPSspectrumID'][i]) and isNaN(combined['MBspectrumID'][i]) and isNaN(combined['HMDBcompoundID'][i]):
        
                #If also SuspectList
                if not isNaN(combined['GLname'][i]):
                    combined.loc[i, 'Annotation'] = 'GNPS, SuspectList'
                elif not isNaN(combined['GNPSSMILES'][i]):
                    combined.loc[i, 'Annotation'] = 'GNPS'
        
            # only MassBank
            if isNaN(combined['GNPSspectrumID'][i]) and not isNaN(combined['MBspectrumID'][i]) and isNaN(combined['HMDBcompoundID'][i]):
                combined.loc[i, 'Annotation'] = 'MassBank'
                #If also SuspectList
                if not isNaN(combined['MLname'][i]):
                    combined.loc[i, 'Annotation'] = 'MassBank, SuspectList'
    
            # only HMDB
            if isNaN(combined['GNPSspectrumID'][i]) and isNaN(combined['MBspectrumID'][i]) and not isNaN(combined['HMDBcompoundID'][i]):
                combined.loc[i, 'Annotation'] = 'HMDB'
                #If also SuspectList
                if not isNaN(combined['HLname'][i]):
                    combined.loc[i, 'Annotation'] = 'HMDB, SuspectList'
            
    else:
        for i, row in combined.iterrows():
    
        ##### When there is an annotaion from all DBs #####
    
        #all entries with a high scoring annotation in all DBs,
            if not isNaN(combined['GNPSspectrumID'][i]) and not isNaN(combined['MBspectrumID'][i]) and not isNaN(combined['HMDBcompoundID'][i]):
        
                # entries with same candidate from all Spectral DBs
                if combined['tanimotoHG'][i] == 1.0 and combined['tanimotoGM'][i] == 1.0 and combined['tanimotoHM'][i] == 1.0:
                    combined.loc[i, 'Annotation'] = 'GNPS, HMDB, MassBank'
                
                # same candidate from GNPS and HMDB        
                if combined['tanimotoHG'][i] == 1.0 and combined['tanimotoGM'][i] != 1.0 and combined['tanimotoHM'][i] != 1.0:
                    combined.loc[i, 'Annotation'] = 'GNPS, HMDB'
        
                # same candidate from GNPS and MassBank        
                if combined['tanimotoHG'][i] != 1.0 and combined['tanimotoGM'][i] == 1.0 and combined['tanimotoHM'][i] != 1.0:
                    combined.loc[i, 'Annotation'] = 'GNPS, MassBank'
                
                # same candidate from MassBank and HMDB        
                if combined['tanimotoHG'][i] != 1.0 and combined['tanimotoGM'][i] != 1.0 and combined['tanimotoHM'][i] == 1.0:
                    combined.loc[i, 'Annotation'] = 'GNPS, HMDB'
                
                
            # all different annotations, take GNPS
            else:
                if not isNaN(combined['GNPSSMILES'][i]):
                    combined.loc[i, 'Annotation'] = 'GNPS'
                else:
                    combined.loc[i, 'Annotation'] = 'MassBank'
            
        
            ##### When there is an annotation from two DBs #####
    
    
            # only GNPS and HMDB
            if not isNaN(combined['GNPSspectrumID'][i]) and isNaN(combined['MBspectrumID'][i]) and not isNaN(combined['HMDBcompoundID'][i]):

                if not isNaN(combined['GNPSSMILES'][i]):
                    combined.loc[i, 'Annotation'] = 'GNPS'
                else:
                    combined.loc[i, 'Annotation'] = 'HMDB'
        
    
            # only GNPS and MassBank
            if not isNaN(combined['GNPSspectrumID'][i]) and not isNaN(combined['MBspectrumID'][i]) and isNaN(combined['HMDBcompoundID'][i]):

                if not isNaN(combined['GNPSSMILES'][i]):
                    combined.loc[i, 'Annotation'] = 'GNPS'
                else:
                    combined.loc[i, 'Annotation'] = 'MassBank'
    
            # only MassBank and HMDB
            if isNaN(combined['GNPSspectrumID'][i]) and not isNaN(combined['MBspectrumID'][i]) and not isNaN(combined['HMDBcompoundID'][i]):
                combined.loc[i, 'Annotation'] = 'MassBank'
    
    
    
            ##### When there is an annotation from one DBs #####
    
    
            # only GNPS
            if not isNaN(combined['GNPSspectrumID'][i]) and isNaN(combined['MBspectrumID'][i]) and isNaN(combined['HMDBcompoundID'][i]):
                if not isNaN(combined['GNPSSMILES'][i]):
                    combined.loc[i, 'Annotation'] = 'GNPS'
        
            # only MassBank
            if isNaN(combined['GNPSspectrumID'][i]) and not isNaN(combined['MBspectrumID'][i]) and isNaN(combined['HMDBcompoundID'][i]):
                combined.loc[i, 'Annotation'] = 'MassBank'
    
            # only HMDB
            if isNaN(combined['GNPSspectrumID'][i]) and isNaN(combined['MBspectrumID'][i]) and not isNaN(combined['HMDBcompoundID'][i]):
                combined.loc[i, 'Annotation'] = 'HMDB'
                
    combined.to_csv(input_dir + "MetabolomicsResults/curatedSDB.csv")
    return(combined)
    


# In[30]:


#print(specDB_Curation.__doc__)


# # combine curated SDB and CDB (S+M)

# In[ ]:





# In[31]:


def combine_CuratedR(input_dir, combinedSDBs, combinedSMs):
    
    """combine_CuratedR prioritizes in the following manner: gnps>
    mbank>suspectlist>sirius>hmdb>metfrag

    Parameters:
    input_dir (str): This is the input directory where all the .mzML 
    files and their respective result directories are stored.
    
    curatedSDB: df from specDB_Curation
    combinedSM: df from combineSM
    
    Returns:
    dataframe: with curated Spectral DB results and CDB (S+M) results
    csv: "MetabolomicsResults/final_curation_without_classes.csv"
    
    Usage:
    combine_CuratedR(input_dir = "usr/project/", curatedSDB, combinedSM)

    """

    combinedSDB = pd.read_csv(combinedSDBs)
    combinedSM = pd.read_csv(combinedSMs)
    
    mega = pd.concat([combinedSM, combinedSDB], axis = 1, levels = ["id_X"])
    
    for i, row in mega.iterrows():
    
        #if only compound database results
        if isNaN(mega['Annotation'][i]) and not isNaN(mega['Annotation_C'][i]):
            mega.loc[i, "Annotation_Source"] = mega['Annotation_C'][i]
        
        # if only spectral db results
        if not isNaN(mega['Annotation'][i]) and isNaN(mega['Annotation_C'][i]):
            mega.loc[i, "Annotation_Source"] = mega['Annotation'][i]
        
        
        # if both have results
        if not isNaN(mega['Annotation'][i]) and not isNaN(mega['Annotation_C'][i]):
        
            ########THREE OR FOUR SDB SOURCES########
        
            #if three sdb sources or more
            # prioritize Spectral DBs
            if len(mega['Annotation'][i].split()) >= 3 and 'SIRIUS' in mega['Annotation_C'][i]:
            
                if 'MassBank' in mega['Annotation'][i]:
                    SKms = [Chem.MolFromSmiles(mega['MBSMILES'][i]), Chem.MolFromSmiles(mega['SMILES'][i])]
                    SKfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in SKms]
                    SKtn = DataStructs.FingerprintSimilarity(SKfps[0],SKfps[1])
                    if SKtn == 1.0:
                        print(SKtn)
                        mega.loc[i, "Annotation_Source"] = mega['Annotation'][i] + ', ' + mega['Annotation_C'][i]
                    else:
                        mega.loc[i, "Annotation_Source"] = mega['Annotation'][i]
            
                elif 'HMDB' in mega['Annotation'][i]:
                    SKms = [Chem.MolFromSmiles(mega['HMDBSMILES'][i]), Chem.MolFromSmiles(mega['SMILES'][i])]
                    SKfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in SKms]
                    SKtn = DataStructs.FingerprintSimilarity(SKfps[0],SKfps[1])
                    if SKtn == 1.0:
                        mega.loc[i, "Annotation_Source"] = mega['Annotation'][i] + ', ' + mega['Annotation_C'][i]
                    else:
                        mega.loc[i, "Annotation_Source"] = mega['Annotation'][i]
                    
            elif len(mega['Annotation'][i].split()) >= 3 and 'SIRIUS' not in mega['Annotation_C'][i]:
                mega.loc[i, "Annotation_Source"] = mega['Annotation'][i]
            
            
            
            ########TWO SDB SOURCES########
            
            #if two sdb sources not HMDB
            #still prioritize Spectral DBs
            elif len(mega['Annotation'][i].split()) == 2 and 'HMDB' not in mega['Annotation'][i] and 'SIRIUS' in mega['Annotation_C'][i]:
                if 'MassBank' in mega['Annotation'][i]:
                    SKms = [Chem.MolFromSmiles(mega['MBSMILES'][i]), Chem.MolFromSmiles(mega['SMILES'][i])]
                    SKfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in SKms]
                    SKtn = DataStructs.FingerprintSimilarity(SKfps[0],SKfps[1])
                    if SKtn == 1.0:
    
                        mega.loc[i, "Annotation_Source"] = mega['Annotation'][i] + ', ' + mega['Annotation_C'][i]
                    else:
                        mega.loc[i, "Annotation_Source"] = mega['Annotation'][i]
            
                #for GNPS also check if there are any smiles given
                elif 'GNPS' in mega['Annotation'][i] and not isNaN(mega['GNPSSMILES'][i]):
                    SKms = [Chem.MolFromSmiles(mega['GNPSSMILES'][i]), Chem.MolFromSmiles(mega['SMILES'][i])]
                    SKfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in SKms]
                    SKtn = DataStructs.FingerprintSimilarity(SKfps[0],SKfps[1])
                    if SKtn == 1:
                        print(SKtn)
                        mega.loc[i, "Annotation_Source"] = mega['Annotation'][i] + ', ' + mega['Annotation_C'][i]
                    else:
                        mega.loc[i, "Annotation_Source"] = mega['Annotation'][i]
                    
                elif 'GNPS' in mega['Annotation'][i] and isNaN(mega['GNPSSMILES'][i]):
                    mega.loc[i, "Annotation_Source"] = mega['Annotation'][i]
                
            elif len(mega['Annotation'][i].split()) == 2 and 'HMDB' not in mega['Annotation'][i] and 'SIRIUS' not in mega['Annotation_C'][i]:
                mega.loc[i, "Annotation_Source"] = mega['Annotation'][i]
            
            
            #if two sdb sources but with HMDB
            #prioritize SIRIUS
            elif len(mega['Annotation'][i].split()) == 2 and 'HMDB' in mega['Annotation'][i] and 'SIRIUS' in mega['Annotation_C'][i]:
                mega.loc[i, "Annotation_Source"] = mega['Annotation_C'][i]
            
            
            # if two sdb sources with HMDB but cant prioritize SIRIUS
            #prioritize Spectral DB
            elif len(mega['Annotation'][i].split()) == 2 and 'HMDB' not in mega['Annotation'][i] and 'SIRIUS' in mega['Annotation_C'][i]:

                mega.loc[i, "Annotation_Source"] = mega['Annotation'][i]
        
        
        
            #######ONE SDB SOURCE#########
        
            #if only one sdb, but no SIRIUS
            #prioritize Spectral DB
            elif len(mega['Annotation'][i].split()) == 1 and 'SIRIUS' not in mega['Annotation_C'][i]:

                mega.loc[i, "Annotation_Source"] = mega['Annotation'][i]
        
            #if only one sdb
            #prioritize SIRIUS
            elif len(mega['Annotation'][i].split()) == 1 and 'SIRIUS' in mega['Annotation_C'][i]:
                mega.loc[i, "Annotation_Source"] = mega['Annotation_C'][i]
            
            
        # if only spectral db results
        if isNaN(mega['Annotation'][i]) and isNaN(mega['Annotation_C'][i]) and not isNaN(mega['Formula'][i]):
            mega.loc[i, "Annotation_Source"] = 'SIRIUS_Formula'
    
    
    
    bef_mega = mega.loc[:,~mega.columns.duplicated()]
    
    for i, row in bef_mega.iterrows():
        if not isNaN(bef_mega['Annotation_Source'][i]):
        
            if 'SIRIUS' in bef_mega['Annotation_Source'][i] and 'SIRIUS_Formula' not in bef_mega['Annotation_Source'][i]:
                bef_mega.loc[i, 'SMILES_final'] = bef_mega['SMILES'][i]
                bef_mega.loc[i,"CompoundNames"] = bef_mega['name'][i]
                bef_mega['PC_MCSS_SMILES'][i] = np.nan
                bef_mega['KG_MCSS_SMILES'][i] = np.nan
            
            
            elif 'MassBank' in bef_mega['Annotation_Source'][i]:
                bef_mega.loc[i, 'SMILES_final'] = bef_mega['MBSMILES'][i]
                bef_mega.loc[i, 'CompoundNames'] = bef_mega['MBcompound_name'][i]
                bef_mega['most_specific_class'][i] = np.nan
                bef_mega['level _5'][i] = np.nan
                bef_mega['subclass'][i] = np.nan
                bef_mega['class'][i] = np.nan
                bef_mega['superclass'][i] = np.nan
                bef_mega['all_classifications'][i] = np.nan
                bef_mega['Classification_Source'][i] = np.nan
                bef_mega['MCSS_SMILES'][i] = np.nan
                bef_mega['PC_MCSS_SMILES'][i] = np.nan
                bef_mega['KG_MCSS_SMILES'][i] = np.nan
                bef_mega['Formula'][i] = np.nan
            
            elif 'HMDB' in bef_mega['Annotation_Source'][i]:
                bef_mega.loc[i, 'SMILES_final'] = bef_mega['HMDBSMILES'][i]
                bef_mega.loc[i, 'CompoundNames'] = bef_mega['HMDBcompound_name'][i]
                bef_mega['most_specific_class'][i] = np.nan
                bef_mega['level _5'][i] = np.nan
                bef_mega['subclass'][i] = np.nan
                bef_mega['class'][i] = np.nan
                bef_mega['superclass'][i] = np.nan
                bef_mega['all_classifications'][i] = np.nan
                bef_mega['Classification_Source'][i] = np.nan
                bef_mega['MCSS_SMILES'][i] = np.nan
                bef_mega['PC_MCSS_SMILES'][i] = np.nan
                bef_mega['KG_MCSS_SMILES'][i] = np.nan
                bef_mega['Formula'][i] = np.nan
            
            elif 'GNPS, SuspectList' in bef_mega['Annotation_Source'][i]:
                bef_mega.loc[i,'SMILES_final'] = bef_mega['GLsmiles'][i]
                bef_mega.loc[i, 'CompoundNames'] = bef_mega['GLname'][i]
                bef_mega.loc[i, 'CompoundNames']
                bef_mega['most_specific_class'][i] = np.nan
                bef_mega['level _5'][i] = np.nan
                bef_mega['subclass'][i] = np.nan
                bef_mega['class'][i] = np.nan
                bef_mega['superclass'][i] = np.nan
                bef_mega['all_classifications'][i] = np.nan
                bef_mega['MCSS_SMILES'][i] = np.nan
                bef_mega['Classification_Source'][i] = np.nan
                bef_mega['PC_MCSS_SMILES'][i] = np.nan
                bef_mega['KG_MCSS_SMILES'][i] = np.nan
                bef_mega['Formula'][i] = np.nan
        
            elif 'GNPS' in bef_mega['Annotation_Source'][i]:
                bef_mega.loc[i,'SMILES_final'] = bef_mega['GNPSSMILES'][i]
                bef_mega.loc[i, 'CompoundNames'] = bef_mega['GNPScompound_name'][i]
                bef_mega['most_specific_class'][i] = np.nan
                bef_mega['level _5'][i] = np.nan
                bef_mega['subclass'][i] = np.nan
                bef_mega['class'][i] = np.nan
                bef_mega['superclass'][i] = np.nan
                bef_mega['all_classifications'][i] = np.nan
                bef_mega['MCSS_SMILES'][i] = np.nan
                bef_mega['Classification_Source'][i] = np.nan
                bef_mega['PC_MCSS_SMILES'][i] = np.nan
                bef_mega['KG_MCSS_SMILES'][i] = np.nan
                bef_mega['Formula'][i] = np.nan
            
            elif 'PubChem' in bef_mega['Annotation_Source'][i]:
                bef_mega.loc[i, 'SMILES_final'] = bef_mega['PC_SMILES'][i]
                bef_mega.loc[i, 'CompoundNames'] = bef_mega['PC_Name'][i]
                bef_mega['most_specific_class'][i] = np.nan
                bef_mega['level _5'][i] = np.nan
                bef_mega['subclass'][i] = np.nan
                bef_mega['class'][i] = np.nan
                bef_mega['superclass'][i] = np.nan
                bef_mega['all_classifications'][i] = np.nan
                bef_mega['Classification_Source'][i] = np.nan
                bef_mega['MCSS_SMILES'][i] = np.nan
                bef_mega['KG_MCSS_SMILES'][i] = np.nan
                bef_mega['Formula'][i] = np.nan

            elif 'KEGG' in bef_mega['Annotation_Source'][i]:
                bef_mega.loc[i, 'SMILES_final'] = bef_mega['KG_SMILES'][i]
                bef_mega.loc[i, 'CompoundNames'] = bef_mega['KG_Name'][i]
                bef_mega['most_specific_class'][i] = np.nan
                bef_mega['level _5'][i] = np.nan
                bef_mega['subclass'][i] = np.nan
                bef_mega['class'][i] = np.nan
                bef_mega['superclass'][i] = np.nan
                bef_mega['all_classifications'][i] = np.nan
                bef_mega['Classification_Source'][i] = np.nan
                bef_mega['MCSS_SMILES'][i] = np.nan
                bef_mega['PC_MCSS_SMILES'][i] = np.nan
                bef_mega['Formula'][i] = np.nan
            
            elif 'SIRIUS_Formula' in bef_mega['Annotation_Source'][i]:
                bef_mega['PC_MCSS_SMILES'][i] = np.nan
                bef_mega['KG_MCSS_SMILES'][i] = np.nan
                
    
    
    
    bef_megaA = bef_mega[['id_X', 
                          'premz', 
                          'rtmed', 
                          'rtmean',
                          'int', 
                          'col_eng', 
                          'pol', 
                          'SMILES_final', 
                          'CompoundNames', 
                          'MCSS_SMILES', 
                          'PC_MCSS_SMILES', 
                          'KG_MCSS_SMILES', 
                          'subclass', 
                          'class', 
                          'superclass', 
                          'Classification_Source', 
                          'Annotation_Source'
                         ]]
    bef_megaA.to_csv(input_dir + "MetabolomicsResults/final_curation_without_classes.csv")
    return(bef_megaA)
    


# In[32]:


#print(combine_CuratedR.__doc__)


# In[33]:


def checkSMILES_validity(input_dir, resultcsv):
    
    """checkSMILES_validity does exactly as the name says, using
    RDKit, whether the SMILES are invalid or have invalid 
    chemistry

    Parameters:
    input_dir (str): This is the input directory where all the .mzML 
    files and their respective result directories are stored.
    
    results: df from combine_CuratedR
    
    Returns:
    dataframe: with valid SMILES
    csv: "MetabolomicsResults/final_curation_with_validSMILES.csv"
    
    Usage:
    checkSMILES_validity(input_dir = "usr/project/", results)

    """
    results = pd.read_csv(resultcsv)
    # check validity of SMILES
    for i, row in results.iterrows():
        if not isNaN(results['SMILES_final'][i]):
            m = Chem.MolFromSmiles(results['SMILES_final'][i] ,sanitize=False)
            if m is None:
                results['SMILES_final'][i] = 'invalid_SMILES'
            else:
                try:
                    Chem.SanitizeMol(m)
                except:
                    results['SMILES_final'][i] = 'invalid_chemistry'
    results.to_csv(input_dir + "MetabolomicsResults/final_curation_with_validSMILES.csv")
    return(results)


# In[34]:


#print(checkSMILES_validity.__doc__)


# In[35]:


def classification(input_dir, resultcsv):
    
    """classification function uses ClassyFire ChemONT

    Parameters:
    input_dir (str): This is the input directory where all the .mzML 
    files and their respective result directories are stored.
    
    resultcsv: csv of df from combine_CuratedR or checkSMILES_validity
    
    Returns:
    dataframe: with classification
    csv: "MetabolomicsResults/final_curationList.csv"
    
    Usage:
    checkSMILES_validity(input_dir = "usr/project/", frame)

    """
    frame = pd.read_csv(resultcsv)
    inchis = []
    for i, row in frame.iterrows():
        if not isNaN(frame['SMILES_final'][i]) and isNaN(frame['Classification_Source'][i]):
            try:
                InChI = Chem.MolToInchi(Chem.MolFromSmiles(frame["SMILES_final"][i]))
                InChIKey = Chem.inchi.InchiToInchiKey(InChI)
                inchis.append({
                    'index': i,
                    'smiles':frame["SMILES_final"][i],
                    'inchi': InChI,
                    'inchikey': InChIKey
                })
            except:
                pass
    inchis = pd.DataFrame(inchis)
    inchis = inchis.loc[-isNaN(inchis['inchikey'])]
    ## Retrieve ClassyFire classifications ##
    
    # This first step is done using inchikey and interrogation of the gnps classified structures
    gnps_proxy = True 
    url = "http://classyfire.wishartlab.com"
    proxy_url =  "https://gnps-classyfire.ucsd.edu"
    chunk_size = 1000
    sleep_interval = 12
    
    all_inchi_keys = list(inchis['inchikey'].drop_duplicates())

    resolved_ik_number_list = [0, 0]
    total_inchikey_number = len(all_inchi_keys)

    while True:
    
        #start_time = time.time()
    
        #print('%s inchikey to resolve' % total_inchikey_number )
        get_classifications_cf_mod(all_inchi_keys, par_level = 6)
    
        cleanse('all_json.json', 'all_json.json')
    
        with open("all_json.json") as tweetfile:
            jsondic = json.loads(tweetfile.read())

        df = json_normalize(jsondic)
        df = df.drop_duplicates( 'inchikey' )
        resolved_ik_number = len( df.drop_duplicates('inchikey').inchikey )
        resolved_ik_number_list.append( resolved_ik_number )
        #print('%s resolved inchikeys' % resolved_ik_number )
        #print("done in --- %s seconds ---" % (time.time() - start_time))
    
        if resolved_ik_number_list[-1] < resolved_ik_number_list[-2] or resolved_ik_number_list[-1] == resolved_ik_number_list[-3]:
            break
        cleanse('all_json.json', 'all_json_cleaned.json')
        
        with open("all_json_cleaned.json") as tweetfile:
            jsondic = json.loads(tweetfile.read())
            
    flattened_classified_json = json_normalize(jsondic)
    flattened_df = flattened_classified_json.drop_duplicates('inchikey')
    flattened_df['inchikey'] = flattened_df['inchikey'].str.replace(r'InChIKey=', '')
    df_merged = pd.merge(inchis, flattened_df, left_on='inchikey', right_on='inchikey', how='left')
    
    for p, rowp in df_merged.iterrows():
        for q, rowq in frame.iterrows():
            if df_merged["smiles_x"][p] is frame["SMILES_final"][q]:
                frame.loc[q, 'subclass'] = df_merged["subclass.name"][p]
                frame.loc[q, 'class'] = df_merged["class.name"][p]
                frame.loc[q, 'superclass'] = df_merged["superclass.name"][p]
                frame.loc[q, 'Classification_Source'] = "ClassyFire"
    #frame.to_csv(input_dir, '/SIRIUS_combined.csv')
    


    frame.to_csv(input_dir + "MetabolomicsResults/final_curationList.csv")
    return(frame)


# In[36]:


#print(classification.__doc__)


# # Comparison with a list of SMILES from any Source

# In[37]:


def SMILESscreening(input_dir, results, listname):
    
    """SMILESscreening takes a list of SMILES

    Parameters:
    input_dir (str): This is the input directory where all the .mzML 
    files and their respective result directories are stored.
    
    results: df from combine_CuratedR or checkSMILES_validity or classification
    
    Returns:
    dataframe: comparison with another list of compounds
    csv: "MetabolomicsResults/final_curation_with_validSMILES.csv"
    
    Usage:
    checkSMILES_validity(input_dir = "usr/project/", results)

    """
    
    for i, row in results.iterrows():
        if not isNaN(results['SMILES_final'][i]):
            if 'invalid_SMILES' not in results['SMILES_final'][i] and 'invalid_chemistry' not in results['SMILES_final'][i]:
                for j in cd:
                    if not isNaN(j):
                        CGms = [Chem.MolFromSmiles(results['SMILES_final'][i]), Chem.MolFromSmiles(j)]
                        CGfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=1024) for x in CGms]
                        CGtn = DataStructs.FingerprintSimilarity(CGfps[0],CGfps[1])
                        if CGtn == 1 and 'CompoundDiscoverer' not in results['Annotation_Source'][i]:
                            results['Annotation_Source'][i] = results['Annotation_Source'][i] + ', ' + listname
                            results['Occurence'][i] = results['Occurence'][i] + 1
    return(frame)

    frame.to_csv(input_dir + "MetabolomicsResults/final_curationListVS"+listname+".csv")


# In[38]:


#print(SMILESscreening.__doc__)


# In[ ]:





# In[ ]:





# In[ ]:




