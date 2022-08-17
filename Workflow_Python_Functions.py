#!/usr/bin/env python
# coding: utf-8

# ## SIRIUS_Metfrag_SList

# In[93]:


from platform import python_version

#print(python_version())


# In[94]:


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
import wget
import string
import urllib.parse
import openpyxl
import statistics
import sys
from itertools import chain


# In[95]:


from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem import PandasTools


# In[4]:


input_dir = "/Users/mahnoorzulfiqar/OneDriveUNI/MAW-data/StandardSMarinoi_Data"
input_dir


# # Suspect List for MetFrag and SIRIUS

# In[5]:


# make sure your Smiles entries in the suspect list csv are in a column named "SMILES"
def slist_metfrag(input_dir, slist_csv, name):
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
        if i is not None:
            mols = Chem.MolFromSmiles(sl['SMILES'][i])
            try:
                sl.loc[i, 'InChIKey'] = Chem.inchi.MolToInchiKey(mols)
                sl_mtfrag.append(sl['InChIKey'][i])
            except Exception as e:
                print(e)
    
    with open((input_dir + "/SL_"+ name + '.txt'), 'w') as filehandle:
        for listitem in sl_mtfrag:
            filehandle.write('%s\n' % listitem)
    return(sl_mtfrag)


# In[6]:


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


# In[7]:


#print(slist_metfrag.__doc__)
#print(slist_sirius.__doc__)


# ## Spectral DB dereplication Results PostProcessing

# In[8]:


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
    
    def isNaN(string):
        return string != string
    
    
    def HMDB_Scoring(db, i):
        if db['HMDBintScore'][i] >= 0.50 and db['HMDBmzScore'][i] >= 0.50 and db['HQMatchingPeaks'][i]/db['hQueryTotalPeaks'][i] >= 0.50:
            return True
        else:
            return False


    def GNPS_Scoring(db, i):
        if db['GNPSintScore'][i] >= 0.50 and db['GNPSmzScore'][i] >= 0.50 and db['GQMatchingPeaks'][i]/db['gQueryTotalPeaks'][i] >= 0.50:
            return True
        else:
            return False


    def MB_Scoring(db, i):
        if db['MBintScore'][i] >= 0.50 and db['MBmzScore'][i] >= 0.50 and db['MQMatchingPeaks'][i]/db['mQueryTotalPeaks'][i] >= 0.50:
            return True
        else:
            return False
    
    
    
    #list all files and directories
    for entry in os.listdir(input_dir):
        if os.path.isdir(os.path.join(input_dir, entry)):
            #print(entry)
            msp_file = (glob.glob(input_dir + "/" + entry + '/spectral_dereplication' +'/*.csv'))
            #print(msp_file)
            if len(msp_file) > 0:
                if os.path.exists(msp_file[0]):
                    msp = pd.read_csv(msp_file[0])
                    # enter the directory with /spectral_dereplication/ results
                    
                    # enter the directory with /spectral_dereplication/ results
                    if Source == "gnps" or Source == "all":
                        #currently only these subsets are removed from the names from GNPS
                        matches = ["M+","[M", "M-", "2M", "M*" "20.0", "50.0", "30.0", "40.0", "60.0", "70.0", "eV", "Massbank"
                                   , "Spectral", "Match", "to", "from", "NIST14", "MoNA", '[IIN-based:',  '[IIN-based', 'on:', 'CCMSLIB00003136269]']

                        #open another csv path holding empty list, which will be filled 
                        #with post processed csv results
                        GNPScsvfiles2 = []
                        #print(entry)
                        # enter the directory with /spectral_dereplication/ results
                        sub_dir = input_dir + "/" + entry + '/spectral_dereplication/GNPS/'
                        if os.path.exists(sub_dir):
                            files = (glob.glob(sub_dir+'/*D.csv'))
                            #print(files)
                            for mz, row in msp.iterrows():
                                #print(msp["id_X"][mz])
                                for fls_g in files:
                                    if msp["id_X"][mz] in fls_g:
                                        gnps_df = pd.read_csv(fls_g)
                                        gnps_df = gnps_df.drop_duplicates(subset=['GNPSspectrumID'])
                                        if len(gnps_df)>0:

                                            for i, row in gnps_df.iterrows():
                                                # if compound name is present
                                                if GNPS_Scoring(gnps_df, i):
                                                    if not isNaN(gnps_df['GNPScompound_name'][i]):
                                                        # split if there is a gap in the names
                                                        string_chng = (gnps_df['GNPScompound_name'][i].split(" "))
                                                        # create an empty list
                                                        newstr = []

                                                        # for each part of the string in the names
                                                        chng = []

                                                        for j in range(len(string_chng)):
                                                            # check if the substrings are present in the matches and no - is present
                                                            if not any(x in string_chng[j] for x in matches): #and not '-' == string_chng[j]:
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
                                                else:
                                                    gnps_df.drop([i], axis=0, inplace=True)

                                            for k, row in gnps_df.iterrows():
                                                if isNaN(gnps_df['GNPSSMILES'][k]):
                                                    if "[" in gnps_df['GNPScompound_name'][k].split(" ")[-1]:
                                                        string_chng = (gnps_df['GNPScompound_name'][k].split("["))
                                                        #print(gnps_df['GNPScompound_name'][i])
                                                        keep_names = []
                                                        for j in range(len(string_chng)-1):
                                                            gnps_df.loc[k, "corr_names"] == string_chng[j]
                                                            s = pcp.get_compounds(string_chng[j], 'name')

                                                            if s:
                                                                for comp in s:
                                                                    gnps_df["GNPSSMILES"][k] = comp.isomeric_smiles
                                                                    gnps_df.loc[k, "GNPSformula"] = comp.molecular_formula
                                                                    gnps_df.loc[k, "GNPSinchi"] = Chem.MolToInchi(Chem.MolFromSmiles(comp.isomeric_smiles))
                                                                    
                                                                    
                                                            else:
                                                                gnps_df["GNPSSMILES"][k] = ''
                                                                gnps_df.loc[k, "GNPSformula"] = ''
                                                                gnps_df.loc[k, "GNPSinchi"] = ''
                                                if not isNaN(gnps_df['GNPSSMILES'][k]):
                                                    try:
                                                        
                                                        sx = pcp.get_compounds(gnps_df['GNPSSMILES'][k], 'smiles')
                                                        gnps_df.loc[k, "GNPSinchi"] = Chem.MolToInchi(Chem.MolFromSmiles(comp.isomeric_smiles))
                                                        if sx:
                                                            sx = str(sx)
                                                            comp = pcp.Compound.from_cid([int(x) for x in re.findall(r'\b\d+\b', sx)])
                                                            gnps_df.loc[k, 'GNPSformula'] = comp.molecular_formula
                                                            
                                                    except:
                                                        gnps_df.loc[k, "GNPSformula"] = ''
                                                        gnps_df.loc[k, "GNPSinchi"] = ''



                                        csvname = (os.path.splitext(fls_g)[0])+"proc"+".csv"
                                        gnps_results_csv = csvname.replace(input_dir, ".")
                                        msp.loc[mz, "gnps_results_csv"] = gnps_results_csv 
                                        gnps_df.to_csv(csvname)
                                        GNPScsvfiles2.append(csvname)
                                    dict1 = {'GNPSr': GNPScsvfiles2} 
                                    df = pd.DataFrame(dict1)
                                    #return(df)
                                    
                    msp.to_csv(msp_file[0])                
                                    
                    if Source == "hmdb" or Source == "all":                 

                        if not os.path.exists(input_dir+"/structures.sdf"):
                            #download SDF structures
                            os.system("wget -P " + input_dir + " https://hmdb.ca/system/downloads/current/structures.zip")
                            os.system("unzip "+ input_dir + "structures.zip" + " -d " + input_dir)
                        # Load the sdf
                        dframe = PandasTools.LoadSDF((input_dir+"/structures.sdf"),
                                                     idName='HMDB_ID',smilesName='SMILES',
                                                     molColName='Molecule', includeFingerprints=False)


                        HMDBcsvfiles2 = []
                        #print(entry)
                        # enter the directory with /spectral_dereplication/ results
                        sub_dir = input_dir + "/" + entry + '/spectral_dereplication/HMDB/'


                        if os.path.exists(sub_dir):

                            #print(sub_dir)
                            files = (glob.glob(sub_dir+'/*D.csv'))
                            #print(files)
                            for mz, row in msp.iterrows():
                                #print(msp["id_X"][mz])
                                for fls_h in files:
                                    if msp["id_X"][mz] in fls_h:
                                        hmdb_df = pd.read_csv(fls_h)
                                        hmdb_df = hmdb_df.drop_duplicates(subset=['HMDBcompoundID'])

                                        if len(hmdb_df)>0:
                                            print(entry)
                                            # merge on basis of id, frame and hmdb result files
                                            SmilesHM = pd.merge(hmdb_df, dframe, left_on=hmdb_df.HMDBcompoundID, right_on=dframe.DATABASE_ID)


                                            for i, row in hmdb_df.iterrows():
                                                if HMDB_Scoring(hmdb_df, i):

                                                    for j, row in SmilesHM.iterrows():

                                                        # where index for both match, add the name and SMILES
                                                        if hmdb_df['HMDBcompoundID'][i]== SmilesHM['HMDBcompoundID'][j]:
                                                            hmdb_df.loc[i, 'HMDBSMILES'] = SmilesHM['SMILES'][j]#add SMILES
                                                            hmdb_df.loc[i, 'HMDBcompound_name'] = SmilesHM["GENERIC_NAME"][j]#add name
                                                            hmdb_df.loc[i, 'HMDBformula'] = SmilesHM["FORMULA"][j]#add formula
                                                            hmdb_df.loc[i, 'HMDBinchi'] = Chem.MolToInchi(Chem.MolFromSmiles(SmilesHM['SMILES'][j]))
                                                else:
                                                    hmdb_df.drop([i], axis=0, inplace=True)

                                        csvname = (os.path.splitext(fls_h)[0])+"proc"+".csv" # name for writing it in a new file
                                        hmdb_results_csv = csvname.replace(input_dir, ".")
                                        msp.loc[mz, "hmdb_results_csv"] = hmdb_results_csv 
                                        hmdb_df.to_csv(csvname) #write
                                        HMDBcsvfiles2.append(csvname)# add to a list
                                    dict1 = {'HMDBr': HMDBcsvfiles2} 
                                    df = pd.DataFrame(dict1)
                                    #return(df)
                    
                    msp.to_csv(msp_file[0])
                    # enter the directory with /spectral_dereplication/ results
                    if Source == "mbank" or Source == "all":
                        #open another csv path holding empty list, which will be filled 
                        #with post processed csv results
                        MassBankcsvfiles2 = []
                        #print(entry)
                        # enter the directory with /spectral_dereplication/ results
                        sub_dir = input_dir + "/" + entry + '/spectral_dereplication/MassBank/'
                        if os.path.exists(sub_dir):
                            files = (glob.glob(sub_dir+'/*D.csv'))
                            #print(files)
                            for mz, row in msp.iterrows():
                                #print(msp["id_X"][mz])
                                for fls_m in files:
                                    if msp["id_X"][mz] in fls_m:
                                        print(fls_m)
                                        mbank_df = pd.read_csv(fls_m)
                                        mbank_df = mbank_df.drop_duplicates(subset=['MBspectrumID'])
                                        if len(mbank_df)>0:


                                            for i, row in mbank_df.iterrows():
                                                if MB_Scoring(mbank_df, i):

                                                    inchiK = str(mbank_df["MBinchiKEY"][i])

                                                    #extract inchikeys
                                                    y = pcp.get_compounds(inchiK, 'inchikey')#compound based on inchikey

                                                    for compound in y:

                                                        #add smiles
                                                        smles = compound.isomeric_smiles   
                                                        mbank_df.loc[i, 'MBSMILES'] = smles
                                                        mbank_df.loc[i, 'MBinchi'] =Chem.MolToInchi(Chem.MolFromSmiles(smles))
                                                else:
                                                    mbank_df.drop([i], axis=0, inplace=True)


                                        csvname = (os.path.splitext(fls_m)[0])+"proc"+".csv"
                                        mbank_results_csv = csvname.replace(input_dir, ".")
                                        msp.loc[mz, "mbank_results_csv"] = mbank_results_csv 
                                        mbank_df.to_csv(csvname)
                                        MassBankcsvfiles2.append(csvname)

                                    dict1 = {'MBr': MassBankcsvfiles2} 
                                    df = pd.DataFrame(dict1)
                                    #return(df)                
                    
                    msp.to_csv(msp_file[0])

    if Source == "all":

        dict1 = {'GNPSr': GNPScsvfiles2, 'HMDBr': HMDBcsvfiles2, 'MBr': MassBankcsvfiles2} 
        df = pd.DataFrame(dict1)

        return(df)


# In[9]:


#spec_postproc(input_dir, Source = "mbank")


# # SIRIUS Post Processing

# In[10]:


def sirius_postproc(input_dir, exp_int = 0.90, csi_score = -150):
    def isNaN(string):
        return string != string
    def str_can_score(db, i):
        if db['explainedIntensity'][i] >= exp_int and db['CSI:FingerIDScore'][i] >= csi_score:
            return True
        else:
            return False
    
    # entry is all files and folders in input_dir
    for entry in os.listdir(input_dir):
        #if the entry is also a directory
        if os.path.isdir(os.path.join(input_dir, entry)):
            sub_dir = input_dir + "/" + entry + '/insilico/SIRIUS/'
            msp_csv = input_dir + "/" + entry + "/insilico/MS1DATA.csv"
            if os.path.exists(msp_csv) and os.path.exists(sub_dir):
                #output json files from SIRIUS
                files_S = (glob.glob(sub_dir+'/*.json'))
                #list of precursor m/z
                msp = pd.read_csv(msp_csv)

                # for each mz
                for mz, row in msp.iterrows():
                    # make a list of files with this mz
                    files_for_mz = []

                    for file in files_S:
                        if str(msp["premz"][mz]) in file:
                            files_for_mz.append(file)

                    # in case if SL and ALL are given
                    if len(files_for_mz)==2:

                        # if the SL file is before the ALL file
                        if len(files_for_mz[0])>len(files_for_mz[1]):

                            # extract the formula and structure files
                            json_dirSL = next(os.walk(files_for_mz[0]))[1]
                            sub_sub_dirSL_structure_can = files_for_mz[0] + "/" + json_dirSL[0]  + "/structure_candidates.tsv"                   
                            sub_sub_dirSL_formula_can = files_for_mz[0] + "/" + json_dirSL[0]  + "/formula_candidates.tsv" 
                            SL_Canopus_csv = files_for_mz[0] + "/canopus_summary.tsv"


                            # extract the formula and structure files
                            json_dirALL = next(os.walk(files_for_mz[1]))[1]
                            sub_sub_dirALL_structure_can = files_for_mz[1] + "/" + json_dirALL[0] +"/structure_candidates.tsv"
                            sub_sub_dirALL_formula_can = files_for_mz[1] + "/" + json_dirALL[0] +"/formula_candidates.tsv"
                            ALL_Canopus_csv = files_for_mz[1] + "/canopus_summary.tsv"


                        # if the ALL file is before the SL file
                        elif len(files_for_mz[1]) >len(files_for_mz[0]):


                            # extract the formula and structure files
                            json_dirALL = next(os.walk(files_for_mz[0]))[1]
                            sub_sub_dirALL_structure_can = files_for_mz[0] + "/" + json_dirALL[0]  + "/structure_candidates.tsv"                   
                            sub_sub_dirALL_formula_can = files_for_mz[0] + "/" + json_dirALL[0]  + "/formula_candidates.tsv"
                            ALL_Canopus_csv = files_for_mz[0] + "/canopus_summary.tsv"


                            # extract the formula and structure files
                            json_dirSL = next(os.walk(files_for_mz[1]))[1]
                            sub_sub_dirSL_structure_can = files_for_mz[1] + "/" + json_dirSL[0] +"/structure_candidates.tsv"
                            sub_sub_dirSL_formula_can = files_for_mz[1] + "/" + json_dirALL[0]  + "/formula_candidates.tsv"
                            SL_Canopus_csv = files_for_mz[1] + "/canopus_summary.tsv"

                        # if both structure files exist and they have more than 0 rows
                        if os.path.exists(sub_sub_dirSL_structure_can) and len(pd.read_csv(sub_sub_dirSL_structure_can, sep = "\t"))>0 and os.path.exists(sub_sub_dirALL_structure_can) and len(pd.read_csv(sub_sub_dirALL_structure_can, sep = "\t")):

                            # read strcuture and formula tsv files for both SL and ALL
                            SL_structure_csv = pd.read_csv(sub_sub_dirSL_structure_can, sep = "\t")
                            SL_formula_csv = pd.read_csv(sub_sub_dirSL_formula_can, sep = "\t")
                            SL_Canopus = pd.read_csv(SL_Canopus_csv, sep = "\t")


                            ALL_structure_csv = pd.read_csv(sub_sub_dirALL_structure_can, sep = "\t")
                            ALL_formula_csv = pd.read_csv(sub_sub_dirALL_formula_can, sep = "\t")
                            ALL_Canopus = pd.read_csv(ALL_Canopus_csv, sep = "\t")

                            # Add the structure and formula files together
                            for structure, rows in ALL_structure_csv.iterrows():
                                for formula, rows in ALL_formula_csv.iterrows():
                                    if ALL_structure_csv["formulaRank"][structure] == ALL_formula_csv["rank"][formula]:
                                        ALL_structure_csv.loc[structure, 'SiriusScore'] = ALL_formula_csv['SiriusScore'][formula]
                                        ALL_structure_csv.loc[structure, 'numExplainedPeaks'] = ALL_formula_csv['numExplainedPeaks'][formula]
                                        ALL_structure_csv.loc[structure, 'explainedIntensity'] = ALL_formula_csv['explainedIntensity'][formula]
                                        ALL_structure_csv.loc[structure, "SuspectListEntry"] = "FALSE"
                                        if ALL_formula_csv["molecularFormula"][formula] == ALL_Canopus["molecularFormula"][0]:
                                            ALL_structure_csv.loc[structure, 'superclass'] = ALL_Canopus['superclass'][0]
                                            ALL_structure_csv.loc[structure, 'class'] = ALL_Canopus['class'][0]
                                            ALL_structure_csv.loc[structure, 'subclass'] = ALL_Canopus['subclass'][0]

                            # Add the structure and formula files together
                            for structure_sl, rows in SL_structure_csv.iterrows():
                                for formula_sl, rows in SL_formula_csv.iterrows():
                                    if SL_structure_csv["formulaRank"][structure_sl] == SL_formula_csv["rank"][formula_sl]:
                                        SL_structure_csv.loc[structure_sl, 'SiriusScore'] = SL_formula_csv['SiriusScore'][formula_sl]
                                        SL_structure_csv.loc[structure_sl, 'numExplainedPeaks'] = SL_formula_csv['numExplainedPeaks'][formula_sl]
                                        SL_structure_csv.loc[structure_sl, 'explainedIntensity'] = SL_formula_csv['explainedIntensity'][formula_sl]
                                        SL_structure_csv.loc[structure_sl, "SuspectListEntry"] = "TRUE"   
                                        if SL_formula_csv["molecularFormula"][formula_sl] == SL_Canopus["molecularFormula"][0]:
                                            SL_structure_csv.loc[structure_sl, 'superclass'] = SL_Canopus['superclass'][0]
                                            SL_structure_csv.loc[structure_sl, 'class'] = SL_Canopus['class'][0]
                                            SL_structure_csv.loc[structure_sl, 'subclass'] = SL_Canopus['subclass'][0] 


                            # after formula and structure have been merged, merge SL and ALL results
                            all_sl_db = pd.concat([ALL_structure_csv, SL_structure_csv], ignore_index=True)
                            for str_sirius, row in all_sl_db.iterrows():
                                if not str_can_score(all_sl_db, str_sirius):
                                    all_sl_db = all_sl_db.drop(str_sirius, inplace=False)

                            if not os.path.exists(sub_dir+"results_for_"+json_dirALL[0].split("_")[-1]):
                                os.mkdir(sub_dir+"results_for_"+json_dirALL[0].split("_")[-1])

                            result_sirius_name = (sub_dir+"results_for_"+json_dirALL[0].split("_")[-1]+"_"+"structure_"+json_dirALL[0].split("_")[-1] + ".csv")
                            msp.loc[mz, "sirius_result_dir"] = result_sirius_name.replace(input_dir, ".")
                            all_sl_db.to_csv(sub_dir+"results_for_"+json_dirALL[0].split("_")[-1]+"_"+"structure_"+json_dirALL[0].split("_")[-1] + ".csv")

                        # if only ALL structure file exists and they have more than 0 rows
                        elif not (os.path.exists(sub_sub_dirSL_structure_can) and len(pd.read_csv(sub_sub_dirSL_structure_can, sep = "\t"))>0) and os.path.exists(sub_sub_dirALL_structure_can) and len(pd.read_csv(sub_sub_dirALL_structure_can, sep = "\t")):
                            ALL_structure_csv = pd.read_csv(sub_sub_dirALL_structure_can, sep = "\t")
                            ALL_formula_csv = pd.read_csv(sub_sub_dirALL_formula_can, sep = "\t")
                            ALL_Canopus = pd.read_csv(ALL_Canopus_csv, sep = "\t")

                            # Add the structure and formula files together
                            for structure, rows in ALL_structure_csv.iterrows():
                                for formula, rows in ALL_formula_csv.iterrows():
                                    if ALL_structure_csv["formulaRank"][structure] == ALL_formula_csv["rank"][formula]:
                                        ALL_structure_csv.loc[structure, 'SiriusScore'] = ALL_formula_csv['SiriusScore'][formula]
                                        ALL_structure_csv.loc[structure, 'numExplainedPeaks'] = ALL_formula_csv['numExplainedPeaks'][formula]
                                        ALL_structure_csv.loc[structure, 'explainedIntensity'] = ALL_formula_csv['explainedIntensity'][formula]
                                        ALL_structure_csv.loc[structure, "SuspectListEntry"] = "FALSE"
                                        if ALL_formula_csv["molecularFormula"][formula] == ALL_Canopus["molecularFormula"][0]:
                                            ALL_structure_csv.loc[structure, 'superclass'] = ALL_Canopus['superclass'][0]
                                            ALL_structure_csv.loc[structure, 'class'] = ALL_Canopus['class'][0]
                                            ALL_structure_csv.loc[structure, 'subclass'] = ALL_Canopus['subclass'][0]

                            for str_siriusA, row in ALL_structure_csv.iterrows():
                                if not str_can_score(ALL_structure_csv, str_siriusA):
                                    ALL_structure_csv = ALL_structure_csv.drop(str_siriusA, inplace=False)


                            if not os.path.exists(sub_dir+"results_for_"+json_dirALL[0].split("_")[-1]):
                                os.mkdir(sub_dir+"results_for_"+json_dirALL[0].split("_")[-1])

                            result_sirius_name = (sub_dir+"results_for_"+json_dirALL[0].split("_")[-1]+"_"+"structure_"+json_dirALL[0].split("_")[-1] + ".csv")
                            msp.loc[mz, "sirius_result_dir"] = result_sirius_name.replace(input_dir, ".")
                            ALL_structure_csv.to_csv(sub_dir+"results_for_"+json_dirALL[0].split("_")[-1]+"_"+"structure_"+json_dirALL[0].split("_")[-1] + ".csv")

                        elif not(os.path.exists(sub_sub_dirSL_structure_can) and len(pd.read_csv(sub_sub_dirSL_structure_can, sep = "\t"))>0 and os.path.exists(sub_sub_dirALL_structure_can) and len(pd.read_csv(sub_sub_dirALL_structure_can, sep = "\t"))):
                            if os.path.exists(sub_sub_dirALL_formula_can) and pd.read_csv(sub_sub_dirALL_formula_can, sep = "\t"):
                                ALL_formula_csv = pd.read_csv(sub_sub_dirALL_formula_can, sep = "\t")
                                ALL_Canopus = pd.read_csv(ALL_Canopus_csv, sep = "\t")
                                for formula, rows in ALL_formula_csv.iterrows():
                                    ALL_formula_csv.loc[formula, 'superclass'] = ALL_Canopus['superclass'][0]
                                    ALL_formula_csv.loc[formula, 'class'] = ALL_Canopus['class'][0]
                                    ALL_formula_csv.loc[formula, 'subclass'] = ALL_Canopus['subclass'][0]

                                for for_siriusA, row in ALL_formula_csv.iterrows():
                                    if not ALL_formula_csv['explainedIntensity'][for_siriusA] >= exp_int:
                                        ALL_formula_csv = ALL_formula_csv.drop(for_siriusA, inplace=False)
                                if not os.path.exists(sub_dir+"results_for_"+json_dirALL[0].split("_")[-1]):
                                    os.mkdir(sub_dir+"results_for_"+json_dirALL[0].split("_")[-1])

                                result_sirius_name = (sub_dir+"results_for_"+json_dirALL[0].split("_")[-1]+"_"+"formula_"+json_dirALL[0].split("_")[-1] + ".csv")
                                msp.loc[mz, "sirius_result_dir"] = result_sirius_name.replace(input_dir, ".")

                                ALL_formula_csv.to_csv(sub_dir+"results_for_"+json_dirALL[0].split("_")[-1]+"_"+"formula_"+json_dirALL[0].split("_")[-1] + ".csv")

                            else:
                                pass
                        else:
                            pass

                    elif len(files_for_mz) == 1:
                        # extract the formula and structure files
                        json_dirALL = next(os.walk(files_for_mz[0]))[1]
                        sub_sub_dirALL_structure_can = files_for_mz[0] + "/" + json_dirALL[0]  + "/structure_candidates.tsv"                   
                        sub_sub_dirALL_formula_can = files_for_mz[0] + "/" + json_dirALL[0]  + "/formula_candidates.tsv" 
                        ALL_Canopus_csv = files_for_mz[1] + "/canopus_summary.tsv"

                        # if both structure files exist
                        if os.path.exists(sub_sub_dirALL_structure_can) and len(pd.read_csv(sub_sub_dirALL_structure_can, sep = "\t"))>0:
                            ALL_structure_csv = pd.read_csv(sub_sub_dirALL_structure_can, sep = "\t")
                            ALL_formula_csv = pd.read_csv(sub_sub_dirALL_formula_can, sep = "\t")
                            ALL_Canopus = pd.read_csv(ALL_Canopus_csv, sep = "\t")
                            # Add the structure and formula files together
                            for structure, rows in ALL_structure_csv.iterrows():
                                for formula, rows in ALL_formula_csv.iterrows():
                                    if ALL_structure_csv["formulaRank"][structure] == ALL_formula_csv["rank"][formula]:
                                        ALL_structure_csv.loc[structure, 'SiriusScore'] = ALL_formula_csv['SiriusScore'][formula]
                                        ALL_structure_csv.loc[structure, 'numExplainedPeaks'] = ALL_formula_csv['numExplainedPeaks'][formula]
                                        ALL_structure_csv.loc[structure, 'explainedIntensity'] = ALL_formula_csv['explainedIntensity'][formula]
                                        ALL_structure_csv.loc[structure, "SuspectListEntry"] = "FALSE"

                                        if ALL_formula_csv["molecularFormula"][formula] == ALL_Canopus["molecularFormula"][0]:
                                            ALL_structure_csv.loc[structure, 'superclass'] = ALL_Canopus['superclass'][0]
                                            ALL_structure_csv.loc[structure, 'class'] = ALL_Canopus['class'][0]
                                            ALL_structure_csv.loc[structure, 'subclass'] = ALL_Canopus['subclass'][0]

                            for str_siriusA, row in ALL_structure_csv.iterrows():
                                if not str_can_score(ALL_structure_csv, str_siriusA):
                                    ALL_structure_csv = ALL_structure_csv.drop(str_siriusA, inplace=False)
                            if not os.path.exists(sub_dir+"results_for_"+json_dirALL[0].split("_")[-1]):
                                os.mkdir(sub_dir+"results_for_"+json_dirALL[0].split("_")[-1])

                            result_sirius_name = (sub_dir+"results_for_"+json_dirALL[0].split("_")[-1]+"_"+"structure_"+json_dirALL[0].split("_")[-1] + ".csv")
                            msp.loc[mz, "sirius_result_dir"] = result_sirius_name.replace(input_dir, ".")

                            ALL_structure_csv.to_csv(sub_dir+"results_for_"+json_dirALL[0].split("_")[-1]+"_"+"structure_"+json_dirALL[0].split("_")[-1] + ".csv")


                        elif not (os.path.exists(sub_sub_dirALL_structure_can) and len(pd.read_csv(sub_sub_dirALL_structure_can, sep = "\t"))):
                            if os.path.exists(sub_sub_dirALL_formula_can) and pd.read_csv(sub_sub_dirALL_formula_can, sep = "\t"):
                                ALL_formula_csv = pd.read_csv(sub_sub_dirALL_formula_can, sep = "\t")
                                ALL_Canopus = pd.read_csv(ALL_Canopus_csv, sep = "\t")
                                for formula, rows in ALL_formula_csv.iterrows():
                                    ALL_formula_csv.loc[formula, 'superclass'] = ALL_Canopus['superclass'][0]
                                    ALL_formula_csv.loc[formula, 'class'] = ALL_Canopus['class'][0]
                                    ALL_formula_csv.loc[formula, 'subclass'] = ALL_Canopus['subclass'][0]

                                for for_siriusA, row in ALL_formula_csv.iterrows():
                                    if not ALL_formula_csv['explainedIntensity'][for_siriusA] >= exp_int:
                                        ALL_formula_csv = ALL_formula_csv.drop(for_siriusA, inplace=False)
                                if not os.path.exists(sub_dir+"results_for_"+json_dirALL[0].split("_")[-1]):
                                    os.mkdir(sub_dir+"results_for_"+json_dirALL[0].split("_")[-1])

                                result_sirius_name = (sub_dir+"results_for_"+json_dirALL[0].split("_")[-1]+"_"+"formula_"+json_dirALL[0].split("_")[-1] + ".csv")
                                msp.loc[mz, "sirius_result_dir"] = result_sirius_name.replace(input_dir, ".")

                                ALL_formula_csv.to_csv(sub_dir+"results_for_"+json_dirALL[0].split("_")[-1]+"_"+"formula_"+json_dirALL[0].split("_")[-1] + ".csv")

                            else:
                                pass
                        else:
                            pass
                msp.to_csv(msp_csv)
                


# In[ ]:





# In[11]:


#sirius_postproc(input_dir, exp_int = 0.90, csi_score = -150)


# # MCSS for SpecDBs

# In[12]:


def MCSS_for_SpecDB(input_dir, Source):
    def isNaN(string):
        return string != string
    # Describe the heavy atoms to be considered for MCSS
    heavy_atoms = ['C', 'N', 'P', 'O', 'S']
    #list all files and directories
    for entry in os.listdir(input_dir):

        if os.path.isdir(os.path.join(input_dir, entry)):

            # for specdb
            specdb_msp_file = (glob.glob(input_dir + "/" + entry + '/spectral_dereplication' +'/*.csv'))

            if len(specdb_msp_file) > 0:

                if os.path.exists(specdb_msp_file[0]):

                    spec_msp = pd.read_csv(specdb_msp_file[0])

                    for mz, row in spec_msp.iterrows():

                        if Source == "gnps" or Source == "specdb" or Source == "all":


                            sub_dir = input_dir + "/" + entry + '/spectral_dereplication/GNPS/'
                            if os.path.exists(sub_dir):
                                gnps_files = (glob.glob(sub_dir+'/*proc.csv'))

                                for files in gnps_files:
                                    if spec_msp["id_X"][mz] in files:
                                        gnpsproc = pd.read_csv(files) 


                                        if len(gnpsproc)>0:
                                            G_Smiles = gnpsproc["GNPSSMILES"]
                                            G_Smiles = list(filter(None, G_Smiles))
                                            #print(G_Smiles)
                                            #create empty list of GNPS top smiles
                                            GNPS_Mol = []
                                            # extract only the InChI of the top 5
                                            for j in list(G_Smiles):
                                                if not isNaN(j):
                                                    print(type(j))
                                                    mol2 = Chem.MolFromSmiles(j)
                                                    GNPS_Mol.append(mol2)

                                            if len(GNPS_Mol) >= 2:
                                                res = rdFMCS.FindMCS(GNPS_Mol)
                                                sm_res = Chem.MolToSmiles(Chem.MolFromSmarts(res.smartsString))
                                                # if there are atleast 3 heavy atoms in the MCSS, then add it to the result file
                                                elem = [ele for ele in heavy_atoms if(ele in sm_res)]
                                                if elem and len(sm_res)>=3:
                                                    spec_msp.loc[mz, 'GNPS_MCSSstring'] = res.smartsString
                                                    spec_msp.loc[mz, 'GNPS_MCSS_SMILES'] = Chem.MolToSmiles(Chem.MolFromSmarts(res.smartsString))
                        if Source == "hmdb" or Source == "specdb" or Source == "all":
                            sub_dir = input_dir + "/" + entry + '/spectral_dereplication/HMDB/'
                            if os.path.exists(sub_dir):
                                hmdb_files = (glob.glob(sub_dir+'/*proc.csv'))

                                for files in hmdb_files:
                                    if spec_msp["id_X"][mz] in files:
                                        hmdbproc = pd.read_csv(files) 

                                        if len(hmdbproc)>0:
                                            H_Smiles = hmdbproc["HMDBSMILES"]
                                            H_Smiles = list(filter(None, H_Smiles))

                                            HMDB_Mol = []
                                            # extract only the InChI of the top 5
                                            for j in list(H_Smiles):
                                                if not isNaN(j):

                                                    mol2 = Chem.MolFromSmiles(j)
                                                    HMDB_Mol.append(mol2)

                                            if len(HMDB_Mol) >= 2:
                                                res = rdFMCS.FindMCS(HMDB_Mol)
                                                sm_res = Chem.MolToSmiles(Chem.MolFromSmarts(res.smartsString))
                                                # if there are atleast 3 heavy atoms in the MCSS, then add it to the result file
                                                elem = [ele for ele in heavy_atoms if(ele in sm_res)]
                                                if elem and len(sm_res)>=3:
                                                    spec_msp.loc[mz, 'HMDB_MCSSstring'] = res.smartsString
                                                    spec_msp.loc[mz, 'HMDB_MCSS_SMILES'] = Chem.MolToSmiles(Chem.MolFromSmarts(res.smartsString))
                        if Source == "mbank" or Source == "specdb" or Source == "all":
                            sub_dir = input_dir + "/" + entry + '/spectral_dereplication/MassBank/'
                            if os.path.exists(sub_dir):
                                mbank_files = (glob.glob(sub_dir+'/*proc.csv'))

                                for files in mbank_files:
                                    if spec_msp["id_X"][mz] in files:
                                        mbankproc = pd.read_csv(files) 

                                        if len(mbankproc)>0:
                                            M_Smiles = mbankproc["MBSMILES"]
                                            M_Smiles = list(filter(None, M_Smiles))

                                            MB_Mol = []
                                            # extract only the InChI of the top 5
                                            for j in list(M_Smiles):
                                                if not isNaN(j):

                                                    mol2 = Chem.MolFromSmiles(j)
                                                    MB_Mol.append(mol2)

                                            if len(MB_Mol) >= 2:
                                                res = rdFMCS.FindMCS(MB_Mol)
                                                sm_res = Chem.MolToSmiles(Chem.MolFromSmarts(res.smartsString))
                                                # if there are atleast 3 heavy atoms in the MCSS, then add it to the result file
                                                elem = [ele for ele in heavy_atoms if(ele in sm_res)]
                                                if elem and len(sm_res)>=3:
                                                    spec_msp.loc[mz, 'MB_MCSSstring'] = res.smartsString
                                                    spec_msp.loc[mz, 'MB_MCSS_SMILES'] = Chem.MolToSmiles(Chem.MolFromSmarts(res.smartsString))
                    spec_msp.to_csv(specdb_msp_file[0])
            return(specdb_msp_file)


# In[13]:


#MCSS_for_SpecDB(input_dir, Source = "all")


# # MCSS for SIRIUS

# In[14]:


def MCSS_for_SIRIUS(input_dir):
    def isNaN(string):
        return string != string
    # Describe the heavy atoms to be considered for MCSS
    heavy_atoms = ['C', 'N', 'P', 'O', 'S']
    #list all files and directories
    for entry in os.listdir(input_dir):

        if os.path.isdir(os.path.join(input_dir, entry)):
            #for sirius
            sirius_msp_csv = input_dir + "/" + entry + "/insilico/MS1DATA.csv"
            sub_dir = input_dir + "/" + entry + '/insilico/SIRIUS/'
            if os.path.exists(sirius_msp_csv) and os.path.exists(sub_dir):
                sirius_msp = pd.read_csv(sirius_msp_csv) 
                sirius_files = (glob.glob(sub_dir))
                for mz, row in sirius_msp.iterrows():
                    for sir_file in sirius_files:
                        r = [s for s in os.listdir(sir_file) if "results" in s][0]
                        sirius_f = (glob.glob((sir_file+r)+'/*.csv'))
                        if str(sirius_msp["id_X"][mz].split("_")[1]) in sirius_f[0]:
                            print(sirius_f)
                            if len(sirius_f) == 1:
                                s_f = pd.read_csv(sirius_f[0])
                                if len(s_f)>0 and 'smiles' in s_f.columns.values.tolist():

                                    S_Smiles = s_f["smiles"]
                                    #create empty list of MB top smiles
                                    SIRIUS_Mol = []

                                    # extract only the InChI of the top 5
                                    for j in list(S_Smiles):
                                        mol2 = Chem.MolFromSmiles(j)
                                        SIRIUS_Mol.append(mol2)
                                    if len(SIRIUS_Mol) >= 2:
                                        res = rdFMCS.FindMCS(SIRIUS_Mol)
                                        sm_res = Chem.MolToSmiles(Chem.MolFromSmarts(res.smartsString))
                                        # if there are atleast 3 heavy atoms in the MCSS, then add it to the result file
                                        elem = [ele for ele in heavy_atoms if(ele in sm_res)]
                                        if elem and len(sm_res)>=3:
                                            sirius_msp.loc[mz, 'SIRIUS_MCSSstring'] = res.smartsString
                                            sirius_msp.loc[mz, 'SIRIUS_MCSS_SMILES'] = Chem.MolToSmiles(Chem.MolFromSmarts(res.smartsString))
                sirius_msp.to_csv(sirius_msp_csv)


# In[ ]:





# # Suspect list Screening with Spec DBs

# In[15]:


def SuspectListScreening(input_dir, SuspectListPath, tanimoto, Source):
    def isNaN(string):
        return string != string
    
    SuspectList = pd.read_csv(SuspectListPath)
    
    
    for entry in os.listdir(input_dir):

        if os.path.isdir(os.path.join(input_dir, entry)):
            
            if Source == "gnps" or Source == "specdb" or Source == "all":
                
                sub_dir = input_dir + "/" + entry + '/spectral_dereplication/GNPS/'
                if os.path.exists(sub_dir):
                    gnps_files = (glob.glob(sub_dir+'/*proc.csv'))
                    for file in gnps_files:
                        gnpsproc = pd.read_csv(file) 
                        if len(gnpsproc)>0:
                            for g, row in gnpsproc.iterrows():
                                for s, row in SuspectList.iterrows():
                                    if not isNaN(gnpsproc['GNPSSMILES'][g]) and gnpsproc['GNPSSMILES'][g] != " ":
                                        if not isNaN(SuspectList['SMILES'][s]) and SuspectList['SMILES'][s] != " ":
                                            LHms2 = [Chem.MolFromSmiles(gnpsproc['GNPSSMILES'][g]), Chem.MolFromSmiles(SuspectList['SMILES'][s])]
                                            LHfps2 = [AllChem.GetMorganFingerprintAsBitVect(x2,2, nBits=2048) for x2 in LHms2]
                                            LHtn2 = DataStructs.FingerprintSimilarity(LHfps2[0],LHfps2[1])
                                            if LHtn2 >= tanimoto:
                                                gnpsproc.loc[g, 'SLGsmiles'] = SuspectList['SMILES'][s]
                                                gnpsproc.loc[g, 'SLGname'] = SuspectList['Name'][s]
                                                gnpsproc.loc[g, 'SLGtanimoto'] = LHtn2
                        gnpsproc.to_csv(file)                       
                        return(gnpsproc)
            if Source == "hmdb" or Source == "specdb" or Source == "all":
                
                sub_dir = input_dir + "/" + entry + '/spectral_dereplication/HMDB/'
                if os.path.exists(sub_dir):
                    hmdb_files = (glob.glob(sub_dir+'/*proc.csv'))
                    for file in hmdb_files:

                        hmdbproc = pd.read_csv(file) 
                        if len(hmdbproc)>0:
                            for h, row in hmdbproc.iterrows():
                                for s, row in SuspectList.iterrows():
                                    if not isNaN(hmdbproc['HMDBSMILES'][h]) and hmdbproc['HMDBSMILES'][h] != " ":
                                        if not isNaN(SuspectList['SMILES'][s]) and SuspectList['SMILES'][s] != " ":
                                            LHms2 = [Chem.MolFromSmiles(hmdbproc['HMDBSMILES'][h]), Chem.MolFromSmiles(SuspectList['SMILES'][s])]
                                            LHfps2 = [AllChem.GetMorganFingerprintAsBitVect(x2,2, nBits=2048) for x2 in LHms2]
                                            LHtn2 = DataStructs.FingerprintSimilarity(LHfps2[0],LHfps2[1])
                                            if LHtn2 >= tanimoto:
                                                hmdbproc.loc[h, 'SLHsmiles'] = SuspectList['SMILES'][s]
                                                hmdbproc.loc[h, 'SLHname'] = SuspectList['Name'][s]
                                                hmdbproc.loc[h, 'SLHtanimoto'] = LHtn2


                        hmdbproc.to_csv(file)                       
                        return(hmdbproc)
                        
                        
            if Source == "mbank" or Source == "specdb" or Source == "all":
                
                sub_dir = input_dir + "/" + entry + '/spectral_dereplication/MassBank/'
                if os.path.exists(sub_dir):
                    mbank_files = (glob.glob(sub_dir+'/*proc.csv'))
                    for file in mbank_files:

                        mbankproc = pd.read_csv(file) 

                        if len(mbankproc)>0:
                            for m, row in mbankproc.iterrows():
                                for s, row in SuspectList.iterrows():
                                    if not isNaN(mbankproc['MBSMILES'][m]) and mbankproc['MBSMILES'][m] != " ":
                                        if not isNaN(SuspectList['SMILES'][s]) and SuspectList['SMILES'][s] != " ":
                                            LHms2 = [Chem.MolFromSmiles(mbankproc['MBSMILES'][m]), Chem.MolFromSmiles(SuspectList['SMILES'][s])]
                                            LHfps2 = [AllChem.GetMorganFingerprintAsBitVect(x2,2, nBits=2048) for x2 in LHms2]
                                            LHtn2 = DataStructs.FingerprintSimilarity(LHfps2[0],LHfps2[1])
                                            if LHtn2 >= tanimoto:
                                                mbankproc.loc[m, 'SLMsmiles'] = SuspectList['SMILES'][s]
                                                mbankproc.loc[m, 'SLMname'] = SuspectList['Name'][s]
                                                mbankproc.loc[m, 'SLMtanimoto'] = LHtn2


                        mbankproc.to_csv(file)                       
                        return(mbankproc)
                            
                            
                            


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# ### GNPS, MassBank and HMDB Results post processing

# In[ ]:





# In[16]:


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
    def isNaN(string):
        return string != string

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
                            
    
    if Source == "hmdb" or Source == "all":

        if not os.path.exists(input_dir+"structures.sdf"):
            #download SDF structures
            os.system("wget -P " + input_dir + " https://hmdb.ca/system/downloads/current/structures.zip")
            os.system("unzip "+ input_dir + "structures.zip" + " -d " + input_dir)
            
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
            dict1 = {'HMDBr': HMDBcsvfiles2} 
            df = pd.DataFrame(dict1)
        
    #MassBank CSV Result file pre_processing
    
    if Source == "mbank" or Source == "all":
        
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
            
            dict1 = {'MBr': MassBankcsvfiles2} 
            df = pd.DataFrame(dict1)
    
    # GNPS CSV Result file pre_processing
    if Source == "gnps" or Source == "all":
        #open another csv path holding empty list, which will be filled 
        #with post processed csv results
        GNPScsvfiles2 = []
        #currently only these subsets are removed from the names from GNPS
        matches = ["M+","[M", "M-", "2M", "M*" "20.0", "50.0", "30.0", "40.0", "60.0", "70.0", "eV", "Massbank"
               , "Spectral", "Match", "to", "from", "NIST14", "MoNA", '[IIN-based:',  '[IIN-based', 'on:', 'CCMSLIB00003136269]']

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

                        # check if the substrings are present in the matches and no - is present
                        if not any(x in string_chng[j] for x in matches): #and not '-' == string_chng[j]:

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
                if isNaN(gnps_df['GNPSSMILES'][i]):
                    if "[" in gnps_df['GNPScompound_name'][i].split(" ")[-1]:
                        string_chng = (gnps_df['GNPScompound_name'][i].split("["))
                        #print(gnps_df['GNPScompound_name'][i])
                        keep_names = []
                        for j in range(len(string_chng)-1):
                            gnps_df.loc[i, "corr_names"] == string_chng[j]
                            s = pcp.get_compounds(string_chng[j], 'name')

                            if s:
                                for comp in s:
                                    gnps_df["GNPSSMILES"][i] = comp.isomeric_smiles
                            else:
                                gnps_df["GNPSSMILES"][i] = ''
                if not isNaN(gnps_df['GNPSSMILES'][i]):
                    try:
                        sx = pcp.get_compounds(gnps_df['GNPSSMILES'][i], 'smiles')
                        if sx:
                            sx = str(sx)
                            comp = pcp.Compound.from_cid([int(x) for x in re.findall(r'\b\d+\b', sx)])
                            gnps_df.loc[i, 'GNPSformula'] = comp.molecular_formula
                    except:
                        gnps_df.loc[i, 'GNPSformula'] = ''

            csvname = (os.path.splitext(l)[0])+"proc"+".csv"
            gnps_df.to_csv(csvname)
            GNPScsvfiles2.append(csvname)
            dict1 = {'GNPSr': GNPScsvfiles2} 
            df = pd.DataFrame(dict1)
        

    if Source == "all":
        
        dict1 = {'GNPSr': GNPScsvfiles2, 'HMDBr': HMDBcsvfiles2, 'MBr': MassBankcsvfiles2} 
        df = pd.DataFrame(dict1)

        return(df)


# ### SIRIUS post processing

# ### SIRIUS Result Post Processing

# In[ ]:





# In[17]:


#print(sirius_postProc2.__doc__)


# In[18]:


def sirius_postProc2(input_dir, input_tablecsv):
    
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
    

    Returns:
    csv: a result file with additional columns such as those for suspect
    list if one is used. It also adds columns on MCSS., named as 
    "input_dir/ResultFileName/insilico/SiriusResults.csv"
    
    
    Usage:
    sirius_postProc2(input_dir = "/user/project/", 
    input_table = "/user/project/suspectlist.csv")


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
                    mol = []
                    for j in top_smiles:
                        n = Chem.MolFromSmiles(j)
                        mol.append(n)
                    # list of mol used to calaculate the MCSS
                    res = rdFMCS.FindMCS(mol)
                    sm_res = Chem.MolToSmiles(Chem.MolFromSmarts(res.smartsString))
                    # Check if the MCSS has one of the heavy atoms and whether they are
                    # more than 3
                    elem = [ele for ele in heavy_atoms if(ele in sm_res)]
                    if elem and len(sm_res)>=3:
                        file1.loc[i, 'MCSSstring'] = res.smartsString
                        file1.loc[i, 'MCSS_SMILES'] = Chem.MolToSmiles(Chem.MolFromSmarts(res.smartsString))
                        
                        
            if file1["FormulaRank"][i] == 1.0:
                sep = 'json/'
                strpd = file1["dir"][i].split(sep, 1)[0] +"json/canopus_summary.tsv"
                if os.path.isfile(strpd):

                    canopus = pd.read_csv(strpd, sep='\t')
                    if len(canopus) > 0:
                        #file1.loc[i, 'most_specific_class'] = canopus["most specific class"][0]
                        #file1.loc[i, 'level _5'] = canopus["level 5"][0]
                        file1.loc[i, 'subclass'] = canopus["subclass"][0]
                        file1.loc[i, 'class'] = canopus["class"][0]
                        file1.loc[i, 'superclass'] = canopus["superclass"][0]
                        #file1.loc[i, 'all_classifications'] = canopus["all classifications"][0]
                        file1.loc[i, 'Classification_Source'] = 'CANOPUS'
                    
        
        file1.to_csv(input_dir + (input_table['ResultFileNames'][m] + '/insilico/SiriusResults.csv').replace("./", ""))


# In[19]:


# check for one file
#sirius_postProc2(input_dir= "/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/", 
#                 input_tablecsv = "/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/input_table.csv")


# In[20]:


# check for multiple file
#sirius_postProc2(input_dir= "/Users/mahnoorzulfiqar/Downloads/MAW-main/", 
#                 input_tablecsv = "/Users/mahnoorzulfiqar/Downloads/MAW-main/input_table.csv")


# ### MetFrag Result Post Processing

# In[21]:


def metfrag_postproc(input_dir, input_tablecsv, sl= True):
    
    
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


# In[22]:


# check for multiple file
#metfrag_postproc(input_dir= "/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/", 
#                 input_tablecsv = "/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/input_table.csv")


# In[23]:


# check for multiple file
#metfrag_postproc(input_dir= "/Users/mahnoorzulfiqar/Downloads/MAW-main/", 
#                 input_tablecsv = "/Users/mahnoorzulfiqar/Downloads/MAW-main/input_table.csv")


# In[24]:


#print(metfrag_postproc.__doc__)


# ### COMBINE IN SILICO -All files with SIRIUS results separate and with MetFragresults separate

# In[25]:


def combine_insilico(input_dir, input_tablecsv, Source = "all_insilico"):
    
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
    
    input_table = pd.read_csv(input_tablecsv)
    # create a new directory to store all results /MetabolomicsResults/
    path = os.path.join(input_dir, "MetabolomicsResults")
    if not os.path.isdir(path):
        os.mkdir(path)    
    # if Sirius results are to be combined
    if Source == "all_insilico" or Source == "SIRIUS":
        
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
        frame.to_csv(input_dir + '/MetabolomicsResults/SIRIUS_combined.csv')       
    
    # if MetFrag results are to be combined
    if Source == "all_insilico" or Source == "MetFrag":
        
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
        


# In[26]:


# check for one file
#combine_insilico(input_dir= "/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/", 
 #                input_tablecsv = "/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/input_table.csv")


# In[27]:


# check for multiple file
##combine_insilico(input_dir= "/Users/mahnoorzulfiqar/Downloads/MAW-main/", 
#                 input_tablecsv = "/Users/mahnoorzulfiqar/Downloads/MAW-main/input_table.csv")


# In[28]:


#print(combine_insilico.__doc__)


# In[ ]:





# In[29]:


# check for one file
#spec_postproc(input_dir= "/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/",
#             Source = "gnps")


# In[30]:


# check for multiple file
#spec_postproc(input_dir= "/Users/mahnoorzulfiqar/Downloads/MAW-main/", 
#             Source = "gnps")


# In[31]:


#print(spec_postproc.__doc__)


# ### Combine_all Spectral DBs for one file

# In[32]:


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
    def isNaN(string):
        return string != string

    
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
                    if 'gnpsproc.' in f: 
                        GNPScsvfiles2.append(f)
                    if 'hmdbproc.' in f: 
                        HMDBcsvfiles2.append(f)
                    if 'mbankproc.' in f: 
                        MassBankcsvfiles2.append(f)
   
    # if all results present
    if len(GNPScsvfiles2)>0 and len(HMDBcsvfiles2)>0 and len(MassBankcsvfiles2)>0:
        
        dict1 = {'GNPSr': GNPScsvfiles2, 'HMDBr': HMDBcsvfiles2, 'MBr': MassBankcsvfiles2} 
        df = pd.DataFrame(dict1)
    
        Merged_Result_df = []
        for i, row in df.iterrows():
            CSVfileG = pd.read_csv(df["GNPSr"][i])
            CSVfileH = pd.read_csv(df["HMDBr"][i])
            CSVfileM = pd.read_csv(df["MBr"][i])
            if os.path.exists(df["MBr"][i]) and os.path.exists(df["HMDBr"][i]) and os.path.exists(df["GNPSr"][i]):
                # merge on the basis of Idx
                MergedRE = CSVfileG.merge(CSVfileH,on='id_X').merge(CSVfileM,on='id_X')
                csvname = (df["GNPSr"][i]).replace("gnpsproc", "mergedR")
                MergedRE.to_csv(csvname)
                Merged_Result_df.append(csvname)
                
                
    # if only GNPS and MassBank           
    if len(GNPScsvfiles2)>0 and len(HMDBcsvfiles2)==0 and len(MassBankcsvfiles2)>0:
            dict1 = {'GNPSr': GNPScsvfiles2, 'MBr': MassBankcsvfiles2} 
            df = pd.DataFrame(dict1)
            Merged_Result_df = []
            for i, row in df.iterrows():
                CSVfileG = pd.read_csv(df["GNPSr"][i])
                CSVfileM = pd.read_csv(df["MBr"][i])
                if os.path.exists(df["MBr"][i]) and os.path.exists(df["GNPSr"][i]):
                    # merge on the basis of Idx
                    MergedRE = CSVfileG.merge(CSVfileM,on='id_X')
                    csvname = (df["MBr"][i]).replace("mbankproc", "mergedR")
                    MergedRE.to_csv(csvname)
                    Merged_Result_df.append(csvname)
            
            
            
            
    # if only GNPS and Hmdb
    if not isNaN(GNPScsvfiles2) and not isNaN(HMDBcsvfiles2) and isNaN(MassBankcsvfiles2):
            dict1 = {'GNPSr': GNPScsvfiles2, 'HMDBr': MassBankcsvfiles2} 
            df = pd.DataFrame(dict1)
            Merged_Result_df = []
            for i, row in df.iterrows():
                CSVfileG = pd.read_csv(df["GNPSr"][i])
                CSVfileH = pd.read_csv(df["HMDBr"][i])
                if os.path.exists(df["HMDBr"][i]) and os.path.exists(df["GNPSr"][i]):
                    # merge on the basis of Idx
                    MergedRE = CSVfileG.merge(CSVfileH,on='id_X')
                    csvname = (df["GNPSr"][i]).replace("gnpsproc", "mergedR")
                    MergedRE.to_csv(csvname)
                    Merged_Result_df.append(csvname)
                
                
                
    # if only MBANK and Hmdb
    if not isNaN(GNPScsvfiles2) and isNaN(HMDBcsvfiles2) and isNaN(MassBankcsvfiles2):
            dict1 = {'GNPSr': GNPScsvfiles2, 'HMDBr': MassBankcsvfiles2} 
            df = pd.DataFrame(dict1)   
            dict1 = {'GNPSr': GNPScsvfiles2, 'HMDBr': MassBankcsvfiles2} 
            df = pd.DataFrame(dict1)
            Merged_Result_df = []
            for i, row in df.iterrows():
                CSVfileG = pd.read_csv(df["MBr"][i])
                CSVfileH = pd.read_csv(df["HMDBr"][i])
                if os.path.exists(df["MBr"][i]) and os.path.exists(df["HMDBr"][i]):
                    # merge on the basis of Idx
                    MergedRE = CSVfileM.merge(CSVfileH,on='id_X')
                    csvname = (df["MBr"][i]).replace("mbankproc", "mergedR")
                    MergedRE.to_csv(csvname)
                    Merged_Result_df.append(csvname)


# In[33]:


# check for multiple file
#combine_specdb(input_dir= "/Users/mahnoorzulfiqar/Downloads/MAW-main/")


# In[34]:


# check for one file
#combine_specdb(input_dir= "/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/")


# In[ ]:





# In[35]:


#print(combine_specdb.__doc__)


# ### Combine all files for spectral db dereplication

# In[36]:


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
    def isNaN(string):
        return string != string
    # create a new directory to store all results /MetabolomicsResults/
    path = os.path.join(input_dir, "MetabolomicsResults")
    if not os.path.isdir(path):
        os.mkdir(path)
        
        
    Mergedcsvfiles = []
    single_file = []
    
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
                    else:
                        single_file.append(f)
    
    if len(Mergedcsvfiles)>0:
        combined_csv = pd.concat([pd.read_csv(l) for l in Mergedcsvfiles], ignore_index=True)
        combined_csv.to_csv(input_dir + 'MetabolomicsResults/SD_post_processed_combined_results.csv')
        return(combined_csv)
    else:
        single_csv = pd.read_csv(single_file[0])
        single_csv.to_csv(input_dir + 'MetabolomicsResults/SD_post_processed_combined_results.csv')
        return(single_csv)
    
    #for i, row in combined_csv.iterrows():
        #if combined_csv['GNPSSMILES'][i] == ' ' or isNaN(combined_csv['GNPSSMILES'][i]):
            #combined_csv['GNPSSMILES'][i] = ''
            
    #for i, row in combined_csv.iterrows():
        #if not isNaN(combined_csv['MBinchiKEY'][i]):
            #try:
                #y = pcp.get_compounds(combined_csv['MBinchiKEY'][i], 'inchikey')
                #if len(y)>1:
                    #combined_csv['MBSMILES'][i] = y[0].isomeric_smiles
            #except:
                #pass
            
    


# In[37]:


# check for one file
#combine_allspec(input_dir= "/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/")


# In[38]:


# check for multiple file
#combine_allspec(input_dir= "/Users/mahnoorzulfiqar/Downloads/MAW-main/")


# In[39]:


#print(combine_allspec.__doc__)


# ### Scoring Scheme for Spectral DB Dereplication

# In[ ]:





# In[40]:


def scoring_spec(input_dir, spec_file):
    
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
    def isNaN(string):
        return string != string
    # the scoring highly depends on the following information:
    # similarity scores should be higher than 0.75
    # intScore >=0.50
    # mzScore >= 0.50
    # ratio of the matchingpeaks by the totalpeaks in the query >= 0.50
    
    combined = pd.read_csv(spec_file)
    
    def HMDB_Scoring(db, i):
        if db['HMDBmax_similarity'][i] >= 0.75 and db['HMDBintScore'][i] >= 0.50 and db['HMDBmzScore'][i] >= 0.50 and db['HQMatchingPeaks'][i]/db['hQueryTotalPeaks'][i] >= 0.50:
            return True
        else:
            return False
    
    
    def GNPS_Scoring(db, i):
        if db['GNPSmax_similarity'][i] >= 0.90 and db['GNPSintScore'][i] >= 0.50 and db['GNPSmzScore'][i] >= 0.50 and db['GQMatchingPeaks'][i]/db['gQueryTotalPeaks'][i] >= 0.50:
            return True
        else:
            return False
    
    
    def MB_Scoring(db, i):
        if db['MBmax_similarity'][i] >= 0.50 and db['MBintScore'][i] >= 0.50 and db['MBmzScore'][i] >= 0.50 and db['MQMatchingPeaks'][i]/db['mQueryTotalPeaks'][i] >= 0.50:
            return True
        else:
            return False
    
    
    for i, row in combined.iterrows():
        
        
        if 'HMDBSMILES' in combined.columns and 'MBSMILES' in combined.columns and 'GNPSSMILES' in combined.columns:
            
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
        
        if 'HMDBSMILES' not in combined.columns and 'MBSMILES' in combined.columns and 'GNPSSMILES' in combined.columns:

            # if MassBank and GNPS show good candidates accorindg to the scoring
            if GNPS_Scoring(combined, i) and MB_Scoring(combined, i) and not isNaN(combined['MBSMILES'][i]) and not isNaN(combined['GNPSSMILES'][i]):
                GMms = [Chem.MolFromSmiles(combined['GNPSSMILES'][i]), Chem.MolFromSmiles(combined['MBSMILES'][i])]
                GMfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in GMms]
                GMtn = DataStructs.FingerprintSimilarity(GMfps[0],GMfps[1])
        
                combined.loc[i, 'annotation'] = 'GNPS, MassBank'
                combined.loc[i, 'tanimotoGM'] = GMtn
                combined.loc[i, 'occurence'] = 2
            # only GNPS
            if GNPS_Scoring(combined, i) and not MB_Scoring(combined, i):
        
                combined.loc[i, 'annotation'] = 'GNPS'
                combined.loc[i, 'tanimotoGM'] = np.nan
                combined.loc[i, 'occurence'] = 1
        
            # only MassBank
            if not GNPS_Scoring(combined, i) and MB_Scoring(combined, i):
        
                combined.loc[i, 'annotation'] = 'MassBank'
                combined.loc[i, 'tanimotoGM'] = np.nan
                combined.loc[i, 'occurence'] = 1
                
            # none
            if not GNPS_Scoring(combined, i) and not MB_Scoring(combined, i):
                combined.loc[i, 'annotation'] = 'none'
                combined.loc[i, 'tanimotoGM'] = np.nan
                combined.loc[i, 'occurence'] = 0
                
                
                
        if 'HMDBSMILES' in combined.columns and 'MBSMILES' not in combined.columns and 'GNPSSMILES' in combined.columns:
            # if HMDB and GNPS show good candidates accorindg to the scoring
            if HMDB_Scoring(combined, i) and GNPS_Scoring(combined, i) and not isNaN(combined['GNPSSMILES'][i]) and not isNaN(combined['HMDBSMILES'][i]):
                HGms = [Chem.MolFromSmiles(combined['HMDBSMILES'][i]), Chem.MolFromSmiles(combined['GNPSSMILES'][i])]
                HGfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in HGms]
                HGtn = DataStructs.FingerprintSimilarity(HGfps[0],HGfps[1])
        
                combined.loc[i, 'annotation'] = 'HMDB, GNPS'
                combined.loc[i, 'tanimotoHG'] = HGtn
                combined.loc[i, 'occurence'] = 2
        
            # only HMDB
            if HMDB_Scoring(combined, i) and not GNPS_Scoring(combined, i):
        
                combined.loc[i, 'annotation'] = 'HMDB'
                combined.loc[i, 'tanimotoHG'] = np.nan
                combined.loc[i, 'occurence'] = 1
            
            # only GNPS
            if not HMDB_Scoring(combined, i) and GNPS_Scoring(combined, i):
        
                combined.loc[i, 'annotation'] = 'GNPS'
                combined.loc[i, 'tanimotoHG'] = np.nan
                combined.loc[i, 'occurence'] = 1
            # none
            if not HMDB_Scoring(combined, i) and not GNPS_Scoring(combined, i):
                combined.loc[i, 'annotation'] = 'none'
                combined.loc[i, 'tanimotoHG'] = np.nan
                combined.loc[i, 'occurence'] = 0
    
        if 'HMDBSMILES' in combined.columns and 'MBSMILES' in combined.columns and 'GNPSSMILES' not in combined.columns:
            
            # if MassBank and HMDB show good candidates accorindg to the scoring
            if HMDB_Scoring(combined, i) and MB_Scoring(combined, i) and not isNaN(combined['MBSMILES'][i]) and not isNaN(combined['HMDBSMILES'][i]):
                HMms = [Chem.MolFromSmiles(combined['HMDBSMILES'][i]), Chem.MolFromSmiles(combined['MBSMILES'][i])]
                HMfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in HMms]
                HMtn = DataStructs.FingerprintSimilarity(HMfps[0],HMfps[1])
        
                combined.loc[i, 'annotation'] = 'HMDB, MassBank'
                combined.loc[i, 'tanimotoHM'] = HMtn
                combined.loc[i, 'occurence'] = 2
                
            # only HMDB
            if HMDB_Scoring(combined, i) and not MB_Scoring(combined, i):
        
                combined.loc[i, 'annotation'] = 'HMDB'
                combined.loc[i, 'tanimotoHM'] = np.nan
                combined.loc[i, 'occurence'] = 1
            
            # only MassBank
            if not HMDB_Scoring(combined, i) and MB_Scoring(combined, i):
        
                combined.loc[i, 'annotation'] = 'MassBank'
                combined.loc[i, 'tanimotoHM'] = np.nan
                combined.loc[i, 'occurence'] = 1
        
            # none
            if not HMDB_Scoring(combined, i) and not MB_Scoring(combined, i):
                combined.loc[i, 'annotation'] = 'none'
                combined.loc[i, 'tanimotoHM'] = np.nan
                combined.loc[i, 'occurence'] = 0
        
        
        #If only HMDB was used
        
        if 'HMDBSMILES' in combined.columns and 'MBSMILES' not in combined.columns and 'GNPSSMILES' not in combined.columns:
            # only HMDB
            if HMDB_Scoring(combined, i):
        
                combined.loc[i, 'annotation'] = 'HMDB'
                combined.loc[i, 'occurence'] = 1
            
            # none
            if not HMDB_Scoring(combined, i):
                combined.loc[i, 'annotation'] = 'none'
                combined.loc[i, 'occurence'] = 0
                
                
        #If only MassBank was used      
                
        if 'HMDBSMILES' not in combined.columns and 'MBSMILES' in combined.columns and 'GNPSSMILES' not in combined.columns:
            # only MassBank
            if MB_Scoring(combined, i):
        
                combined.loc[i, 'annotation'] = 'MassBank'
                combined.loc[i, 'occurence'] = 1
            
            # none
            if not MB_Scoring(combined, i):
                combined.loc[i, 'annotation'] = 'none'
                combined.loc[i, 'occurence'] = 0
        
        
        
        #If only GNPS was used
        
        if 'HMDBSMILES' not in combined.columns and 'MBSMILES' not in combined.columns and 'GNPSSMILES' in combined.columns:
            # only GNPS
            if GNPS_Scoring(combined, i):
        
                combined.loc[i, 'annotation'] = 'GNPS'
                combined.loc[i, 'occurence'] = 1
            
            # none
            if not GNPS_Scoring(combined, i):
                combined.loc[i, 'annotation'] = 'none'
                combined.loc[i, 'occurence'] = 0
                
                
    combined.to_csv(input_dir + "MetabolomicsResults/scoredSpecDB.csv")
    return(combined)


# In[41]:


#scoring_spec(input_dir = "/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/",
#             spec_file = "/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/MetabolomicsResults/SD_post_processed_combined_results.csv")


# In[42]:


##scoring_spec(input_dir = "/Users/mahnoorzulfiqar/Downloads/MAW-main/",
 #            spec_file = "/Users/mahnoorzulfiqar/Downloads/MAW-main/MetabolomicsResults/SD_post_processed_combined_results.csv")


# In[43]:


#print(scoring_spec.__doc__)


# ### Suspect List Screening

# In[44]:


def suspectListScreening(input_dir, slistcsv, SpectralDB_Results, db = "all"):
    
    """suspectListScreening runs tanoimoto similarity score to between
    compounds from the results from spectral DBs and suspect list

    Parameters:
    input_dir (str): This is the input directory where all the .mzML 
    files and their respective result directories are stored.
    slistcsv (str): path to suspect list
    SpectralDB_Results (dataframe): dataframe from scoring_spec
    db(str): can be all, gnps, mbank, hmdb, gm, hg, hm
    
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
    Suspect_list = pd.read_csv(slistcsv)
    
    def isNaN(string):
        return string != string
    if db == "hmdb" or db == "hm" or db == "hg" or db == "all":
        
        # add columns to the result from scoring_spec
        # these columns are for high similiarity canidtes between the databases and suspect list
        SpectralDB_Results['HLsmiles'] = np.nan
        SpectralDB_Results['HLname'] = np.nan

        for i, row in SpectralDB_Results.iterrows():
            if not isNaN(SpectralDB_Results['HMDBSMILES'][i]) and SpectralDB_Results['HMDBSMILES'][i] != " ":
                for j, row in Suspect_list.iterrows():
                    LHms2 = [Chem.MolFromSmiles(SpectralDB_Results['HMDBSMILES'][i]), Chem.MolFromSmiles(Suspect_list['SMILES'][j])]
                    LHfps2 = [AllChem.GetMorganFingerprintAsBitVect(x2,2, nBits=2048) for x2 in LHms2]
                    LHtn2 = DataStructs.FingerprintSimilarity(LHfps2[0],LHfps2[1])
                    if LHtn2 >= 0.9:
                        SpectralDB_Results.loc[i, 'HLsmiles'] = Suspect_list['SMILES'][j]
                        SpectralDB_Results.loc[i, 'HLname'] = Suspect_list['Name'][j]

        # add annotations and occurences
        for i, row in SpectralDB_Results.iterrows():
            if not isNaN(SpectralDB_Results['HLname'][i]):
                SpectralDB_Results['occurence'][i] = SpectralDB_Results['occurence'][i] + 1
                if SpectralDB_Results['annotation'][i] == "none":
                    SpectralDB_Results['annotation'][i] = 'Suspect_List'
                else:
                    SpectralDB_Results['annotation'][i] = SpectralDB_Results['annotation'][i] + ', Suspect_List'
    
    if db == "gnps" or db == "gm" or db == "hg" or db == "all":

        # add columns to the result from scoring_spec
        # these columns are for high similiarity canidtes between the databases and suspect list
        SpectralDB_Results['GLsmiles'] = np.nan
        SpectralDB_Results['GLname'] = np.nan


        for i, row in SpectralDB_Results.iterrows():

            if not isNaN(SpectralDB_Results['GNPSSMILES'][i]) and SpectralDB_Results['GNPSSMILES'][i] != " ":
                for k, row in Suspect_list.iterrows():
                    LGms2 = [Chem.MolFromSmiles(SpectralDB_Results['GNPSSMILES'][i]), Chem.MolFromSmiles(Suspect_list['SMILES'][k])]
                    LGfps2 = [AllChem.GetMorganFingerprintAsBitVect(x2,2, nBits=2048) for x2 in LGms2]
                    LGtn2 = DataStructs.FingerprintSimilarity(LGfps2[0],LGfps2[1])
                    if LGtn2 >= 0.9:
                        SpectralDB_Results.loc[i, 'GLsmiles'] = Suspect_list['SMILES'][k]
                        SpectralDB_Results.loc[i, 'GLname'] = Suspect_list['Name'][k]
        # add annotations and occurences
        for i, row in SpectralDB_Results.iterrows():
            if not isNaN(SpectralDB_Results['GLname'][i]):
                SpectralDB_Results['occurence'][i] = SpectralDB_Results['occurence'][i] + 1
                if SpectralDB_Results['annotation'][i] == "none":
                    SpectralDB_Results['annotation'][i] = 'Suspect_List'
                else:
                    SpectralDB_Results['annotation'][i] = SpectralDB_Results['annotation'][i] + ', Suspect_List'
    
    if db == "mbank" or db == "gm" or db == "hm" or db == "all":

        # add columns to the result from scoring_spec
        # these columns are for high similiarity canidtes between the databases and suspect list
        SpectralDB_Results['MLsmiles'] = np.nan
        SpectralDB_Results['MLname'] = np.nan


        for i, row in SpectralDB_Results.iterrows():
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
            if not isNaN(SpectralDB_Results['MLname'][i]):
                SpectralDB_Results['occurence'][i] = SpectralDB_Results['occurence'][i] + 1
                if SpectralDB_Results['annotation'][i] == "none":
                    SpectralDB_Results['annotation'][i] = 'Suspect_List'
                else:
                    SpectralDB_Results['annotation'][i] = SpectralDB_Results['annotation'][i] + ', Suspect_List'
                
    SpectralDB_Results.to_csv(input_dir + "MetabolomicsResults/SpecDBvsSL.csv")
    return(SpectralDB_Results)


# In[45]:


#suspectListScreening(input_dir = "/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/", 
                     #slistcsv = "/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/SkeletonemaSuspectListV1.csv", 
                     #SpectralDB_Results = "/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/MetabolomicsResults/scoredSpecDB.csv", 
                     #db = "all")


# In[46]:


#suspectListScreening(input_dir = "/Users/mahnoorzulfiqar/Downloads/MAW-main/", 
                     #slistcsv = "/Users/mahnoorzulfiqar/Downloads/MAW-main/SkeletonemaSuspectListV1.csv", 
                     ##SpectralDB_Results = "/Users/mahnoorzulfiqar/Downloads/MAW-main/MetabolomicsResults/scoredSpecDB.csv", 
                     #db = "gm")


# In[ ]:





# In[47]:


#print(suspectListScreening.__doc__)


# In[ ]:





# # Final Candidate List Curation

# ## MetFrag Curation

# In[ ]:





# In[48]:


def metfrag_curation(input_dir, metfragcsv, sl = True):
    def isNaN(string):
        return string != string
    
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
        
        
        # If only KEGG
        if not isNaN(metfrag['KG_SMILES'][i]) and isNaN(metfrag['PC_SMILES'][i]):
            metfrag.loc[i, 'Annotation_M'] = 'KEGG'
            if sl:
                if metfrag['KGSL_Score'][i]>=0.9:
                    metfrag.loc[i, 'Annotation_M'] = 'KEGG, SuspectList'
                else:
                    metfrag.loc[i, 'Annotation_M'] = 'KEGG'
    
        # If only Pubchem
        if not isNaN(metfrag['PC_SMILES'][i]) and isNaN(metfrag['KG_SMILES'][i]):
            metfrag.loc[i, 'Annotation_M'] = 'PubChem'
            if sl:
                if metfrag['PCSL_Score'][i]>=0.9:
                    metfrag.loc[i, 'Annotation_M'] = 'PubChem, SuspectList'
                else:
                    metfrag.loc[i, 'Annotation_M'] = 'PubChem'           
        
    
        # If both, calculate the similarity
        if not isNaN(metfrag['PC_SMILES'][i]) and not isNaN(metfrag['KG_SMILES'][i]):
        
            PKms = [Chem.MolFromSmiles(metfrag['KG_SMILES'][i]), Chem.MolFromSmiles(metfrag['PC_SMILES'][i])]
            PKfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in PKms]
            PKtn = DataStructs.FingerprintSimilarity(PKfps[0],PKfps[1])
        
            # if both are similar, add both
            if PKtn == 1:
                metfrag.loc[i, 'Annotation_M'] = 'KEGG, PubChem'
                if sl:
                    if metfrag['KGSL_Score'][i]>=0.9 and metfrag['PCSL_Score'][i]>=0.9:
                        metfrag.loc[i, 'Annotation_M'] = metfrag['Annotation_M'][i] + ", SuspectList"
        
            # if not similar:
            # check Suspect list score and Fragmenter Score
            
            else:
                if not isNaN(metfrag["KG_Score"][i]):
                    metfrag.loc[i, 'Annotation_M'] = 'KEGG'
                else:
                    metfrag.loc[i, 'Annotation_M'] = 'PubChem'
                    
                                
    metfrag.to_csv(input_dir + "MetabolomicsResults/metfrag_curated.csv")  
    return(metfrag)
    


# In[49]:


#print(metfrag_curation.__doc__)


# In[50]:


#metfrag_curation(input_dir = "/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/",
                 ##metfragcsv = "/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/MetabolomicsResults/MetFrag_combined.csv", 
                 #sl = True)


# ## SIRIUS Results Curation

# In[51]:


#metfrag_curation(input_dir = "/Users/mahnoorzulfiqar/Downloads/MAW-main/", 
##                 metfragcsv = "/Users/mahnoorzulfiqar/Downloads/MAW-main/MetabolomicsResults/MetFrag_combined.csv", 
#                 sl = True)


# In[52]:


def sirius_curation(input_dir, siriuscsv, sl = True):
    def isNaN(string):
        return string != string
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
        if sirius['exp_int'][i] >= 0.70 and "SIRIUS_SL" not in sirius['Result'][i]:
            sirius.loc[i, 'Annotation_S'] = 'SIRIUS'
            #sirius.loc[i, 'SMILES_final'] = sirius['SMILES'][i]
        else:
            if sl:
                
                #If the explained intensity is greater than 0.70 and there is an entry from suspect list
                if sirius['exp_int'][i] >= 0.70 and "SIRIUS_SL" in sirius['Result'][i]:
                    sirius.loc[i, 'Annotation_S'] = 'SIRIUS, SuspectList'
                    #sirius.loc[i, 'SMILES_final'] = sirius['SMILES'][i]
    
                # if the intensity is less thna 0.70 but it still is similar to an entry in Suspect list,
                elif sirius['exp_int'][i] < 0.70 and "SIRIUS_SL" in sirius['Result'][i]:
                    sirius.loc[i, 'Annotation_S'] = 'SIRIUS, SuspectList'
                    #sirius.loc[i, 'SMILES_final'] = sirius['SMILES'][i]
        
    sirius.to_csv(input_dir + "MetabolomicsResults/sirius_curated.csv")
    return(sirius)


# In[53]:


#print(sirius_curation.__doc__)


# In[54]:


#sirius_curation(input_dir = "/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/",
 #                siriuscsv = "/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/MetabolomicsResults/SIRIUS_combined.csv", 
 #                sl = True)


# In[55]:


#sirius_curation(input_dir = "/Users/mahnoorzulfiqar/Downloads/MAW-main/", 
              #   siriuscsv = "/Users/mahnoorzulfiqar/Downloads/MAW-main/MetabolomicsResults/SIRIUS_combined.csv", 
              #   sl = True)


# ## combine curated S and M results

# In[56]:


def combineSM(input_dir, metfragcsv, siriuscsv):
    
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
    
    def isNaN(string):
        return string != string
    
    metfrag = pd.read_csv(metfragcsv)
    sirius = pd.read_csv(siriuscsv)
    S_M_CSV = pd.concat([sirius, metfrag], axis = 1, levels = ["id_X"])
    
    for i, rows in S_M_CSV.iterrows():
        # if results has Sirius Structure annotation, and the explained inetnsity is >= 0.70, keep the annotation as is.
        if S_M_CSV["Result"][i] == "SIRIUS_STR" and S_M_CSV['exp_int'][i] >= 0.70:
            S_M_CSV.loc[i, 'Annotation_C'] = S_M_CSV['Annotation_S'][i]
            
            # to add to that annotation
            if not isNaN(S_M_CSV["Annotation_M"][i]):
                # if annotation has PubChem, by default add SIRIUS
                if S_M_CSV["Annotation_M"][i] == "KEGG":
                    SKms = [Chem.MolFromSmiles(S_M_CSV['SMILES'][i]), Chem.MolFromSmiles(S_M_CSV['KG_SMILES'][i])]
                    SKfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in SKms]
                    SKtn = DataStructs.FingerprintSimilarity(SKfps[0],SKfps[1])

                    if SKtn >= 0.75:

                        S_M_CSV.loc[i, 'Annotation_C'] = S_M_CSV['Annotation_S'][i] +', KEGG'

                    else:
                        S_M_CSV.loc[i, 'Annotation_C'] = S_M_CSV['Annotation_S'][i]
                        
                # if annotation has PubChem, by default add SIRIUS
                if S_M_CSV["Annotation_M"][i] == "PubChem":
                    PSms = [Chem.MolFromSmiles(S_M_CSV['SMILES'][i]), Chem.MolFromSmiles(S_M_CSV['PC_SMILES'][i])]
                    PSfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in PSms]
                    PStn = DataStructs.FingerprintSimilarity(PSfps[0],PSfps[1])

                    # if similar strcutres, then add Pubchme and sirius
                    if PStn >= 0.7:

                        S_M_CSV.loc[i, 'Annotation_C'] = S_M_CSV['Annotation_S'][i] + ', PubChem'

                    # if not then just keep sirius
                    else:
                        S_M_CSV.loc[i, 'Annotation_C'] = S_M_CSV['Annotation_S'][i]
                        
                        
                if S_M_CSV["Annotation_M"][i] == "KEGG, PubChem":
                    SKms = [Chem.MolFromSmiles(S_M_CSV['SMILES'][i]), Chem.MolFromSmiles(S_M_CSV['KG_SMILES'][i])]
                    SKfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in SKms]
                    SKtn = DataStructs.FingerprintSimilarity(SKfps[0],SKfps[1])
                    if SKtn >= 0.7:

                        S_M_CSV.loc[i, 'Annotation_C'] = S_M_CSV['Annotation_S'][i] +', KEGG, PubChem'

                    else:
                        S_M_CSV.loc[i, 'Annotation_C'] = S_M_CSV['Annotation_S'][i]
    S_M_CSV.to_csv(input_dir + "MetabolomicsResults/combinedSM.csv")
    return(S_M_CSV)


# In[57]:


#print(combineSM.__doc__)


# In[58]:


#combineSM(input_dir = "/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/",
#          metfragcsv = "/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/MetabolomicsResults/metfrag_curated.csv", 
 #         siriuscsv = "/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/MetabolomicsResults/sirius_curated.csv")


# In[59]:


#combineSM(input_dir = "/Users/mahnoorzulfiqar/Downloads/MAW-main/",
 #         metfragcsv = "/Users/mahnoorzulfiqar/Downloads/MAW-main/MetabolomicsResults/metfrag_curated.csv", 
 #         siriuscsv = "/Users/mahnoorzulfiqar/Downloads/MAW-main/MetabolomicsResults/sirius_curated.csv")


# ## Spec DB Curation

# In[60]:


def specDB_Curation(input_dir, combinedx, sl = True, db = "all"):
    
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
    def isNaN(string):
        return string != string
    def HMDB_Scoring(db, i):
        if db['HMDBmax_similarity'][i] >= 0.75 and db['HMDBintScore'][i] >= 0.50 and db['HMDBmzScore'][i] >= 0.50 and db['HQMatchingPeaks'][i]/db['hQueryTotalPeaks'][i] >= 0.50:
            return True
        else:
            return False
    
    def GNPS_Scoring(db, i):
        if db['GNPSmax_similarity'][i] >= 0.90 and db['GNPSintScore'][i] >= 0.50 and db['GNPSmzScore'][i] >= 0.50 and db['GQMatchingPeaks'][i]/db['gQueryTotalPeaks'][i] >= 0.50:
            return True
        else:
            return False
    
    
    def MB_Scoring(db, i):
        if db['MBmax_similarity'][i] >= 0.50 and db['MBintScore'][i] >= 0.50 and db['MBmzScore'][i] >= 0.50 and db['MQMatchingPeaks'][i]/db['mQueryTotalPeaks'][i] >= 0.50:
            return True
        else:
            return False
    
    combined = pd.read_csv(combinedx)
    
    
    # remove the similarity scores from low scoring candidates
    for i, row in combined.iterrows():
        if db == "all" or db == "hg" or db == "hm" or db == "hmdb":
            if not HMDB_Scoring(combined, i):
                combined['HMDBcompoundID'][i] = np.nan
        if db == "all" or db == "hg" or db == "gm" or db == "gnps":
            if not GNPS_Scoring(combined, i):
                combined['GNPSspectrumID'][i] = np.nan
        if db == "all" or db == "gm" or db == "hm" or db == "mbank":
            if not MB_Scoring(combined, i):
                combined['MBspectrumID'][i] = np.nan
    
    # if sl = True
    if sl:
        for i, row in combined.iterrows():
            # if all databases are used to generate results
            if db == "all":
                
                # if all dbs have results
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
                        combined.loc[i, 'Annotation'] = 'MassBank, HMDB'
                        # if its present in Suspect List
                        if combined['MLname'][i] == combined['HLname'][i]:
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
                        #elif not isNaN(combined['HLname'][i]):
                            #combined.loc[i, 'Annotation'] = 'HMDB, SuspectList'
            
        
                        # all different annotations, take GNPS
                        else:
                            if not isNaN(combined['GNPSSMILES'][i]):
                                combined.loc[i, 'Annotation'] = 'GNPS'
                            else:
                                combined.loc[i, 'Annotation'] = 'MassBank'
    
                #### When there is an annotation from two DBs #####

                # only GNPS and HMDB
                if not isNaN(combined['GNPSspectrumID'][i]) and isNaN(combined['MBspectrumID'][i]) and not isNaN(combined['HMDBcompoundID'][i]):
                    if combined['tanimotoHG'][i] == 1.0:
                        combined.loc[i, 'Annotation'] = 'GNPS, HMDB'
                        if not isNaN(combined['GLname'][i]) and not isNaN(combined['HLname'][i]):
                            if combined['GLname'][i] == combined['HLname'][i]:
                                combined.loc[i, 'Annotation'] = 'GNPS, HMDB, SuspectList'
                    else:
                        combined.loc[i, 'Annotation'] = 'GNPS'
                        if not isNaN(combined['GLname'][i]):
                            combined.loc[i, 'Annotation'] = 'GNPS, SuspectList'
                        #elif not isNaN(combined['HLname'][i]):
                            #combined.loc[i, 'Annotation'] = 'HMDB, SuspectList'


                # only GNPS and MassBank
                if not isNaN(combined['GNPSspectrumID'][i]) and not isNaN(combined['MBspectrumID'][i]) and isNaN(combined['HMDBcompoundID'][i]):
                    if combined['tanimotoGM'][i] == 1.0:
                        combined.loc[i, 'Annotation'] = 'GNPS, MassBank'
                        if not isNaN(combined['GLname'][i]) and not isNaN(combined['MLname'][i]):
                            if combined['GLname'][i] == combined['MLname'][i]:
                                combined.loc[i, 'Annotation'] = 'GNPS, MassBank, SuspectList'
                                
                    else:
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
                    if combined['tanimotoHM'][i] == 1.0:
                        combined.loc[i, 'Annotation'] = 'HMDB, MassBank'
                        if not isNaN(combined['HLname'][i]) and not isNaN(combined['MLname'][i]):
                            if combined['HLname'][i] == combined['MLname'][i]:
                                combined.loc[i, 'Annotation'] = 'HMDB, MassBank, SuspectList'
                                
                    else:
                        if not isNaN(combined['MLname'][i]):
                            combined.loc[i, 'Annotation'] = 'MassBank, SuspectList'
                        #elif not isNaN(combined['MLname'][i]):
                            #combined.loc[i, 'Annotation'] = 'MassBank, SuspectList'
                        #elif not isNaN(combined['GNPSSMILES'][i]):
                            #combined.loc[i, 'Annotation'] = 'GNPS'
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
                #if isNaN(combined['GNPSspectrumID'][i]) and isNaN(combined['MBspectrumID'][i]) and not isNaN(combined['HMDBcompoundID'][i]):
                    #combined.loc[i, 'Annotation'] = 'HMDB'
                    #If also SuspectList
                    #if not isNaN(combined['HLname'][i]):
                        #combined.loc[i, 'Annotation'] = 'HMDB, SuspectList'
                        
                        
                        
            
            
            # if GNPS AND MassBank databases are used to generate results
            if db == "gm":
                
                # only GNPS and MassBank
                if not isNaN(combined['GNPSspectrumID'][i]) and not isNaN(combined['MBspectrumID'][i]):
                    if combined['tanimotoGM'][i] == 1.0:
                        combined.loc[i, 'Annotation'] = 'GNPS, MassBank'
                        if not isNaN(combined['GLname'][i]) and not isNaN(combined['MLname'][i]):
                            if combined['GLname'][i] == combined['MLname'][i]:
                                combined.loc[i, 'Annotation'] = 'GNPS, MassBank, SuspectList'
                                
                    else:
                        if not isNaN(combined['GLname'][i]):
                            combined.loc[i, 'Annotation'] = 'GNPS, SuspectList'
                        elif not isNaN(combined['MLname'][i]):
                            combined.loc[i, 'Annotation'] = 'MassBank, SuspectList'
                        elif not isNaN(combined['GNPSSMILES'][i]):
                            combined.loc[i, 'Annotation'] = 'GNPS'
                        else:
                            combined.loc[i, 'Annotation'] = 'MassBank'

                
                
                # only GNPS
                if not isNaN(combined['GNPSspectrumID'][i]) and isNaN(combined['MBspectrumID'][i]):

                    #If also SuspectList
                    if not isNaN(combined['GLname'][i]):
                        combined.loc[i, 'Annotation'] = 'GNPS, SuspectList'
                    elif not isNaN(combined['GNPSSMILES'][i]):
                        combined.loc[i, 'Annotation'] = 'GNPS'

                # only MassBank
                if isNaN(combined['GNPSspectrumID'][i]) and not isNaN(combined['MBspectrumID'][i]):
                    combined.loc[i, 'Annotation'] = 'MassBank'
                    #If also SuspectList
                    if not isNaN(combined['MLname'][i]):
                        combined.loc[i, 'Annotation'] = 'MassBank, SuspectList'
            
            
            
            # if GNPS AND HMDB databases are used to generate results
            if db == "hg":
                
                # only GNPS and HMDB
                if not isNaN(combined['GNPSspectrumID'][i]) and not isNaN(combined['HMDBcompoundID'][i]):
                    if combined['tanimotoHG'][i] == 1.0:
                        combined.loc[i, 'Annotation'] = 'GNPS, HMDB'
                        if not isNaN(combined['GLname'][i]) and not isNaN(combined['HLname'][i]):
                            if combined['GLname'][i] == combined['HLname'][i]:
                                combined.loc[i, 'Annotation'] = 'GNPS, HMDB, SuspectList'
                    else:
                        combined.loc[i, 'Annotation'] = 'GNPS'
                        if not isNaN(combined['GLname'][i]):
                            combined.loc[i, 'Annotation'] = 'GNPS, SuspectList'
                        #elif not isNaN(combined['HLname'][i]):
                            #combined.loc[i, 'Annotation'] = 'HMDB, SuspectList'

                # only GNPS
                if not isNaN(combined['GNPSspectrumID'][i]) and isNaN(combined['HMDBcompoundID'][i]):

                    #If also SuspectList
                    if not isNaN(combined['GLname'][i]):
                        combined.loc[i, 'Annotation'] = 'GNPS, SuspectList'
                    elif not isNaN(combined['GNPSSMILES'][i]):
                        combined.loc[i, 'Annotation'] = 'GNPS'
                # only HMDB
                #if isNaN(combined['GNPSspectrumID'][i]) and not isNaN(combined['HMDBcompoundID'][i]):
                    #combined.loc[i, 'Annotation'] = 'HMDB'
                    #If also SuspectList
                    #if not isNaN(combined['HLname'][i]):
                        #combined.loc[i, 'Annotation'] = 'HMDB, SuspectList'
            
            # if MassBank AND HMDB databases are used to generate results
            if db == "hm":
                
                # only MassBank and HMDB

                if not isNaN(combined['MBspectrumID'][i]) and not isNaN(combined['HMDBcompoundID'][i]):
                    if combined['tanimotoHM'][i] == 1.0:
                        combined.loc[i, 'Annotation'] = 'HMDB, MassBank'
                        if not isNaN(combined['HLname'][i]) and not isNaN(combined['MLname'][i]):
                            if combined['HLname'][i] == combined['MLname'][i]:
                                combined.loc[i, 'Annotation'] = 'HMDB, MassBank, SuspectList'
                                
                    else:
                        if not isNaN(combined['MLname'][i]):
                            combined.loc[i, 'Annotation'] = 'MassBank, SuspectList'
                        else:
                            combined.loc[i, 'Annotation'] = 'MassBank'
                
                
                
                # only MassBank
                if not isNaN(combined['MBspectrumID'][i]) and isNaN(combined['HMDBcompoundID'][i]):
                    combined.loc[i, 'Annotation'] = 'MassBank'
                    #If also SuspectList
                    if not isNaN(combined['MLname'][i]):
                        combined.loc[i, 'Annotation'] = 'MassBank, SuspectList'

                # only HMDB
                #if isNaN(combined['MBspectrumID'][i]) and not isNaN(combined['HMDBcompoundID'][i]):
                    #combined.loc[i, 'Annotation'] = 'HMDB'
                    #If also SuspectList
                    #if not isNaN(combined['HLname'][i]):
                        #combined.loc[i, 'Annotation'] = 'HMDB, SuspectList'
            if db == "gnps":
                if not isNaN(combined['GNPSspectrumID'][i]):

                    #If also SuspectList
                    if not isNaN(combined['GLname'][i]):
                        combined.loc[i, 'Annotation'] = 'GNPS, SuspectList'
                    elif not isNaN(combined['GNPSSMILES'][i]):
                        combined.loc[i, 'Annotation'] = 'GNPS'
            if db == "mbank":
                # only MassBank
                if not isNaN(combined['MBspectrumID'][i]):
                    combined.loc[i, 'Annotation'] = 'MassBank'
                    #If also SuspectList
                    if not isNaN(combined['MLname'][i]):
                        combined.loc[i, 'Annotation'] = 'MassBank, SuspectList'
            #if db == "hmdb":
                # only HMDB
                #if not isNaN(combined['HMDBcompoundID'][i]):
                    #combined.loc[i, 'Annotation'] = 'HMDB'
                    #If also SuspectList
                    #if not isNaN(combined['HLname'][i]):
                        #combined.loc[i, 'Annotation'] = 'HMDB, SuspectList'
                        
                        
                        
                        
    else:
        for i, row in combined.iterrows():
            #if all databases were used
            if db == "all":
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
                        combined.loc[i, 'Annotation'] = 'MassBank, HMDB'
                
                    # all different annotations, take GNPS
                    else:
                        if not isNaN(combined['GNPSSMILES'][i]):
                            combined.loc[i, 'Annotation'] = 'GNPS'
                        else:
                            combined.loc[i, 'Annotation'] = 'MassBank'
                ##### When there is an annotation from two DBs #####
    
    
                # only GNPS and HMDB
                if not isNaN(combined['GNPSspectrumID'][i]) and isNaN(combined['MBspectrumID'][i]) and not isNaN(combined['HMDBcompoundID'][i]):
                    if combined['tanimotoHG'][i] == 1.0:
                        combined.loc[i, 'Annotation'] = 'GNPS, HMDB'
                    else:
                        combined.loc[i, 'Annotation'] = 'GNPS'
                    
                    
                # only GNPS and MassBank
                if not isNaN(combined['GNPSspectrumID'][i]) and not isNaN(combined['MBspectrumID'][i]) and isNaN(combined['HMDBcompoundID'][i]):

                    if combined['tanimotoGM'][i] == 1.0:
                        combined.loc[i, 'Annotation'] = 'GNPS, MassBank'
                    else:
                        combined.loc[i, 'Annotation'] = 'GNPS'
    
                # only MassBank and HMDB
                if isNaN(combined['GNPSspectrumID'][i]) and not isNaN(combined['MBspectrumID'][i]) and not isNaN(combined['HMDBcompoundID'][i]):
                    if combined['tanimotoHM'][i] == 1.0:
                        combined.loc[i, 'Annotation'] = 'HMDB, MassBank'
                    else:
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
                    #if isNaN(combined['GNPSspectrumID'][i]) and isNaN(combined['MBspectrumID'][i]) and not isNaN(combined['HMDBcompoundID'][i]):
                    #combined.loc[i, 'Annotation'] = 'HMDB'
                
            
            #if GNPS and MassBank databases were used
            if db == "gm":
                # only GNPS and MassBank
                if not isNaN(combined['GNPSspectrumID'][i]) and not isNaN(combined['MBspectrumID'][i]):
                    if combined['tanimotoGM'][i] == 1.0:
                        combined.loc[i, 'Annotation'] = 'GNPS, MassBank'
                    else:
                        combined.loc[i, 'Annotation'] = 'GNPS'
                    
                
                ##### When there is an annotation from one DBs #####
                # only GNPS
                if not isNaN(combined['GNPSspectrumID'][i]) and isNaN(combined['MBspectrumID'][i]):
                    if not isNaN(combined['GNPSSMILES'][i]):
                        combined.loc[i, 'Annotation'] = 'GNPS'
        
                # only MassBank
                if isNaN(combined['GNPSspectrumID'][i]) and not isNaN(combined['MBspectrumID'][i]):
                    combined.loc[i, 'Annotation'] = 'MassBank'
                    
                    
            # only GNPS and HMDB   
            if db == "hg":
                ##### When there is an annotation from two DBs #####
    
    
                # only GNPS and HMDB
                if not isNaN(combined['GNPSspectrumID'][i]) and not isNaN(combined['HMDBcompoundID'][i]):
                    if combined['tanimotoHG'][i] == 1.0:
                        combined.loc[i, 'Annotation'] = 'GNPS, HMDB'
                    else:
                        combined.loc[i, 'Annotation'] = 'GNPS'
                
                
                ##### When there is an annotation from one DBs #####
    
                # only GNPS
                if not isNaN(combined['GNPSspectrumID'][i]) and isNaN(combined['HMDBcompoundID'][i]):
                    if not isNaN(combined['GNPSSMILES'][i]):
                        combined.loc[i, 'Annotation'] = 'GNPS'
                # only HMDB
                    #if isNaN(combined['GNPSspectrumID'][i]) and not isNaN(combined['HMDBcompoundID'][i]):
                    #combined.loc[i, 'Annotation'] = 'HMDB'
                    
                    
                    
                    
            # only MassBank and HMDB        
            if db == "hm":
                # only MassBank and HMDB
                if not isNaN(combined['MBspectrumID'][i]) and not isNaN(combined['HMDBcompoundID'][i]):
                    if combined['tanimotoHM'][i] == 1.0:
                        combined.loc[i, 'Annotation'] = 'HMDB, MassBank'
                    else:
                        combined.loc[i, 'Annotation'] = 'MassBank'
                
                
                
                # only MassBank
                if not isNaN(combined['MBspectrumID'][i]) and isNaN(combined['HMDBcompoundID'][i]):
                    combined.loc[i, 'Annotation'] = 'MassBank'
    
                # only HMDB
                    #if isNaN(combined['MBspectrumID'][i]) and not isNaN(combined['HMDBcompoundID'][i]):
                    #combined.loc[i, 'Annotation'] = 'HMDB'
                
            if db == "gnps":
                # only GNPS
                if not isNaN(combined['GNPSspectrumID'][i]):
                    if not isNaN(combined['GNPSSMILES'][i]):
                        combined.loc[i, 'Annotation'] = 'GNPS'
            if db == "mbank":
                # only MassBank
                if not isNaN(combined['MBspectrumID'][i]):
                    combined.loc[i, 'Annotation'] = 'MassBank'
            #if db == "hmdb":
                # only HMDB
                #if not isNaN(combined['HMDBcompoundID'][i]):
                    #combined.loc[i, 'Annotation'] = 'HMDB'
    combined.to_csv(input_dir + "MetabolomicsResults/curatedSDB.csv")
    return(combined)
                


# In[61]:


#print(specDB_Curation.__doc__)


# In[62]:


#specDB_Curation(input_dir = "/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/", 
 #               combinedx = "/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/MetabolomicsResults/SpecDBvsSL.csv",
   #             sl = True)


# In[63]:


#specDB_Curation(input_dir = "/Users/mahnoorzulfiqar/Downloads/MAW-main/",
 #               combinedx = "/Users/mahnoorzulfiqar/Downloads/MAW-main/MetabolomicsResults/SpecDBvsSL.csv",
 #               sl = True,
   #            db = "gm")


# # combine curated SDB and CDB (S+M)

# In[64]:


def combine_CuratedR(input_dir, combinedSDBs, combinedSMs, data_type = "standards"):
    
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
    def isNaN(string):
        return string != string

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
                if 'KEGG' in mega['Annotation_C'][i]:
                    if 'MassBank' in mega['Annotation'][i]:
                        SKms = [Chem.MolFromSmiles(mega['MBSMILES'][i]), Chem.MolFromSmiles(mega['KG_SMILES'][i])]
                        SKfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in SKms]
                        SKtn = DataStructs.FingerprintSimilarity(SKfps[0],SKfps[1])
                        if SKtn == 1.0:
                            mega.loc[i, "Annotation_Source"] = mega['Annotation'][i] + ', ' + mega['Annotation_C'][i]
                        else:
                            mega.loc[i, "Annotation_Source"] = mega['Annotation'][i]
            
                    elif 'HMDB' in mega['Annotation'][i]:
                        SKms = [Chem.MolFromSmiles(mega['HMDBSMILES'][i]), Chem.MolFromSmiles(mega['KG_SMILES'][i])]
                        SKfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in SKms]
                        SKtn = DataStructs.FingerprintSimilarity(SKfps[0],SKfps[1])
                        if SKtn == 1.0:
                            mega.loc[i, "Annotation_Source"] = mega['Annotation'][i] + ', ' + mega['Annotation_C'][i]
                        else:
                            mega.loc[i, "Annotation_Source"] = mega['Annotation'][i]
                else:
                    mega.loc[i, "Annotation_Source"] = mega['Annotation'][i]
            
            
            
            #######TWO OR ONE SDB SOURCE#########
                
            #if both 2 SDBs and results from insilico tools
            elif len(mega['Annotation'][i].split()) <= 2:
                mega.loc[i, "Annotation_Source"] = mega['Annotation_C'][i]
                
                
        # if no results from any databases
        if isNaN(mega['Annotation'][i]) and isNaN(mega['Annotation_C'][i]) and not isNaN(mega['Formula'][i]):
            mega.loc[i, "Annotation_Source"] = 'SIRIUS_Formula'
        
    bef_mega = mega.loc[:,~mega.columns.duplicated()]
    for i, row in bef_mega.iterrows():
        if not isNaN(bef_mega['Annotation_Source'][i]):
            # check if SIRIUS is in the annotation source but keep in mind it shouldnt be SIRIUS_Formula
            if 'SIRIUS' in bef_mega['Annotation_Source'][i] and 'SIRIUS_Formula' not in bef_mega['Annotation_Source'][i]:
                bef_mega.loc[i, 'SMILES_final'] = bef_mega['SMILES'][i]
                bef_mega.loc[i,"CompoundNames"] = bef_mega['name'][i]
                bef_mega['PC_MCSS_SMILES'][i] = np.nan
                bef_mega['KG_MCSS_SMILES'][i] = np.nan
            elif 'KEGG' in bef_mega['Annotation_Source'][i]:
                bef_mega.loc[i, 'SMILES_final'] = bef_mega['KG_SMILES'][i]
                bef_mega.loc[i, 'CompoundNames'] = bef_mega['KG_Name'][i]
                #bef_mega['most_specific_class'][i] = np.nan
                #bef_mega['level _5'][i] = np.nan
                bef_mega['subclass'][i] = np.nan
                bef_mega['class'][i] = np.nan
                bef_mega['superclass'][i] = np.nan
                #bef_mega['all_classifications'][i] = np.nan
                bef_mega['Classification_Source'][i] = np.nan
                bef_mega['MCSS_SMILES'][i] = np.nan
                bef_mega['PC_MCSS_SMILES'][i] = np.nan
                bef_mega['Formula'][i] = np.nan
            
            elif 'GNPS, SuspectList' in bef_mega['Annotation_Source'][i]:
                bef_mega.loc[i,'SMILES_final'] = bef_mega['GLsmiles'][i]
                bef_mega.loc[i, 'CompoundNames'] = bef_mega['GLname'][i]
                bef_mega.loc[i, 'CompoundNames']
                #bef_mega['most_specific_class'][i] = np.nan
                #bef_mega['level _5'][i] = np.nan
                bef_mega['subclass'][i] = np.nan
                bef_mega['class'][i] = np.nan
                bef_mega['superclass'][i] = np.nan
                #bef_mega['all_classifications'][i] = np.nan
                bef_mega['MCSS_SMILES'][i] = np.nan
                bef_mega['Classification_Source'][i] = np.nan
                bef_mega['PC_MCSS_SMILES'][i] = np.nan
                bef_mega['KG_MCSS_SMILES'][i] = np.nan
                bef_mega['Formula'][i] = np.nan
        
            elif 'GNPS' in bef_mega['Annotation_Source'][i]:
                bef_mega.loc[i,'SMILES_final'] = bef_mega['GNPSSMILES'][i]
                bef_mega.loc[i, 'CompoundNames'] = bef_mega['GNPScompound_name'][i]
                #bef_mega['most_specific_class'][i] = np.nan
                #bef_mega['level _5'][i] = np.nan
                bef_mega['subclass'][i] = np.nan
                bef_mega['class'][i] = np.nan
                bef_mega['superclass'][i] = np.nan
                #bef_mega['all_classifications'][i] = np.nan
                bef_mega['MCSS_SMILES'][i] = np.nan
                bef_mega['Classification_Source'][i] = np.nan
                bef_mega['PC_MCSS_SMILES'][i] = np.nan
                bef_mega['KG_MCSS_SMILES'][i] = np.nan
                bef_mega['Formula'][i] = np.nan
            elif 'MassBank' in bef_mega['Annotation_Source'][i]:
                bef_mega.loc[i, 'SMILES_final'] = bef_mega['MBSMILES'][i]
                bef_mega.loc[i, 'CompoundNames'] = bef_mega['MBcompound_name'][i]
                #bef_mega['most_specific_class'][i] = np.nan
                #bef_mega['level _5'][i] = np.nan
                bef_mega['subclass'][i] = np.nan
                bef_mega['class'][i] = np.nan
                bef_mega['superclass'][i] = np.nan
                #bef_mega['all_classifications'][i] = np.nan
                bef_mega['Classification_Source'][i] = np.nan
                bef_mega['MCSS_SMILES'][i] = np.nan
                bef_mega['PC_MCSS_SMILES'][i] = np.nan
                bef_mega['KG_MCSS_SMILES'][i] = np.nan
                bef_mega['Formula'][i] = np.nan
                
                
            elif 'PubChem' in bef_mega['Annotation_Source'][i]:
                bef_mega.loc[i, 'SMILES_final'] = bef_mega['PC_SMILES'][i]
                bef_mega.loc[i, 'CompoundNames'] = bef_mega['PC_Name'][i]
                #bef_mega['most_specific_class'][i] = np.nan
                #bef_mega['level _5'][i] = np.nan
                bef_mega['subclass'][i] = np.nan
                bef_mega['class'][i] = np.nan
                bef_mega['superclass'][i] = np.nan
                #bef_mega['all_classifications'][i] = np.nan
                bef_mega['Classification_Source'][i] = np.nan
                bef_mega['MCSS_SMILES'][i] = np.nan
                bef_mega['KG_MCSS_SMILES'][i] = np.nan
                bef_mega['Formula'][i] = np.nan
            
            elif 'HMDB' in bef_mega['Annotation_Source'][i]:
                bef_mega.loc[i, 'SMILES_final'] = bef_mega['HMDBSMILES'][i]
                bef_mega.loc[i, 'CompoundNames'] = bef_mega['HMDBcompound_name'][i]
                #bef_mega['most_specific_class'][i] = np.nan
                #bef_mega['level _5'][i] = np.nan
                bef_mega['subclass'][i] = np.nan
                bef_mega['class'][i] = np.nan
                bef_mega['superclass'][i] = np.nan
                #bef_mega['all_classifications'][i] = np.nan
                bef_mega['Classification_Source'][i] = np.nan
                bef_mega['MCSS_SMILES'][i] = np.nan
                bef_mega['PC_MCSS_SMILES'][i] = np.nan
                bef_mega['KG_MCSS_SMILES'][i] = np.nan
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
            
    bef_megaA.rename(columns = {'SMILES_final':'SMILES'}, inplace = True)
    
    
    Standards = ['Experimental']
    SpectralDB = ['GNPS', 'HMDB', 'MassBank']
    CompoundDB = ['SuspectList', 'SIRIUS', 'KEGG', 'PubChem']
    Formula = ['SIRIUS_Formula']

    
    #bef_megaA['MSI_Level'] = np.nan
    for i, rows in bef_megaA.iterrows():
        
        
        if not isNaN(bef_megaA['Annotation_Source'][i]):
            
            if data_type == "standards":
                bef_megaA.loc[i, 'Annotation_Source'] = bef_megaA['Annotation_Source'][i] + ', Experimental'

                if any(x in bef_megaA['Annotation_Source'][i] for x in SpectralDB):
                    bef_megaA.loc[i, 'MSI_Level'] = 'Level_1'
                    
                elif any(x in bef_megaA['Annotation_Source'][i] for x in CompoundDB) and not any(x in bef_megaA['Annotation_Source'][i] for x in Formula):
                    bef_megaA.loc[i, 'MSI_Level'] = 'Level_2/Level_3'
                    
                elif any(x in bef_megaA['Annotation_Source'][i] for x in Formula):
                    bef_megaA.loc[i, 'MSI_Level'] = 'Level_4'
                    
            else:

                if any(x in bef_megaA['Annotation_Source'][i] for x in SpectralDB):
                    bef_megaA.loc[i, 'MSI_Level'] = 'Level_2'
                    
                elif any(x in bef_megaA['Annotation_Source'][i] for x in CompoundDB) and not any(x in bef_megaA['Annotation_Source'][i] for x in Formula):
                    bef_megaA.loc[i, 'MSI_Level'] = 'Level_3'
                    
                elif any(x in bef_megaA['Annotation_Source'][i] for x in Formula):
                    bef_megaA.loc[i, 'MSI_Level'] = 'Level_4'
                
        else:
            bef_megaA.loc[i, 'MSI_Level'] = 'Level_5'
            
                
    
            
    bef_megaA.to_csv(input_dir + "MetabolomicsResults/final_curation_without_classes.csv")
    return(bef_megaA)


# In[65]:


#print(combine_CuratedR.__doc__)


# In[66]:


#combine_CuratedR(input_dir = "/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/", 
                 #combinedSDBs = "/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/MetabolomicsResults/curatedSDB.csv", 
                 #combinedSMs = "/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/MetabolomicsResults/combinedSM.csv",
                #data_type = "standards")


# In[67]:


#combine_CuratedR(input_dir = "/Users/mahnoorzulfiqar/Downloads/MAW-main/", 
                 #combinedSDBs = "/Users/mahnoorzulfiqar/Downloads/MAW-main/MetabolomicsResults/curatedSDB.csv", 
                 #combinedSMs ="/Users/mahnoorzulfiqar/Downloads/MAW-main/MetabolomicsResults/combinedSM.csv",
                #data_type = "standards")


# In[68]:


def checkSMILES_validity(input_dir, resultcsv):
    def isNaN(string):
        return string != string
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
        if not isNaN(results['SMILES'][i]):
            m = Chem.MolFromSmiles(results['SMILES'][i] ,sanitize=False)
            if m is None:
                results['SMILES_final'][i] = 'invalid_SMILES'
            else:
                try:
                    Chem.SanitizeMol(m)
                except:
                    results['SMILES_final'][i] = 'invalid_chemistry'
    results.to_csv(input_dir + "MetabolomicsResults/final_curation_with_validSMILES.csv")
    return(results)


# In[69]:


#print(checkSMILES_validity.__doc__)


# In[70]:


#checkSMILES_validity(input_dir = '/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/',
  #                   resultcsv = '/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/MetabolomicsResults/final_curation_without_classes.csv')


# In[71]:


#checkSMILES_validity(input_dir = '/Users/mahnoorzulfiqar/Downloads/MAW-main/', 
      #               resultcsv = '/Users/mahnoorzulfiqar/Downloads/MAW-main/MetabolomicsResults/final_curation_without_classes.csv')


# In[ ]:





# In[72]:


def classification(input_dir, resultcsv):
    def isNaN(string):
        return string != string
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
        if not isNaN(frame['SMILES'][i]) and isNaN(frame['Classification_Source'][i]):
            try:
                InChI = Chem.MolToInchi(Chem.MolFromSmiles(frame["SMILES"][i]))
                InChIKey = Chem.inchi.InchiToInchiKey(InChI)
                inchis.append({
                    'index': i,
                    'smiles':frame["SMILES"][i],
                    'inchi': InChI,
                    'inchikey': InChIKey
                })
            except:
                pass
    inchis = pd.DataFrame(inchis)
    if len(inchis):
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
                if df_merged["smiles_x"][p] is frame["SMILES"][q]:
                    frame.loc[q, 'subclass'] = df_merged["subclass.name"][p]
                    frame.loc[q, 'class'] = df_merged["class.name"][p]
                    frame.loc[q, 'superclass'] = df_merged["superclass.name"][p]
                    frame.loc[q, 'Classification_Source'] = "ClassyFire"



        frame.to_csv(input_dir + "MetabolomicsResults/final_curationList.csv")
        return(frame)


# In[73]:


#print(classification.__doc__)


# In[ ]:





# In[ ]:





# In[74]:


#classification(input_dir = '/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/',
              # resultcsv = '/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/MetabolomicsResults/final_curation_with_validSMILES.csv')


# In[75]:


#classification(input_dir = '/Users/mahnoorzulfiqar/Downloads/MAW-main/',
               #resultcsv = '/Users/mahnoorzulfiqar/Downloads/MAW-main/MetabolomicsResults/final_curation_with_validSMILES.csv')


# # Comparison with a list of SMILES from any Source

# In[76]:


#cd = pd.read_csv('/Users/mahnoorzulfiqar/OneDriveUNI/MZML/CD/CD_Results.csv')


# In[77]:


#def create_SMILES_list(input_dir, compcsv)
    #CDCSV = list(cd[-isNaN(cd['SMILES'])]['SMILES'])
    #f = open("/Users/mahnoorzulfiqar/OneDriveUNI/MZML/CD/CDCSV.txt", "w")
    #for item in CDCSV:
        #f.write(item + "\n")
    #f.close()


# In[78]:


def SMILESscreening(input_dir, resultcsv, complist, listname):
    def isNaN(string):
        return string != string
    """SMILESscreening takes a list of SMILES

    Parameters:
    input_dir (str): This is the input directory where all the .mzML 
    files and their respective result directories are stored.
    
    resultcsv: df from combine_CuratedR or checkSMILES_validity or classification
    complist: list of /n separated txt file conyaining smiles on each line
    listname: name of the list of compounds
    
    Returns:
    dataframe: comparison with another list of compounds
    csv: "MetabolomicsResults/final_curation_with_validSMILES.csv"
    
    Usage:
    checkSMILES_validity(input_dir = "usr/project/", results)

    """
    
    results = pd.read_csv(resultcsv)
    with open(complist, "r") as text_file:
        cd = text_file.read().split('\n')
    
    for i, row in results.iterrows():
        if not isNaN(results['SMILES'][i]):
            if 'invalid_SMILES' not in results['SMILES'][i] and 'invalid_chemistry' not in results['SMILES'][i]:
                for j in cd:
                    if not isNaN(j):
                        CGms = [Chem.MolFromSmiles(results['SMILES'][i]), Chem.MolFromSmiles(j)]
                        CGfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=1024) for x in CGms]
                        CGtn = DataStructs.FingerprintSimilarity(CGfps[0],CGfps[1])
                        if CGtn == 1 and listname not in results['Annotation_Source'][i]:
                            results['Annotation_Source'][i] = results['Annotation_Source'][i] + ', ' + listname
    

    frame.to_csv(input_dir + "MetabolomicsResults/final_curationListVS"+listname+".csv")
    return(frame)


# In[79]:


#print(SMILESscreening.__doc__)


# In[ ]:





# In[80]:


#SMILESscreening(input_dir = '/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/',
                #resultcsv = '/Users/mahnoorzulfiqar/OneDriveUNI/CheckDocker/MetabolomicsResults/final_curation_with_validSMILES.csv',
                #complist = '/Users/mahnoorzulfiqar/OneDriveUNI/MZML/CD/CDCSV.txt',
                #listname = 'CompoundDiscoverer')


# In[81]:


#SMILESscreening(input_dir = '/Users/mahnoorzulfiqar/Downloads/MAW-main/', 
                #resultcsv = '/Users/mahnoorzulfiqar/Downloads/MAW-main/MetabolomicsResults/final_curationList.csv', 
                #complist = '/Users/mahnoorzulfiqar/OneDriveUNI/MZML/CD/CDCSV.txt', 
                #listname = 'CompoundDiscoverer')


# ## NP_Classifier classification

# In[82]:



def Np_pathways(input_dir, resultcsv):
    df = pd.read_csv(resultcsv)
    npresults = []
    for i, row in df.iterrows():
        if not isNaN(df['SMILES'][i]):
            try:
                cvv = Chem.MolFromSmiles(df['SMILES'][i])
                cvv = Chem.MolToSmiles(cvv, isomericSmiles = False)
                c = urllib.parse.quote_plus(cvv, safe=' ')
            
                url = 'https://npclassifier.ucsd.edu/classify?smiles='+c
                names = str(df['id_X'][i])
                outx = str("myFile"+names+".txt")
                file = wget.download(url, out = outx)
                a_dataframe = pd.read_csv(file, delimiter = "]")
                xox = list(a_dataframe.columns.values)
                splitting0 = xox[0].split(':')
                xoc = re.sub('\ |\[|\]|\"', ' ', splitting0[1]).strip()
                splitting1 = xox[1].split(':')
                xos = re.sub('\ |\[|\]|\"', ' ', splitting1[1]).strip()
                #except:
                    #splitting1 = xox[1].split(':')
                    #xos = re.sub('\ |\[|\]|\"', '', splitting1[0])
                splitting2 = xox[2].split(':')
                xop = re.sub('\ |\[|\]|\"', ' ', splitting2[1]).strip()
                #df.loc[i, 'npclass'] = xoc
                #df.loc[i, 'npsuper_class'] = xos
                if not isNaN(df['class'][i]) and df['class'][i] in xoc:
                    df.loc[i, 'np_pathway'] = xop
                os.remove(outx)
                time.sleep(0.5)

                npresults.append({
                    'index':i,
                    #'id': df['file_id'][i],
                    'mz': df['premz'][i],
                    'rt': df['rtmed'][i],
                    'SMILES': df['SMILES'][i],
                    'class': xoc,
                    'subclass': xos,
                    'pathway': xop
                })
            except:
                pass
    np_results = pd.DataFrame(npresults)
    np_results.to_csv(input_dir + "/MetabolomicsResults/NPClassifier_Results.csv")
    df.to_csv(input_dir + "/MetabolomicsResults/final_results_with_Pathways.csv")


# In[83]:


#Np_pathways(input_dir = "/Users/mahnoorzulfiqar/OneDriveUNI/MZML", 
            #resultcsv = "/Users/mahnoorzulfiqar/OneDriveUNI/MZML/MetabolomicsResults/Final_Candidate_List.csv")


# In[84]:


#df = pd.read_csv("/Users/mahnoorzulfiqar/OneDriveUNI/MZML/MetabolomicsResults/Final_Candidate_List.csv")


# In[ ]:





# ## Chemical Similarity MN

# In[86]:


def chemMN(input_dir, resultcsv):
    #read csv
    df = pd.read_csv(resultcsv)
    
    # define empty variable
    dbn= []

    # check the result csv
    for i, row in df.iterrows():
        # to compare each element with each opther element
        for j, row in df.iterrows():

            # if its not same id
            if df['SMILES'][i] != df['SMILES'][j]:

                if not isNaN(df['SMILES'][i]):
                    if not isNaN(df['SMILES'][j]):

                        try:
                            ms = [Chem.MolFromSmiles(df['SMILES'][i]), Chem.MolFromSmiles(df['SMILES'][j])]
                            fps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in ms]
                            tn = DataStructs.FingerprintSimilarity(fps[0],fps[1])
                            dbn.append({
                                'Name_i':df['id_X'][i],
                                'Name_j':df['id_X'][j],
                                'i': df['SMILES'][i],
                                'j': df['SMILES'][j],
                                'Tanimoto': tn
                            })
                        except Exception as e:
                            print(i)
                            print(j)
                            print(e)
    # save chemical similarities                    
    db_edgenode = pd.DataFrame(dbn)

    dfe = []
    heavy_atoms = ['C', 'N', 'P', 'O', 'S']
    for i, row in db_edgenode.iterrows():        
        if 1.0 > db_edgenode['Tanimoto'][i] >= 0.70:
            # list of mol used to calaculate the MCSS
            n = [Chem.MolFromSmiles(db_edgenode['i'][i]),Chem.MolFromSmiles(db_edgenode['j'][i])]
            res = rdFMCS.FindMCS(n)
            sm_res = Chem.MolToSmiles(Chem.MolFromSmarts(res.smartsString))
            # Check if the MCSS has one of the heavy atoms and whether they are
            # more than 3
            elem = [ele for ele in heavy_atoms if(ele in sm_res)]
            if elem and len(sm_res)>=3:
                MCSS_SMILES = Chem.MolToSmiles(Chem.MolFromSmarts(res.smartsString))

            dfe.append({
                'Start':db_edgenode['Name_i'][i],
                'End':db_edgenode['Name_j'][i],
                'Tanimoto':db_edgenode['Tanimoto'][i],
                'Start_SMILES':db_edgenode['i'][i],
                'End_SMILES':db_edgenode['j'][i],
                #'Start_Source':db_edgenode['Source_i'][i],
                #'End_Source':db_edgenode['Source_j'][i],
                'MCSS': MCSS_SMILES
            })

    df_edge = pd.DataFrame(dfe)
    df_edge['Start'] = df_edge['Start'].astype(str)
    df_edge['End'] = df_edge['End'].astype(str)
    df_edge['sorted_row'] = [sorted([a,b]) for a,b in zip(df_edge.Start,df_edge.End)]
    df_edge['sorted_row'] = df_edge['sorted_row'].astype(str)
    df_edge.drop_duplicates(subset=['sorted_row'], inplace=True)

    nodes= []
    for i, row in df.iterrows():
        n = df['id_X'][i]
        nodes.append({
            'nodes':n
        })

    node= pd.DataFrame(nodes)
    
    
    df_edge.to_csv(input_dir + "/MetabolomicsResults/ChemMNedges.tsv", sep='\t')
    node.to_csv(input_dir + "/MetabolomicsResults/ChemMNnodes.csv", index = False)

    newdf = df_edge
    newdf['StartAtt']=np.nan
    newdf['EndAtt']=np.nan
    for i, row in newdf.iterrows():
        for j, row in df.iterrows():
            if newdf['Start'][i]==df['id_X'][j]:
                newdf.loc[i, 'StartAtt'] = df['class'][j]
            if newdf['End'][i]==df['id_X'][j]:
                newdf.loc[i, 'EndAtt'] = df['class'][j]
    newdf.to_csv(input_dir + "/MetabolomicsResults/ChemMNcys.tsv", sep='\t')
    
    return(newdf)

    


# In[101]:


#newdf = chemMN(input_dir = "/Users/mahnoorzulfiqar/OneDriveUNI/MZML", 
            #resultcsv = "/Users/mahnoorzulfiqar/OneDriveUNI/MZML/MetabolomicsResults/Final_Candidate_List.csv")


# In[114]:


#input_dir = "/Users/mahnoorzulfiqar/OneDriveUNI/MZML"


# In[87]:


#save_tuples = []
#my_list = newdf["sorted_row"]
#for i in my_list:
    #save_tuples.append(str(i).replace('[','(').replace(']',')'))
#save_tuples
#list(save_tuples)
#G=nx.from_edgelist(save_tuples)


# ## Molecular Networking

# ## MN with GNPS

# In[88]:


def gnpsMNvsgnpsMAW(input_dir):
    def isNaN(string):
        return string != string
    """gnpsMNvsgnpsMAW checks with tanimoto similarity score, whether
    results from MAW GNPS and GNPS MN Masst results give same candidate

    Parameters:
    input_dir = input directory where you have stored the cytoscape file 
    from GNPS MN results and have exported edge and node tables from cytoscape
    These two csv egde and node files must have "edge" and "node" in their name
    
    Returns:
    GNPS results with cluster index named
    GNPS MN results with a confirmation column if MAW detected same candidate, 
    file named: 
    
    Usage: 
    gnpsMNvsgnpsMAW(input_dir)
    
    """
    # extract files with edges from MN results
    GMNfile_edge = [f for f in os.listdir(input_dir) if "edge" in f]
    # extract files with nodes from MN results
    GMNfile_node = [f for f in os.listdir(input_dir) if "node" in f]
    # read the files
    GMNdf_node = pd.read_csv(GMNfile_node[0])
    GMNdf_edge = pd.read_csv(GMNfile_edge[0])
    
    # extract only important columns from both csv files
    GMNdf_node = GMNdf_node[['precursor mass', 'RTMean', 'UniqueFileSources', 
                   'charge', 'cluster index', 'componentindex', 
                   'Compound_Name', 'Smiles', 'SpectrumID']]
    GMNdf_edge = GMNdf_edge[['cosine_score', 'EdgeAnnotation', 'node1', 'node2',
                     'mass_difference']]
    
    # rename node1 to cluster index to merge nodes and edges results from MN
    GMNdf_edge = GMNdf_edge.rename(columns={'node1': 'cluster index'})
    GMNdf = pd.merge(GMNdf_node, GMNdf_edge, on = "cluster index")
    
    # Read results obtained from scoring_spec, named input_dir/MetabolomicsResults/scoredSpecDB.csv
    SDB = pd.read_csv(input_dir + "/MetabolomicsResults/scoredSpecDB.csv")
    # only keep GNPS resulst and remove other columns
    only_GNPS = SDB[SDB['annotation'].str.contains('GNPS')]
    only_GNPS = only_GNPS[['id_X', 'premz_x', 'rtmean_x', 'GNPSmax_similarity', 
                       'GNPSSMILES', 'GNPSspectrumID', 'GNPScompound_name', 
                       'GNPSmirrorSpec']]
    
    # from GNPS MAW results and GNPS MN results, calculate how many MAW results are same as MN:
    for i, row in only_GNPS.iterrows():
        for j, row in GMNdf.iterrows():
            if not isNaN(only_GNPS["GNPSSMILES"][i]) and not isNaN(GMNdf["Smiles"][j]):
                SKms = [Chem.MolFromSmiles(only_GNPS['GNPSSMILES'][i]), Chem.MolFromSmiles(GMNdf['Smiles'][j])]
                SKfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in SKms]
                SKtn = DataStructs.FingerprintSimilarity(SKfps[0],SKfps[1])
                if SKtn == 1.0:
                    GMNdf.loc[j, "gnps_maw"] = "confirmed"
                    only_GNPS.loc[i, "index_MN_nodes"] = j
                elif SKtn < 1.0 and SKtn < 0.75:
                    GMNdf.loc[j, "gnps_maw"] = "similar"
                    only_GNPS.loc[i, "index_MN_nodes"] = j
    only_GNPS.to_csv(input_dir + "/MetabolomicsResults/only_GNPS.csv")
    GMNdf.to_csv(input_dir + "/MetabolomicsResults/GMNdf.csv")


# In[89]:


#gnpsMNvsgnpsMAW(input_dir = "/Users/mahnoorzulfiqar/Downloads/MAW-main")


# In[90]:


# get the compounds in the cluster extracted out also in the same fucntion and store in a csv file


# In[91]:


# calculate the MN vs MCSS results to see any


# In[ ]:





# ### MN vs MCSS

# In[5]:


#final_list = pd.read_csv("/Users/mahnoorzulfiqar/OneDriveUNI/MZML/MetabolomicsResults/Final_Candidate_List.csv")


# In[13]:


#final_list = final_list.rename(columns = {'SMILES_final':'SMILES'})


# In[14]:


#final_list.to_csv("/Users/mahnoorzulfiqar/OneDriveUNI/MZML/MetabolomicsResults/Final_Candidate_List.csv")


# In[ ]:




