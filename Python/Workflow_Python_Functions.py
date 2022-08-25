#!/usr/bin/env python
# coding: utf-8

# In[1]:


from platform import python_version

#print(python_version())


# In[2]:


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


# In[3]:


from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem import PandasTools


# In[4]:


input_dir = "/Users/mahnoorzulfiqar/OneDriveUNI/MAW-data/StandardSMarinoi_Data"
input_dir


# # Suspect List for SIRIUS

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


# In[6]:


#print(slist_metfrag.__doc__)
#print(slist_sirius.__doc__)


# ## Spectral DB dereplication Results PostProcessing

# In[7]:


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
    
    # Define scoring for all DBs
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
                    # GNPS Results
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
                            files = (glob.glob(sub_dir+'/*.csv'))
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
                    
                    # HMDB Results
                    if Source == "hmdb" or Source == "all":                 

                        if not os.path.exists(input_dir+"/structures.sdf"):
                            #download SDF structures
                            os.system("wget -P " + input_dir + " https://hmdb.ca/system/downloads/current/structures.zip")
                            os.system("unzip "+ input_dir + "/structures.zip" + " -d " + input_dir)
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
                            files = (glob.glob(sub_dir+'/*.csv'))
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
                    
                    # MASSBANK Results
                    
                    # enter the directory with /spectral_dereplication/ results
                    if Source == "mbank" or Source == "all":
                        #open another csv path holding empty list, which will be filled 
                        #with post processed csv results
                        MassBankcsvfiles2 = []
                        #print(entry)
                        # enter the directory with /spectral_dereplication/ results
                        sub_dir = input_dir + "/" + entry + '/spectral_dereplication/MassBank/'
                        if os.path.exists(sub_dir):
                            files = (glob.glob(sub_dir+'/*.csv'))
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


# In[8]:


#spec_postproc(input_dir, Source = "mbank")


# # SIRIUS Post Processing

# In[7]:


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
                        if os.path.exists(sub_sub_dirSL_structure_can) and len(pd.read_csv(sub_sub_dirSL_structure_can, sep = "\t"))>0 and os.path.exists(sub_sub_dirALL_structure_can) and len(pd.read_csv(sub_sub_dirALL_structure_can, sep = "\t"))>0:

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
                                        if len(ALL_Canopus)>0:
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
                                        if len(SL_Canopus)>0:
                                            if SL_formula_csv["molecularFormula"][formula_sl] == SL_Canopus["molecularFormula"][0]:
                                                SL_structure_csv.loc[structure_sl, 'superclass'] = SL_Canopus['superclass'][0]
                                                SL_structure_csv.loc[structure_sl, 'class'] = SL_Canopus['class'][0]
                                                SL_structure_csv.loc[structure_sl, 'subclass'] = SL_Canopus['subclass'][0] 


                            # after formula and structure have been merged, merge SL and ALL results
                            all_sl_db = pd.concat([ALL_structure_csv, SL_structure_csv], ignore_index=True)
                            for str_sirius, row in all_sl_db.iterrows():
                                if not str_can_score(all_sl_db, str_sirius):
                                    all_sl_db = all_sl_db.drop(str_sirius, inplace=False)

                            result_sirius_name = (sub_dir+"results_for_"+json_dirALL[0].split("_")[-1]+"_"+"structure.csv")
                            msp.loc[mz, "sirius_result_dir"] = result_sirius_name.replace(input_dir, ".")
                            all_sl_db.to_csv(result_sirius_name)

                        # if only ALL structure file exists and they have more than 0 rows
                        elif not (os.path.exists(sub_sub_dirSL_structure_can) and len(pd.read_csv(sub_sub_dirSL_structure_can, sep = "\t"))>0) and os.path.exists(sub_sub_dirALL_structure_can) and len(pd.read_csv(sub_sub_dirALL_structure_can, sep = "\t"))>0:
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
                                        if len(ALL_Canopus)>0:
                                            if ALL_formula_csv["molecularFormula"][formula] == ALL_Canopus["molecularFormula"][0]:
                                                ALL_structure_csv.loc[structure, 'superclass'] = ALL_Canopus['superclass'][0]
                                                ALL_structure_csv.loc[structure, 'class'] = ALL_Canopus['class'][0]
                                                ALL_structure_csv.loc[structure, 'subclass'] = ALL_Canopus['subclass'][0]
                            # scoring function
                            for str_siriusA, row in ALL_structure_csv.iterrows():
                                if not str_can_score(ALL_structure_csv, str_siriusA):
                                    ALL_structure_csv = ALL_structure_csv.drop(str_siriusA, inplace=False)

                            result_sirius_name = (sub_dir+"results_for_"+json_dirALL[0].split("_")[-1]+"_"+"structure.csv")
                            msp.loc[mz, "sirius_result_dir"] = result_sirius_name.replace(input_dir, ".")
                            ALL_structure_csv.to_csv(result_sirius_name)

                        elif not(os.path.exists(sub_sub_dirSL_structure_can) and len(pd.read_csv(sub_sub_dirSL_structure_can, sep = "\t"))>0 and os.path.exists(sub_sub_dirALL_structure_can) and len(pd.read_csv(sub_sub_dirALL_structure_can, sep = "\t"))>0):
                            if os.path.exists(sub_sub_dirALL_formula_can) and len(pd.read_csv(sub_sub_dirALL_formula_can, sep = "\t"))>0:
                                ALL_formula_csv = pd.read_csv(sub_sub_dirALL_formula_can, sep = "\t")
                                ALL_Canopus = pd.read_csv(ALL_Canopus_csv, sep = "\t")
                                if len(ALL_Canopus)>0:
                                    for formula, rows in ALL_formula_csv.iterrows():
                                        ALL_formula_csv.loc[formula, 'superclass'] = ALL_Canopus['superclass'][0]
                                        ALL_formula_csv.loc[formula, 'class'] = ALL_Canopus['class'][0]
                                        ALL_formula_csv.loc[formula, 'subclass'] = ALL_Canopus['subclass'][0]

                                for for_siriusA, row in ALL_formula_csv.iterrows():
                                    if not ALL_formula_csv['explainedIntensity'][for_siriusA] >= exp_int:
                                        ALL_formula_csv = ALL_formula_csv.drop(for_siriusA, inplace=False)

                                result_sirius_name = (sub_dir+"results_for_"+json_dirALL[0].split("_")[-1]+"_"+"formula.csv")
                                msp.loc[mz, "sirius_result_dir"] = result_sirius_name.replace(input_dir, ".")

                                ALL_formula_csv.to_csv(result_sirius_name)

                            else:
                                pass
                        else:
                            pass

                    elif len(files_for_mz) == 1:
                        # extract the formula and structure files
                        json_dirALL = next(os.walk(files_for_mz[0]))[1]
                        sub_sub_dirALL_structure_can = files_for_mz[0] + "/" + json_dirALL[0]  + "/structure_candidates.tsv"                   
                        sub_sub_dirALL_formula_can = files_for_mz[0] + "/" + json_dirALL[0]  + "/formula_candidates.tsv" 
                        ALL_Canopus_csv = files_for_mz[0] + "/canopus_summary.tsv"

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
                                        if len(ALL_Canopus)>0:
                                            if ALL_formula_csv["molecularFormula"][formula] == ALL_Canopus["molecularFormula"][0]:
                                                ALL_structure_csv.loc[structure, 'superclass'] = ALL_Canopus['superclass'][0]
                                                ALL_structure_csv.loc[structure, 'class'] = ALL_Canopus['class'][0]
                                                ALL_structure_csv.loc[structure, 'subclass'] = ALL_Canopus['subclass'][0]

                            for str_siriusA, row in ALL_structure_csv.iterrows():
                                if not str_can_score(ALL_structure_csv, str_siriusA):
                                    ALL_structure_csv = ALL_structure_csv.drop(str_siriusA, inplace=False)

                            result_sirius_name = (sub_dir+"results_for_"+json_dirALL[0].split("_")[-1]+"_"+"structure.csv")
                            msp.loc[mz, "sirius_result_dir"] = result_sirius_name.replace(input_dir, ".")

                            ALL_structure_csv.to_csv(result_sirius_name)


                        elif not (os.path.exists(sub_sub_dirALL_structure_can) and len(pd.read_csv(sub_sub_dirALL_structure_can, sep = "\t"))>0):
                            if os.path.exists(sub_sub_dirALL_formula_can) and len(pd.read_csv(sub_sub_dirALL_formula_can, sep = "\t"))>0:
                                ALL_formula_csv = pd.read_csv(sub_sub_dirALL_formula_can, sep = "\t")
                                ALL_Canopus = pd.read_csv(ALL_Canopus_csv, sep = "\t")
                                if len(ALL_Canopus)>0:
                                    for formula, rows in ALL_formula_csv.iterrows():
                                        ALL_formula_csv.loc[formula, 'superclass'] = ALL_Canopus['superclass'][0]
                                        ALL_formula_csv.loc[formula, 'class'] = ALL_Canopus['class'][0]
                                        ALL_formula_csv.loc[formula, 'subclass'] = ALL_Canopus['subclass'][0]

                                for for_siriusA, row in ALL_formula_csv.iterrows():
                                    if not ALL_formula_csv['explainedIntensity'][for_siriusA] >= exp_int:
                                        ALL_formula_csv = ALL_formula_csv.drop(for_siriusA, inplace=False)

                                result_sirius_name = (sub_dir+"results_for_"+json_dirALL[0].split("_")[-1]+"_"+"formula.csv")
                                msp.loc[mz, "sirius_result_dir"] = result_sirius_name.replace(input_dir, ".")

                                ALL_formula_csv.to_csv(sub_dir+"results_for_"+json_dirALL[0].split("_")[-1]+"_"+"formula.csv")

                            else:
                                print("no file for formula")
                        else:
                            print("no file for structure or formula")
                msp.to_csv(msp_csv)
                


# In[ ]:





# In[10]:


#sirius_postproc(input_dir, exp_int = 0.90, csi_score = -150)


# # MCSS for SpecDBs

# In[11]:


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


# In[12]:


#MCSS_for_SpecDB(input_dir, Source = "all")


# # MCSS for SIRIUS

# In[13]:


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
                for sir_file in sirius_files:
                    r = [s for s in os.listdir(sir_file) if "structure" in s]
                    for filenames in r:
                        sirius_f = sir_file+filenames
                        for mz, row in sirius_msp.iterrows():
                            if str(sirius_msp["id_X"][mz].split("_")[1]) in sirius_f:
                                s_f = pd.read_csv(str(sirius_f))
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

# In[1]:


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


def chemMN_CandidateSelection(df, tn_sim = 0.85):
    
    """chemMN_CandidateSelection function is used to generate a Cytoscape readable tsv file.
    This file contains start(starting SMILES) and end(target SMILES) nodes and the tanimoto 
    similarity scores between the nodes. User can visualize the structural similarity 
    between the given SMILES. It provides an "ALL against ALL" network.
    
    Parameters: 
    df: dataframe that contains "SMILES", "ranks", "Source". This function is specifically for 
    candidate selection and so these columns are necessary.
    
    
    Returns:
    dataframe: it returns a df with follwoing columns to be loaded into Cytoscape.
    1. Start, starting node/SMILES
    2. End, ending node/SMILES
    3. Tanimoto, Tanimoto between Start and End node
    4. Start_SMILES
    5. End_SMILES
    6. Start_Source
    7. End_Source
    8. MCSS, Maximum Common Substructure between start and end node/SMILES
    9. sorted_row, contains ids of the start and end nodes as a list
    
    
    Usage:
    chemMN_CandidateSelection(df)

    """
    
    # define an empty variable 
    #one_df = []
    # define empty variable to save the edges
    dbn= []
    # for each entry in the df
    for i, row in df.iterrows():
        # to compare each element with each other element of the df
        for j, row in df.iterrows():
            try:
                # calcultae tanimoto
                ms = [Chem.MolFromSmiles(df['SMILES'][i]), Chem.MolFromSmiles(df['SMILES'][j])]
                fps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in ms]
                tn = DataStructs.FingerprintSimilarity(fps[0],fps[1])
                
                # save all entries to a matrix
                dbn.append({
                    'Name_i':df['ranks'][i],
                    'Name_j':df['ranks'][j],
                    'i': df['SMILES'][i],
                    'j': df['SMILES'][j],
                    'Source_i':df['Source'][i],
                    'Source_j':df['Source'][j],
                    'Tanimoto': tn
                })

            except Exception as e:
                #print(e.string)
                pass

    # save chemical similarities                    
    db_edgenode = pd.DataFrame(dbn)

    # another empty variable to store the results for final tsv file
    dfe = []
    
    # heavy atoms for MCSS Calculation
    heavy_atoms = ['C', 'N', 'P', 'O', 'S']
    
    # for the previous dataframe
    for i, row in db_edgenode.iterrows():      
        # if the tanimoto > 0.85 for high similarity
        if db_edgenode['Tanimoto'][i] >= tn:
            
            # calculate MCSS
            n = [Chem.MolFromSmiles(db_edgenode['i'][i]),Chem.MolFromSmiles(db_edgenode['j'][i])]
            res = rdFMCS.FindMCS(n)
            sm_res = Chem.MolToSmiles(Chem.MolFromSmarts(res.smartsString))
            
            
            # Check if the MCSS has one of the heavy atoms and whether they are
            # more than 3
            elem = [ele for ele in heavy_atoms if(ele in sm_res)]
            if elem and len(sm_res)>=3:
                MCSS_SMILES = Chem.MolToSmiles(Chem.MolFromSmarts(res.smartsString))
            
            # save everything into a dataframe
            dfe.append({
                'Start':db_edgenode['Name_i'][i],
                'End':db_edgenode['Name_j'][i],
                'Tanimoto':db_edgenode['Tanimoto'][i],
                'Start_SMILES':db_edgenode['i'][i],
                'End_SMILES':db_edgenode['j'][i],
                'Start_Source':db_edgenode['Source_i'][i],
                'End_Source':db_edgenode['Source_j'][i],
                'MCSS': MCSS_SMILES
            })
    df_edge = pd.DataFrame(dfe)
    # generate a column called sorted_row which contains ids of the start and end nodes as a list
    df_edge['Start'] = df_edge['Start'].astype(str)
    df_edge['End'] = df_edge['End'].astype(str)
    df_edge['sorted_row'] = [sorted([a,b]) for a,b in zip(df_edge.Start,df_edge.End)]
    df_edge['sorted_row'] = df_edge['sorted_row'].astype(str)
    df_edge.drop_duplicates(subset=['sorted_row'], inplace=True)

    #nodes= []
    #for i, row in df.iterrows():
        #n = df['ranks'][i]
        #nodes.append({
            #'nodes':n
        #})

    #node= pd.DataFrame(nodes)

    return(df_edge)


# In[ ]:


def one_candidate_selection(df, Source = "SGHM", tn_ident = 0.99):
    
    """one_candidate_selection function is used to generate a dataframe that tells, 
    for each candidate SMILES, what was the source or how many sources had the same 
    candidate. The idea is to merge all candidate SMILES into one list, preserving
    the rank and source, and then checking whether these SMILES come from SIRIUS or
    any spectral DB. If a SMILE is repeated in more sources, its confidence score 
    increases and is considered the most likely candidate structure. This function 
    is not stand-alone and is part of the function CandidateSelection_SimilarityandIdentity

    
    Parameters: 
    df: dataframe that contains "SMILES", "ranks", "Source". This function is 
    specifically for candidate selection and so these columns are necessary.
    Source: this depends on how many sources were used. Possiblilities are: 
    1. SGHM (all)
    2. SGM (SIRIUS, GNPS, MassBank)
    3. SHM (SIRIUS, HMDB, MassBank)
    4. SGH (SIRIUS, GNPS, HMDB)
    5. GHM (GNPS, HMDB, MassBank)
    6. SG (SIRIUS, GNPS)
    7. SH (SIRIUS, HMDB)
    8. SM (SIRIUS, MassBank)
    9. GM (GNPS, MassBank)
    10. GH (GNPS, HMDB)
    11. HM (HMDB, MassBank)
    
    Returns:
    dataframe: it returns a df with follwoing columns which can be used to 
    prioritize a database for the final candidate selection.
    1. Source, contains name of the source (SIRIUS, GNPS, HMDB or MassBank)
    2. ranks, contains first letter of the source and a rank number seperated 
    by _ e.g: G_1(GNPS, 1st rank)
    3. SMILES
    4. SIRIUS, the rank again but only when the corresponding row SMILES is 
    also part of SIRIUS results
    5. GNPS , same as SIRIUS but for GNPS
    5. MassBank
    6. HMDB
    
    
    Usage:
    chemMN_CandidateSelection(df, Source = "SGHM")

    """
    
    # define empty columns for each Source to only fill if the corresponding 
    # SMILES is also present in the source
    
    df["SIRIUS"] = np.nan
    df["GNPS"] = np.nan
    df["MassBank"] = np.nan
    df["HMDB"] = np.nan
    
    # for each SMILES in df
    for smiles, rows in df.iterrows():
        
        # If the source contains SIRIUS
        if Source == "SGHM" or Source == "SGM" or Source == "SGH" or Source == "SHM" or Source == "SG" or Source == "SM" or Source == "SH":
            # sirius_df comes from within the function CandidateSelection_SimilarityandIdentity
            for sirius_i, row in sirius_df.iterrows():
                # calculate tanimoto
                try:
                    ms = [Chem.MolFromSmiles(df["SMILES"][smiles]), Chem.MolFromSmiles(sirius_df["smiles"][sirius_i])]
                    fps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in ms]
                    tn = DataStructs.FingerprintSimilarity(fps[0],fps[1])
                    # since we are dealing with idenity here so tanimoto of 0.99 is appropriate
                    if tn >= tn_ident:
                        
                        # if SIRIUS is blank, add the SIRIUS id
                        if isNaN(df["SIRIUS"][smiles]):

                            df.loc[smiles, "SIRIUS"] = sirius_df["rank_ids"][sirius_i]
                        # if not empty, add SIRIUS id, with a comma
                        else:
                            df.loc[smiles, "SIRIUS"] = str(df["SIRIUS"][smiles]) + ", "+ sirius_df["rank_ids"][sirius_i]
                            
                       
                except Exception as e:
                    print(e.string)
                    pass
        
        
        # If the Source contains GNPS
        if Source == "SGHM" or Source == "SGM" or Source == "SGH" or Source == "GHM" or Source == "SG" or Source == "GM" or Source == "GH":

            # gnps_df comes from within the function CandidateSelection_SimilarityandIdentity
            for gnps_i, row in gnps_df.iterrows():
                try:
                    # calculate tanimoto
                    ms = [Chem.MolFromSmiles(df["SMILES"][smiles]), Chem.MolFromSmiles(gnps_df["GNPSSMILES"][gnps_i])]
                    fps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in ms]
                    tn = DataStructs.FingerprintSimilarity(fps[0],fps[1])
                    
                    # since we are dealing with idenity here so tanimoto of 0.99 is appropriate
                    if tn >= tn_ident:
                        
                        # if GNPS is blank, add the GNPS id
                        if isNaN(df["GNPS"][smiles]):

                            df.loc[smiles, "GNPS"] = gnps_df["rank_ids"][gnps_i]
                        # if not empty, add GNPS id, with a comma
                        else:
                            df.loc[smiles, "GNPS"] = str(df["GNPS"][smiles]) + ", "+ gnps_df["rank_ids"][gnps_i]          
                        
                except Exception as e:
                    print(e.string)
                    pass

        # If the source contains MassBank
        if Source == "SGHM" or Source == "SGM" or Source == "SHM" or Source == "GHM" or Source == "SM" or Source == "GM" or Source == "HM":
            # mbank_df comes from within the function CandidateSelection_SimilarityandIdentity
            for mbank_i, row in mbank_df.iterrows():
                try:
                    # calculate tanimoto
                    ms = [Chem.MolFromSmiles(df["SMILES"][smiles]), Chem.MolFromSmiles(mbank_df["MBSMILES"][mbank_i])]
                    fps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in ms]
                    tn = DataStructs.FingerprintSimilarity(fps[0],fps[1])
                    
                    # since we are dealing with idenity here so tanimoto of 0.99 is appropriate
                    if tn >= tn_ident:
                        
                        # if MassBank is blank, add the MassBank id
                        if isNaN(df["MassBank"][smiles]):

                            df.loc[smiles, "MassBank"] = mbank_df["rank_ids"][mbank_i]
                        # if not empty, add MassBank id, with a comma
                        else:
                            df.loc[smiles, "MassBank"] = str(df["MassBank"][smiles]) + ", "+ mbank_df["rank_ids"][mbank_i]
                            
                       
                except Exception as e:
                    print(e.string)
                    pass


        # If the source contains HMDB
        if Source == "SGHM" or Source == "SGH" or Source == "SHM" or Source == "GHM" or Source == "SH" or Source == "GH" or Source == "HM":


            for hmdb_i, row in hmdb_df.iterrows():
                try:
                    # calculate tanimoto
                    ms = [Chem.MolFromSmiles(df["SMILES"][smiles]), Chem.MolFromSmiles(hmdb_df["HMDBSMILES"][hmdb_i])]
                    fps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in ms]
                    tn = DataStructs.FingerprintSimilarity(fps[0],fps[1])
                    
                    # since we are dealing with idenity here so tanimoto of 0.99 is appropriate
                    if tn >= tn_ident:
                        
                        # if HMDB is blank, add the HMDB id
                        if isNaN(df["HMDB"][smiles]):

                            df.loc[smiles, "HMDB"] = hmdb_df["rank_ids"][hmdb_i]
                        # if not empty, add HMDB id, with a comma
                        else:
                            df.loc[smiles, "HMDB"] = str(df["HMDB"][smiles]) + ", "+ hmdb_df["rank_ids"][hmdb_i]
                            
                       
                except Exception as e:
                    print(e.string)
                    pass
 
    return(df)


# In[ ]:


def CandidateSelection_SimilarityandIdentity(input_dir):

    # entry is all files and folders in input_dir
    for entry in os.listdir(input_dir):
         #if the entry is also a directory
        if os.path.isdir(os.path.join(input_dir, entry)):


            # reach spectra_dereplication folder
            sub_dir_spec = input_dir + "/" + entry + '/spectral_dereplication/'
            # reach SIRIUS results
            sub_dir_sir = input_dir + "/" + entry + '/insilico/SIRIUS/'

            #list of all csv files in the spectral dereplication foler
            spec_msp_csv = (glob.glob(input_dir + "/" + entry + '/spectral_dereplication' +'/*.csv'))
            # Sirius csv result file
            sir_msp_csv = input_dir + "/" + entry + "/insilico/MS1DATA.csv"

            # if both exist; which should be the case, even in case of 0 results
            if os.path.exists(sir_msp_csv) and os.path.exists(spec_msp_csv[0]):

                # read both csv files
                spec_msv = pd.read_csv(spec_msp_csv[0])
                sir_msv = pd.read_csv(sir_msp_csv)

                spec_msv = spec_msv[[
                    'id_X',
                    'premz',
                    'rtmin',
                    'rtmax',
                    'rtmed',
                    'rtmean',
                    'col_eng',
                    'pol',
                    'int',
                    'source_file',
                    'mbank_results_csv',
                    ]]
                sir_msv = sir_msv[['id_X',
                    'premz',
                    'rtmed',
                    'rtmean',
                    'int',
                    'col_eng',
                    'pol',
                    'ms2Peaks',
                    'ms1Peaks',
                    'sirius_result_dir']]

                merged_df = sir_msv.merge(spec_msv, 
                              how='inner', 
                              left_on=['premz', 'rtmed','rtmean','int','col_eng','pol'], 
                              right_on=['premz','rtmed','rtmean','int','col_eng','pol'])

                for mer, rows in merged_df.iterrows():
                    print(mer)

                    sirius_csv = merged_df["sirius_result_dir"][mer].replace("./", input_dir+"/")
                    mbank_csv = merged_df["mbank_results_csv"][mer].replace("./", input_dir+"/")
                    gnps_csv = merged_df["mbank_results_csv"][mer].replace("./", input_dir+"/").replace('mbank', 'gnps').replace('MassBank', 'GNPS')
                    hmdb_csv = merged_df["mbank_results_csv"][mer].replace("./", input_dir+"/").replace('mbank', 'hmdb').replace('MassBank', 'HMDB')


                    if os.path.exists(sirius_csv) and os.path.exists(gnps_csv) and os.path.exists(mbank_csv) and os.path.exists(hmdb_csv):


                        sirius_df = pd.read_csv(sirius_csv)
                        sirius_df = sirius_df.drop_duplicates('smiles')

                        mbank_df = pd.read_csv(mbank_csv)
                        mbank_df = mbank_df.drop_duplicates('MBSMILES')

                        gnps_df = pd.read_csv(gnps_csv)
                        gnps_df = gnps_df.drop_duplicates('GNPSSMILES')

                        hmdb_df = pd.read_csv(hmdb_csv)
                        hmdb_df = hmdb_df.drop_duplicates('HMDBSMILES')

                        #1 SGHM
                        if len(sirius_df) > 0 and len(gnps_df) > 0 and len(mbank_df) > 0 and len(hmdb_df) > 0:
                            mbank_df["rank_ids"] = ["M_" + str(s+1) for s in range(len(mbank_df))]

                            gnps_df["rank_ids"] = ["G_" + str(s+1) for s in range(len(gnps_df))]

                            hmdb_df["rank_ids"] = ["H_" + str(s+1) for s in range(len(hmdb_df))]

                            sirius_df["rank_ids"] = ["S_" + str(s) for s in sirius_df["rank"]]
                            sirius_df["Source"] = "SIRIUS"


                            source_l1 = [*(list(sirius_df["Source"])), 
                               *(list(gnps_df["Source"]))
                               ,*(list(mbank_df["Source"])),
                                        *(list(hmdb_df["Source"]))]

                            rank_l2 = [*(list(sirius_df["rank_ids"])), 
                               *(list(gnps_df["rank_ids"]))
                               ,*(list(mbank_df["rank_ids"])), 
                                      *(list(hmdb_df["rank_ids"]))]

                            smiles_l3 = [*(list(sirius_df["smiles"])), 
                               *(list(gnps_df["GNPSSMILES"]))
                               ,*(list(mbank_df["MBSMILES"])),
                                        *(list(hmdb_df["HMDBSMILES"]))]

                            sm = pd.DataFrame(list(zip(source_l1, rank_l2, smiles_l3)), columns = ["Source", "ranks", "SMILES"])   
                            
                            df_edge = chemMN_CandidateSelection(sm)

                            df_edge.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_ChemMNedges.tsv", sep='\t')
                            one_candidate = one_candidate_selection(sm, Source = "SGHM")
                            one_candidate.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_one_candidate_list.tsv", sep='\t')
                            
                            
                            df_sources = pd.DataFrame({"SIRIUS": one_candidate["SIRIUS"], 
                                                       "GNPS": one_candidate["GNPS"],
                                                       "HMDB": one_candidate["HMDB"],
                                                       "MassBank": one_candidate["MassBank"]})
                            
                            
                            # now check which rows have a value

                            index_SIRIUS = [x for x, row in df_sources.iterrows() if not isNaN(df_sources["SIRIUS"][x])]
                            index_GNPS = [x for x, row in df_sources.iterrows() if not isNaN(df_sources["GNPS"][x])]
                            index_MassBank = [x for x, row in df_sources.iterrows() if not isNaN(df_sources["MassBank"][x])]
                            index_HMDB = [x for x, row in df_sources.iterrows() if not isNaN(df_sources["HMDB"][x])]

                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            

                        #2 SGM
                        elif  len(sirius_df) > 0 and len(gnps_df) > 0 and len(mbank_df) > 0 and len(hmdb_df) == 0:

                            mbank_df["rank_ids"] = ["M_" + str(s+1) for s in range(len(mbank_df))]

                            gnps_df["rank_ids"] = ["G_" + str(s+1) for s in range(len(gnps_df))]

                            sirius_df["rank_ids"] = ["S_" + str(s) for s in sirius_df["rank"]]
                            sirius_df["Source"] = "SIRIUS"


                            source_l1 = [*(list(sirius_df["Source"])), 
                               *(list(gnps_df["Source"]))
                               ,*(list(mbank_df["Source"]))]

                            rank_l2 = [*(list(sirius_df["rank_ids"])), 
                               *(list(gnps_df["rank_ids"]))
                               ,*(list(mbank_df["rank_ids"]))]

                            smiles_l3 = [*(list(sirius_df["smiles"])), 
                               *(list(gnps_df["GNPSSMILES"]))
                               ,*(list(mbank_df["MBSMILES"]))]

                            sm = pd.DataFrame(list(zip(source_l1, rank_l2, smiles_l3)), columns = ["Source", "ranks", "SMILES"])


                            df_edge = chemMN_CandidateSelection(sm)

                            df_edge.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_ChemMNedges.tsv", sep='\t')
                            one_candidate = one_candidate_selection(sm, Source = "SGM")
                            one_candidate.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_one_candidate_list.tsv", sep='\t')

                        #3 SHM
                        elif  len(sirius_df) > 0 and len(gnps_df) == 0 and len(mbank_df) > 0 and len(hmdb_df) > 0:

                            mbank_df["rank_ids"] = ["M_" + str(s+1) for s in range(len(mbank_df))]

                            #gnps_df["rank_ids"] = ["G_" + str(s+1) for s in range(len(gnps_df))]

                            hmdb_df["rank_ids"] = ["H_" + str(s+1) for s in range(len(hmdb_df))]

                            sirius_df["rank_ids"] = ["S_" + str(s) for s in sirius_df["rank"]]
                            sirius_df["Source"] = "SIRIUS"


                            source_l1 = [*(list(sirius_df["Source"]))
                               ,*(list(mbank_df["Source"])),
                                        *(list(hmdb_df["Source"]))]

                            rank_l2 = [*(list(sirius_df["rank_ids"]))
                               ,*(list(mbank_df["rank_ids"])), 
                                      *(list(hmdb_df["rank_ids"]))]

                            smiles_l3 = [*(list(sirius_df["smiles"]))
                               ,*(list(mbank_df["MBSMILES"])),
                                        *(list(hmdb_df["HMDBSMILES"]))]

                            sm = pd.DataFrame(list(zip(source_l1, rank_l2, smiles_l3)), columns = ["Source", "ranks", "SMILES"])


                            df_edge = chemMN_CandidateSelection(sm)

                            df_edge.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_ChemMNedges.tsv", sep='\t')
                            one_candidate = one_candidate_selection(sm, Source = "SHM")
                            one_candidate.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_one_candidate_list.tsv", sep='\t')

                        #4 SGH
                        elif  len(sirius_df) > 0 and len(gnps_df) > 0 and len(mbank_df) == 0 and len(hmdb_df) > 0:
                            #mbank_df["rank_ids"] = ["M_" + str(s+1) for s in range(len(mbank_df))]

                            gnps_df["rank_ids"] = ["G_" + str(s+1) for s in range(len(gnps_df))]

                            hmdb_df["rank_ids"] = ["H_" + str(s+1) for s in range(len(hmdb_df))]

                            sirius_df["rank_ids"] = ["S_" + str(s) for s in sirius_df["rank"]]
                            sirius_df["Source"] = "SIRIUS"


                            source_l1 = [*(list(sirius_df["Source"]))
                               ,*(list(gnps_df["Source"])),
                                        *(list(hmdb_df["Source"]))]

                            rank_l2 = [*(list(sirius_df["rank_ids"]))
                               ,*(list(gnps_df["rank_ids"])), 
                                      *(list(hmdb_df["rank_ids"]))]

                            smiles_l3 = [*(list(sirius_df["smiles"]))
                               ,*(list(gnps_df["GNPSSMILES"])),
                                        *(list(hmdb_df["HMDBSMILES"]))]

                            sm = pd.DataFrame(list(zip(source_l1, rank_l2, smiles_l3)), columns = ["Source", "ranks", "SMILES"])

                            df_edge = chemMN_CandidateSelection(sm)

                            df_edge.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_ChemMNedges.tsv", sep='\t')
                            one_candidate = one_candidate_selection(sm, Source = "SGH")
                            one_candidate.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_one_candidate_list.tsv", sep='\t')
                            
                        #5 GHM
                        elif  len(sirius_df) == 0 and len(gnps_df) > 0 and len(mbank_df) > 0 and len(hmdb_df) > 0:
                            mbank_df["rank_ids"] = ["M_" + str(s+1) for s in range(len(mbank_df))]

                            gnps_df["rank_ids"] = ["G_" + str(s+1) for s in range(len(gnps_df))]

                            hmdb_df["rank_ids"] = ["H_" + str(s+1) for s in range(len(hmdb_df))]

                            #sirius_df["rank_ids"] = ["S_" + str(s) for s in sirius_df["rank"]]
                            #sirius_df["Source"] = "SIRIUS"


                            source_l1 = [*(list(gnps_df["Source"]))
                               ,*(list(mbank_df["Source"])),
                                        *(list(hmdb_df["Source"]))]

                            rank_l2 = [*(list(gnps_df["rank_ids"]))
                               ,*(list(mbank_df["rank_ids"])), 
                                      *(list(hmdb_df["rank_ids"]))]

                            smiles_l3 = [*(list(gnps_df["GNPSSMILES"]))
                               ,*(list(mbank_df["MBSMILES"])),
                                        *(list(hmdb_df["HMDBSMILES"]))]

                            sm = pd.DataFrame(list(zip(source_l1, rank_l2, smiles_l3)), columns = ["Source", "ranks", "SMILES"])

                            df_edge = chemMN_CandidateSelection(sm)

                            df_edge.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_ChemMNedges.tsv", sep='\t')
                            one_candidate = one_candidate_selection(sm, Source = "GHM")
                            one_candidate.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_one_candidate_list.tsv", sep='\t')
                            
                        #6 SG
                        elif  len(sirius_df) > 0 and len(gnps_df) > 0 and len(mbank_df) == 0 and len(hmdb_df) == 0:
                            #mbank_df["rank_ids"] = ["M_" + str(s+1) for s in range(len(mbank_df))]

                            gnps_df["rank_ids"] = ["G_" + str(s+1) for s in range(len(gnps_df))]

                            sirius_df["rank_ids"] = ["S_" + str(s) for s in sirius_df["rank"]]
                            sirius_df["Source"] = "SIRIUS"


                            source_l1 = [*(list(sirius_df["Source"])), 
                               *(list(gnps_df["Source"]))]
                               #,*(list(mbank_df["Source"]))]

                            rank_l2 = [*(list(sirius_df["rank_ids"])), 
                               *(list(gnps_df["rank_ids"]))]
                               #,*(list(mbank_df["rank_ids"]))]

                            smiles_l3 = [*(list(sirius_df["smiles"])), 
                               *(list(gnps_df["GNPSSMILES"]))]
                               #,*(list(mbank_df["MBSMILES"]))]

                            sm = pd.DataFrame(list(zip(source_l1, rank_l2, smiles_l3)), columns = ["Source", "ranks", "SMILES"])


                            df_edge = chemMN_CandidateSelection(sm)

                            df_edge.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_ChemMNedges.tsv", sep='\t')
                            
                            one_candidate = one_candidate_selection(sm, Source = "SG")
                            one_candidate.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_one_candidate_list.tsv", sep='\t')

                        #7 SH
                        elif  len(sirius_df) > 0 and len(gnps_df) == 0 and len(mbank_df) == 0 and len(hmdb_df) > 0:
                            #mbank_df["rank_ids"] = ["M_" + str(s+1) for s in range(len(mbank_df))]

                            #gnps_df["rank_ids"] = ["G_" + str(s+1) for s in range(len(gnps_df))]

                            hmdb_df["rank_ids"] = ["H_" + str(s+1) for s in range(len(hmdb_df))]

                            sirius_df["rank_ids"] = ["S_" + str(s) for s in sirius_df["rank"]]
                            sirius_df["Source"] = "SIRIUS"


                            source_l1 = [*(list(sirius_df["Source"])),
                                        *(list(hmdb_df["Source"]))]

                            rank_l2 = [*(list(sirius_df["rank_ids"])), 
                                      *(list(hmdb_df["rank_ids"]))]

                            smiles_l3 = [*(list(sirius_df["smiles"])),
                                        *(list(hmdb_df["HMDBSMILES"]))]

                            sm = pd.DataFrame(list(zip(source_l1, rank_l2, smiles_l3)), columns = ["Source", "ranks", "SMILES"])
                            df_edge = chemMN_CandidateSelection(sm)

                            df_edge.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_ChemMNedges.tsv", sep='\t')
                            
                            one_candidate = one_candidate_selection(sm, Source = "SH")
                            one_candidate.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_one_candidate_list.tsv", sep='\t')
                        #8 SM
                        elif  len(sirius_df) > 0 and len(gnps_df) == 0 and len(mbank_df) > 0 and len(hmdb_df) == 0:
                            mbank_df["rank_ids"] = ["M_" + str(s+1) for s in range(len(mbank_df))]

                            #gnps_df["rank_ids"] = ["G_" + str(s+1) for s in range(len(gnps_df))]

                            #hmdb_df["rank_ids"] = ["H_" + str(s+1) for s in range(len(hmdb_df))]

                            sirius_df["rank_ids"] = ["S_" + str(s) for s in sirius_df["rank"]]
                            sirius_df["Source"] = "SIRIUS"


                            source_l1 = [*(list(sirius_df["Source"]))
                               ,*(list(mbank_df["Source"]))]

                            rank_l2 = [*(list(sirius_df["rank_ids"]))
                               ,*(list(mbank_df["rank_ids"]))]

                            smiles_l3 = [*(list(sirius_df["smiles"]))
                               ,*(list(mbank_df["MBSMILES"]))]

                            sm = pd.DataFrame(list(zip(source_l1, rank_l2, smiles_l3)), columns = ["Source", "ranks", "SMILES"])


                            df_edge = chemMN_CandidateSelection(sm)

                            df_edge.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_ChemMNedges.tsv", sep='\t')
                            
                            one_candidate = one_candidate_selection(sm, Source = "SM")
                            one_candidate.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_one_candidate_list.tsv", sep='\t')
                        #9 GM
                        elif  len(sirius_df) == 0 and len(gnps_df) > 0 and len(mbank_df) > 0 and len(hmdb_df) == 0:
                            mbank_df["rank_ids"] = ["M_" + str(s+1) for s in range(len(mbank_df))]

                            gnps_df["rank_ids"] = ["G_" + str(s+1) for s in range(len(gnps_df))]

                            #hmdb_df["rank_ids"] = ["H_" + str(s+1) for s in range(len(hmdb_df))]

                            #sirius_df["rank_ids"] = ["S_" + str(s) for s in sirius_df["rank"]]
                            #sirius_df["Source"] = "SIRIUS"


                            source_l1 = [*(list(mbank_df["Source"])),
                                        *(list(gnps_df["Source"]))]

                            rank_l2 = [*(list(mbank_df["rank_ids"])), 
                                      *(list(gnps_df["rank_ids"]))]

                            smiles_l3 = [*(list(mbank_df["MBSMILES"])),
                                        *(list(gnps_df["GNPSSMILES"]))]

                            sm = pd.DataFrame(list(zip(source_l1, rank_l2, smiles_l3)), columns = ["Source", "ranks", "SMILES"])


                            df_edge = chemMN_CandidateSelection(sm)

                            df_edge.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_ChemMNedges.tsv", sep='\t')
                            
                            one_candidate = one_candidate_selection(sm, Source = "GM")
                            one_candidate.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_one_candidate_list.tsv", sep='\t')
                        #10 GH
                        elif  len(sirius_df) == 0 and len(gnps_df) > 0 and len(mbank_df) == 0 and len(hmdb_df) > 0:
                            #mbank_df["rank_ids"] = ["M_" + str(s+1) for s in range(len(mbank_df))]

                            gnps_df["rank_ids"] = ["G_" + str(s+1) for s in range(len(gnps_df))]

                            hmdb_df["rank_ids"] = ["H_" + str(s+1) for s in range(len(hmdb_df))]

                            #sirius_df["rank_ids"] = ["S_" + str(s) for s in sirius_df["rank"]]
                            #sirius_df["Source"] = "SIRIUS"


                            source_l1 = [*(list(gnps_df["Source"])),
                                        *(list(hmdb_df["Source"]))]

                            rank_l2 = [*(list(gnps_df["rank_ids"])),
                                      *(list(hmdb_df["rank_ids"]))]

                            smiles_l3 = [*(list(gnps_df["GNPSSMILES"])),
                                        *(list(hmdb_df["HMDBSMILES"]))]

                            sm = pd.DataFrame(list(zip(source_l1, rank_l2, smiles_l3)), columns = ["Source", "ranks", "SMILES"])
                            df_edge = chemMN_CandidateSelection(sm)



                            df_edge.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_ChemMNedges.tsv", sep='\t')
                            one_candidate = one_candidate_selection(sm, Source = "GH")
                            one_candidate.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_one_candidate_list.tsv", sep='\t')
                        #11 HM
                        elif  len(sirius_df) == 0 and len(gnps_df) == 0 and len(mbank_df) > 0 and len(hmdb_df) > 0:
                            mbank_df["rank_ids"] = ["M_" + str(s+1) for s in range(len(mbank_df))]

                            #gnps_df["rank_ids"] = ["G_" + str(s+1) for s in range(len(gnps_df))]

                            hmdb_df["rank_ids"] = ["H_" + str(s+1) for s in range(len(hmdb_df))]

                            #sirius_df["rank_ids"] = ["S_" + str(s) for s in sirius_df["rank"]]
                            #sirius_df["Source"] = "SIRIUS"


                            source_l1 = [*(list(mbank_df["Source"])),
                                        *(list(hmdb_df["Source"]))]

                            rank_l2 = [*(list(mbank_df["rank_ids"])), 
                                      *(list(hmdb_df["rank_ids"]))]

                            smiles_l3 = [*(list(mbank_df["MBSMILES"])),
                                        *(list(hmdb_df["HMDBSMILES"]))]

                            sm = pd.DataFrame(list(zip(source_l1, rank_l2, smiles_l3)), columns = ["Source", "ranks", "SMILES"])


                            df_edge = chemMN_CandidateSelection(sm)

                            df_edge.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_ChemMNedges.tsv", sep='\t')
                            one_candidate = one_candidate_selection(sm, Source = "HM")
                            one_candidate.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_one_candidate_list.tsv", sep='\t')

                        # S
                        #elif  len(sirius_df) > 0 and len(gnps_df) == 0 and len(mbank_df) == 0 and len(hmdb_df) == 0:

                        # G
                        #elif  len(sirius_df) == 0 and len(gnps_df) > 0 and len(mbank_df) == 0 and len(hmdb_df) == 0:

                        # M
                        #elif  len(sirius_df) == 0 and len(gnps_df) == 0 and len(mbank_df) > 0 and len(hmdb_df) == 0:

                        # H
                        #elif  len(sirius_df) == 0 and len(gnps_df) == 0 and len(mbank_df) == 0 and len(hmdb_df) > 0:
                merged_df.to_csv("mergedResults_withCandidates.csv")


# In[ ]:





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


# # Comparison with a list of SMILES from any Source

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




