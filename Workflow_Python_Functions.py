#!/usr/bin/env python
# coding: utf-8

# In[27]:


import pandas as pd
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
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
import pandas.io.formats.style
from rdkit.Chem import PandasTools
import xlrd
import openpyxl
import statistics


# ### SIRIUS post processing

# In[1]:


#input_table = pd.read_csv("/Users/mahnoorzulfiqar/Standards_CodeSet/ input_table.csv")
#input_table


# ### MetFrag Result Post Processing

# In[2]:


def metfrag_postproc(input_table, input_dir):

    for m, row in input_table.iterrows():
        result = input_table.loc[m, "ResultFileNames"] + "/insilico/MetFrag"
        files_met = (glob.glob(result+'/*.xls'))
        file1  = pd.read_csv(input_table['ResultFileNames'][m] + '/insilico/MS1DATA.csv')
        for i, row in file1.iterrows():
            pattern = file1.loc[i, "id_X"]
            results = [i for i in files_met if pattern in i]
            KEGG = [i for i in results if "KEGG" in i]
        
            try:
            
                KEGG_file = pd.read_excel(KEGG[0])
                if len(KEGG_file)>0:
            
                    KEGG_file = KEGG_file[KEGG_file['Score'].notna()]
                    KEGG_file = KEGG_file.drop(KEGG_file[KEGG_file.Score < 0.75].index)
        
                    file1.loc[i, 'KG_ID'] = KEGG_file.loc[0, 'Identifier']
                    file1.loc[i, 'KG_Name'] = KEGG_file.loc[0, 'CompoundName']
                    file1.loc[i, 'KG_Formula'] = KEGG_file.loc[0, 'MolecularFormula']
                    file1.loc[i, 'KG_expPeaks'] = KEGG_file.loc[0, 'NoExplPeaks']
                    file1.loc[i, 'KG_SMILES'] = Chem.MolToSmiles(Chem.MolFromInchi(KEGG_file["InChI"][0]))
                    file1.loc[i, 'KG_file'] = KEGG
                    Kegg_smiles = []
                    Kegg_s = []
                    for j in KEGG_file["InChI"][0:5].tolist():
                        mol = Chem.MolToSmiles(Chem.MolFromInchi(j))
                        Kegg_s.append(mol)
                        for k, row in sl.iterrows():
                            try:
                                SSms = [mol, Chem.MolFromSmiles(sl['SMILES'][k])]
                                SSfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=1024) for x in SSms]
                                SStn = DataStructs.FingerprintSimilarity(SSfps[0],SSfps[1])
                                if SStn >= 0.8:
                                    file1['KG_Top_can_SL'][i] = j
                                    file1['KG_tanimotoSLvsCAN'][i] = SStn
                                    file1['KG_SL_comp'][i] = sl['SMILES'][k]
                            except:
                                pass
                        sm = Chem.MolFromSmiles(mol)
                        Kegg_smiles.append(sm)
                    if len(Kegg_smiles) > 1:
                        res = rdFMCS.FindMCS(Kegg_smiles)
                        file1['KG_MCSSstring'][i] = res.smartsString
                        file1['KG_MCSS_SMILES'][i] = Chem.MolToSmiles(Chem.MolFromSmarts(res.smartsString))
                    
                #file1.at[i, 'KG_SMILESforMCSS'] = str(Kegg_s)
            except:
                pass
            PubChem = [i for i in results if "PubChem" in i]
            try:
            
                PubChem_file = pd.read_excel(PubChem[0])
        
                if len(PubChem_file)>0:
                    PubChem_file = PubChem_file[PubChem_file['Score'].notna()]
                    PubChem_file = PubChem_file.drop(PubChem_file[PubChem_file.Score < 0.75].index)
        
                    file1.loc[i, 'PC_ID'] = PubChem_file.loc[0, 'Identifier']
                    file1.loc[i, 'PC_Formula'] = PubChem_file.loc[0, 'MolecularFormula']
                    file1.loc[i, 'PC_expPeaks'] = PubChem_file.loc[0, 'NoExplPeaks']
                    file1.loc[i, 'PC_SMILES'] = PubChem_file.loc[0, 'SMILES']
                    file1.loc[i, 'PC_file'] = PubChem
            
                    Pubchem_smiles = []
        
                    for l in PubChem_file["SMILES"][0:5].tolist():
                        for n, row in sl.iterrows():
                            try:
                                SSms = [l, Chem.MolFromSmiles(sl['SMILES'][n])]
                                SSfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=1024) for x in SSms]
                                SStn2 = DataStructs.FingerprintSimilarity(SSfps[0],SSfps[1])
                                if SStn2 >= 0.8:
                                    file1['PC_Top_can_SL'][i] = l
                                    file1['PC_tanimotoSLvsCAN'][i] = SStn2
                                    file1['PC_SL_comp'][i] = sl['SMILES'][n]
                            except:
                                pass
                        sm2 = Chem.MolFromSmiles(l)
                        Pubchem_smiles.append(sm2)
                        if len(Pubchem_smiles) > 1:
                            res2 = rdFMCS.FindMCS(Pubchem_smiles)
                            file1['PC_MCSSstring'][i] = res2.smartsString
                            file1['PC_MCSS_SMILES'][i] = Chem.MolToSmiles(Chem.MolFromSmarts(res2.smartsString))
        
                    #file1.at[i, 'PC_SMILESforMCSS'] = PubChem_file["SMILES"][0:5].tolist()
            except:
                pass
    file1.to_csv(input_table['ResultFileNames'][m] + '/insilico/MetFragResults.csv')


# ### SIRIUS Result Post Processing

# In[3]:


def sirius_postProc2(input_table, input_dir):
    
    
    for m, row in input_table.iterrows():
        file1  = pd.read_csv(input_table['ResultFileNames'][m] + '/insilico/MS1DATAsirius.csv')
        file1  = pd.read_csv(input_table['ResultFileNames'][m] + '/insilico/MS1DATAsirius.csv')
        file1['MCSSstring'] = np.nan
        file1['MCSS_SMILES'] = np.nan
        file1['Top_can_SL'] = np.nan
        file1['tanimotoSLvsCAN'] = np.nan
        file1['SL_comp'] = np.nan
        for i, row in file1.iterrows():
            if not isNaN(file1['SMILESforMCSS'][i]):
                top_smiles = file1['SMILESforMCSS'][i].split("|")
                mols = []
                for j in top_smiles:
                    for k, row in sl.iterrows():
                        SSms = [Chem.MolFromSmiles(j), Chem.MolFromSmiles(sl['SMILES'][k])]
                        SSfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=1024) for x in SSms]
                        SStn = DataStructs.FingerprintSimilarity(SSfps[0],SSfps[1])
                        if SStn >= 0.8:
                            file1['Top_can_SL'][i] = j
                            file1['tanimotoSLvsCAN'][i] = SStn
                            file1['SL_comp'][i] = sl['SMILES'][k]
                    sm = Chem.MolFromSmiles(j)
                    mols.append(sm)
                if len(mols) > 1:
                    res = rdFMCS.FindMCS(mols)
                    file1['MCSSstring'][i] = res.smartsString
                    file1['MCSS_SMILES'][i] = Chem.MolToSmiles(Chem.MolFromSmarts(res.smartsString))
        file1.to_csv(input_table['ResultFileNames'][m] + '/insilico/SiriusResults.csv')


# ### COMBINE IN SILICO -All files with SIRIUS results separate and with MetFragresults separate

# In[4]:


def combine_insilico(input_table, Source = "SIRIUS"):
    
    path = os.path.join(input_dir, "MetabolomicsResults")
    if not os.path.isdir(path):
        os.mkdir(path)    
    
    if Source is "SIRIUS":
        
        all_files = []
        for n, row in input_table.iterrows():
            all_files.append(input_table['ResultFileNames'][n] + '/MetabolomicsResults/SiriusResults.csv')
        
        li = []
    
        for filename in all_files:
            df = pd.read_csv(filename, index_col=None, header=0)
            df["ResultFileNames"] = filename
            li.append(df)

        frame = pd.concat(li, axis=0, ignore_index=True)
        frame.to_csv(input_dir, '/SIRIUS_combined.csv')
        return(frame)
    
    
    elif Source is "MetFrag":
        all_files = []
        for m, row in input_table.iterrows():
            all_files.append(input_table['result_dir'][m] + '/insilico/MetFragResults.csv')
        li = []

        for filename in all_files:
            df = pd.read_csv(filename, index_col=None, header=0)
            df["result_dir"] = filename
            li.append(df)

        frame = pd.concat(li, axis=0, ignore_index=True)
        frame.to_csv(input_dir+'MetabolomicsResults/MetFrag_combined.csv')
        return(frame)


# ### Classification using Canopus and ClassyFire

# In[5]:


def classification(frame):
    inchis = []
    for l, row in frame.iterrows():
        if frame["Result"][l] == "SIRIUS_FOR":
            sep = 'json/'
            strpd = frame["dir"][l].split(sep, 1)[0] +"json/canopus_summary.tsv"
            if os.path.isfile(strpd):
                canopus = pd.read_csv(strpd, sep='\t')
                if len(canopus) > 0:
                    frame.loc[l, 'most_specific_class'] = canopus["most specific class"][0]
                    frame.loc[l, 'level _5'] = canopus["level 5"][0]
                    frame.loc[l, 'subclass'] = canopus["subclass"][0]
                    frame.loc[l, 'class'] = canopus["class"][0]
                    frame.loc[l, 'superclass'] = canopus["superclass"][0]
                    frame.loc[l, 'all_classifications'] = canopus["all classifications"][0]
                    frame.loc[l, 'class_source'] = "Canopus"
        elif frame["Result"][l] == "SIRIUS_STR":
            InChI = Chem.MolToInchi(Chem.MolFromSmiles(frame["SMILES"][l]))
            InChIKey = Chem.inchi.InchiToInchiKey(InChI)
            inchis.append({
                'index': l,
                'smiles':frame["SMILES"][l],
                'inchi': InChI,
                'inchikey': InChIKey
            })
    inchis = pd.DataFrame(inchis)
    
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
    
        start_time = time.time()
    
        print('%s inchikey to resolve' % total_inchikey_number )
        get_classifications_cf_mod(all_inchi_keys, par_level = 6)
    
        cleanse('all_json.json', 'all_json.json')
    
        with open("all_json.json") as tweetfile:
            jsondic = json.loads(tweetfile.read())

        df = json_normalize(jsondic)
        df = df.drop_duplicates( 'inchikey' )
        resolved_ik_number = len( df.drop_duplicates('inchikey').inchikey )
        resolved_ik_number_list.append( resolved_ik_number )
        print('%s resolved inchikeys' % resolved_ik_number )
        print("done in --- %s seconds ---" % (time.time() - start_time))
    
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
                frame.loc[q, 'class_source'] = "ClassyFire"
    frame.to_csv(input_dir, '/SIRIUS_combined.csv')
    return(frame)


# In[ ]:





# ### Spectral DB Dereplication Results Post Processing

# In[6]:


def isNaN(string):
    return string != string


# In[7]:


#input_dir = '/Users/mahnoorzulfiqar/Standards_CodeSet/'


# In[8]:


#hmdb_sdf = '/Users/mahnoorzulfiqar/OneDriveUNI/S_CResults/structures.sdf'


# In[9]:


#dframe = PandasTools.LoadSDF(hmdb_sdf,
                                #idName='HMDB_ID',
                                #smilesName='SMILES',
                                #molColName='Molecule',
                                #includeFingerprints=False)


# In[10]:


#dframe.columns


# In[ ]:





# In[11]:


#input_dir = os.getcwd() +"/"


# ### GNPS, MassBank and HMDB Results post processing

# In[12]:


def spec_postproc(input_dir, dframe):
    
    GNPScsvfiles = []
    HMDBcsvfiles = []
    MassBankcsvfiles = []

    #list all files and directories
    for entry in os.listdir(input_dir):
        if os.path.isdir(os.path.join(input_dir, entry)):
            sub_dir = input_dir + entry + '/spectral_dereplication'
            if os.path.exists(sub_dir):
                files = (glob.glob(sub_dir+'/*.csv'))
                #if len(files) is 3:
                for f in files:
                    if 'gnps.' in f: 
                        GNPScsvfiles.append(f)
                    if 'hmdb.' in f: 
                        HMDBcsvfiles.append(f)
                    if 'mbank.' in f: 
                        MassBankcsvfiles.append(f)
                            
    

    #### read sdf file from HMDB to collect names and smiles ####
    
    #HMDB CSV Result file pre_processing

    HMDBcsvfiles2 = []
    for k in HMDBcsvfiles:
    
        # read the cs files
        hmdb_df = pd.read_csv(k)
    
        #merge on basis of id, frame and hmdb result files
        SmilesHM = pd.merge(hmdb_df, dframe, left_on=hmdb_df.HMDBcompoundID, right_on=dframe.DATABASE_ID)
    
        for i, row in hmdb_df.iterrows():
            for j, row in SmilesHM.iterrows():
                #where index for both match, add the name and SMILES
                if hmdb_df['id_X'][i]== SmilesHM['id_X'][j]:
                    hmdb_df.loc[i, 'HMDBSMILES'] = SmilesHM['SMILES'][j]#add SMILES
                    hmdb_df.loc[i, 'HMDBcompound_name'] = SmilesHM["GENERIC_NAME"][j]#add name
                    hmdb_df.loc[i, 'HMDBformula'] = SmilesHM["FORMULA"][j]
                
        csvname = (os.path.splitext(k)[0])+"proc"+".csv" # name for writing it in a new file
        hmdb_df.to_csv(csvname) #write
        HMDBcsvfiles2.append(csvname)# add to a list
        
        
    MassBankcsvfiles2 = []
    for l in MassBankcsvfiles:
        # read mbank csv file
        mbank_df = pd.read_csv(l)
        for i, row in mbank_df.iterrows():
            try:
                #if not isNaN(mbank_df['MBinchiKEY'][i]):
                inchiK = str(mbank_df["MBinchiKEY"][i])#name from suspect list
                y = pcp.get_compounds(inchiK, 'inchikey')#compound based on inchikey
                for compound in y:
                    smles = compound.isomeric_smiles
                    #print(smles)
                    mbank_df.loc[i, 'MBSMILES'] = smles
            except:
                print("error")
        csvname = (os.path.splitext(l)[0])+"proc"+".csv"
        mbank_df.to_csv(csvname)
        MassBankcsvfiles2.append(csvname)
    
    
    matches = ["M+","[M", "M-", "2M", "M*" "20.0", "50.0", "30.0", "40.0", "60.0", "70.0", "eV", "Massbank"
               , "Spectral", "Match", "to", "from", "NIST14", "MoNA", '[IIN-based:',  '[IIN-based', 'on:', 'CCMSLIB00003136269]']

    GNPScsvfiles2 = []
    for l in GNPScsvfiles:
        gnps_df = pd.read_csv(l)
        gnps_df['corr_names'] = np.nan
        for i, row in gnps_df.iterrows():
            if not isNaN(gnps_df['GNPScompound_name'][i]):
                string_chng = (gnps_df['GNPScompound_name'][i].split(" "))
                newstr = []
                # for each part of the string in the names
                chng = []
                for j in range(len(string_chng)):
                    if not any(x in string_chng[j] for x in matches) and not '-' == string_chng[j]:
                        if '|' not in string_chng[j] or '!' not in string_chng[j]:
                            newstr.append(string_chng[j])
                        elif '|' in string_chng[j]:
                            jlen = string_chng[j].split("|")
                            lst = len(jlen)-1
                            chng.append(jlen[lst])
                            break
                chng.append(' '.join(newstr))
                gnps_df["corr_names"][i] = chng[0]
                if isNaN(gnps_df['GNPSSMILES'][i]):
                    try:
                        if chng is '':
                            break
                        else:
                            s = pcp.get_compounds(chng[0], 'name')
                            if s:
                                for comp in s:
                                    gnps_df["GNPSSMILES"][i] = comp.isomeric_smiles
                            else:
                                gnps_df["GNPSSMILES"][i] = ''
                    except Exception:
                        pass
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
    
    dict1 = {'GNPSr': GNPScsvfiles2, 'HMDBr': HMDBcsvfiles2, 'MBr': MassBankcsvfiles2} 
    df = pd.DataFrame(dict1)

    return(df)


# In[13]:


#spec_df = spec_postproc(input_dir, dframe)


# ### Combine_all Spectral DBs

# In[14]:


def combine_specdb(df):
    Merged_Result_df = []
    for i, row in df.iterrows():
        CSVfileG = pd.read_csv(df["GNPSr"][i])
        CSVfileH = pd.read_csv(df["HMDBr"][i])
        CSVfileM = pd.read_csv(df["MBr"][i])

        MergedRE = CSVfileG.merge(CSVfileH,on='premz').merge(CSVfileM,on='premz')

    
        csvname = (os.path.splitext(df["GNPSr"][i])[0]).replace("gnps_csv_with_cor_names", "mergedR.csv")
        MergedRE.to_csv(csvname)
        Merged_Result_df.append(csvname)
    return(Merged_Result_df)


# In[15]:


#Merged_Result_df = combine_specdb(spec_df)


# ### Combine all files for spectral db dereplication

# In[16]:


def combine_allspec(comb_df):
    combined_csv = pd.concat([pd.read_csv(f) for f in comb_df], ignore_index=True)
    for i, row in combined_csv.iterrows():
        if combined_csv['GNPSSMILES'][i] == ' ' or isNaN(combined_csv['GNPSSMILES'][i]):
            combined_csv['GNPSSMILES'][i] = ''
    for i, row in combined_csv.iterrows():
        if not isNaN(combined_csv['MBinchiKEY'][i]):
            try:
                y = pcp.get_compounds(combined_csv['MBinchiKEY'][i], 'inchikey')
                if len(y)>1:
                    combined_csv['MBSMILES'][i] = y[0].isomeric_smiles
            except:
                pass
    combined_csv.to_csv(input_dir + 'MetabolomicsResults/SD_post_processed_combined_results.csv')
    return(combined_csv)


# In[17]:


#combined = combine_allspec(Merged_Result_df)


# ### Scoring Scheme for Spectral DB Dereplication

# In[18]:


def HMDB_Scoring(db, i):
    if db['HMDBmax_similarity'][i] >= 0.75 and db['HMDBintScore'][i] >= 0.50 and db['HMDBmzScore'][i] >= 0.50 and db['HMDBQMatchingPeaks'][i]/db['hmdbQueryTotalPeaks'][i] >= 0.50:
        return True
    else:
        return False


# In[19]:


def GNPS_Scoring(db, i):
    if db['GNPSmax_similarity'][i] >= 0.75 and db['GNPSintScore'][i] >= 0.50 and db['GNPSmzScore'][i] >= 0.50 and db['GQMatchingPeaks'][i]/db['gQueryTotalPeaks'][i] >= 0.50:
        return True
    else:
        return False


# In[20]:


def MB_Scoring(db, i):
    if db['MBmax_similarity'][i] >= 0.75 and db['MBintScore'][i] >= 0.50 and db['MBmzScore'][i] >= 0.50 and db['MQMatchingPeaks'][i]/db['mQueryTotalPeaks'][i] >= 0.50:
        return True
    else:
        return False


# In[ ]:





# In[21]:


def scoring_spec(combined):
    
    for i, row in combined.iterrows():
        if HMDB_Scoring(combined, i) and GNPS_Scoring(combined, i) and MB_Scoring(combined, i) and not isNaN(combined['GNPSSMILES'][i]) and not isNaN(combined['MBSMILES'][i]) and not isNaN(combined['HMDBSMILES'][i]):
            HGms = [Chem.MolFromSmiles(combined['HMDBSMILES'][i]), Chem.MolFromSmiles(combined['GNPSSMILES'][i])]
            HGfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=1024) for x in HGms]
            HGtn = DataStructs.FingerprintSimilarity(HGfps[0],HGfps[1])
            GMms = [Chem.MolFromSmiles(combined['GNPSSMILES'][i]), Chem.MolFromSmiles(combined['MBSMILES'][i])]
            GMfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=1024) for x in GMms]
            GMtn = DataStructs.FingerprintSimilarity(GMfps[0],GMfps[1])
            HMms = [Chem.MolFromSmiles(combined['HMDBSMILES'][i]), Chem.MolFromSmiles(combined['MBSMILES'][i])]
            HMfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=1024) for x in HMms]
            HMtn = DataStructs.FingerprintSimilarity(HMfps[0],HMfps[1])
            
            combined.loc[i, 'annotation'] = 'HMDB, GNPS, MassBank'
            combined.loc[i, 'tanimotoHG'] = HGtn
            combined.loc[i, 'tanimotoGM'] = GMtn
            combined.loc[i, 'tanimotoHM'] = HMtn
            combined.loc[i, 'occurance'] = 3
        
        if HMDB_Scoring(combined, i) and GNPS_Scoring(combined, i) and not MB_Scoring(combined, i) and not isNaN(combined['GNPSSMILES'][i]) and not isNaN(combined['HMDBSMILES'][i]):
            HGms = [Chem.MolFromSmiles(combined['HMDBSMILES'][i]), Chem.MolFromSmiles(combined['GNPSSMILES'][i])]
            HGfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=1024) for x in HGms]
            HGtn = DataStructs.FingerprintSimilarity(HGfps[0],HGfps[1])
        
            combined.loc[i, 'annotation'] = 'HMDB, GNPS'
            combined.loc[i, 'tanimotoHG'] = HGtn
            combined.loc[i, 'tanimotoGM'] = np.nan
            combined.loc[i, 'tanimotoHM'] = np.nan
            combined.loc[i, 'occurance'] = 2
        
        if not HMDB_Scoring(combined, i) and GNPS_Scoring(combined, i) and MB_Scoring(combined, i) and not isNaN(combined['MBSMILES'][i]) and not isNaN(combined['GNPSSMILES'][i]):
            GMms = [Chem.MolFromSmiles(combined['GNPSSMILES'][i]), Chem.MolFromSmiles(combined['MBSMILES'][i])]
            GMfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=1024) for x in GMms]
            GMtn = DataStructs.FingerprintSimilarity(GMfps[0],GMfps[1])
        
            combined.loc[i, 'annotation'] = 'GNPS, MassBank'
            combined.loc[i, 'tanimotoHG'] = np.nan
            combined.loc[i, 'tanimotoGM'] = GMtn
            combined.loc[i, 'tanimotoHM'] = np.nan
            combined.loc[i, 'occurance'] = 2
        
        if HMDB_Scoring(combined, i) and not GNPS_Scoring(combined, i) and MB_Scoring(combined, i) and not isNaN(combined['MBSMILES'][i]) and not isNaN(combined['HMDBSMILES'][i]):
            HMms = [Chem.MolFromSmiles(combined['HMDBSMILES'][i]), Chem.MolFromSmiles(combined['MBSMILES'][i])]
            HMfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=1024) for x in HMms]
            HMtn = DataStructs.FingerprintSimilarity(HMfps[0],HMfps[1])
        
            combined.loc[i, 'annotation'] = 'HMDB, MassBank'
            combined.loc[i, 'tanimotoHG'] = np.nan
            combined.loc[i, 'tanimotoGM'] = np.nan
            combined.loc[i, 'tanimotoHM'] = HMtn
            combined.loc[i, 'occurance'] = 2
        
        if HMDB_Scoring(combined, i) and not GNPS_Scoring(combined, i) and not MB_Scoring(combined, i):
        
            combined.loc[i, 'annotation'] = 'HMDB'
            combined.loc[i, 'tanimotoHG'] = np.nan
            combined.loc[i, 'tanimotoGM'] = np.nan
            combined.loc[i, 'tanimotoHM'] = np.nan
            combined.loc[i, 'occurance'] = 1
        
        if not HMDB_Scoring(combined, i) and GNPS_Scoring(combined, i) and not MB_Scoring(combined, i):
        
            combined.loc[i, 'annotation'] = 'GNPS'
            combined.loc[i, 'tanimotoHG'] = np.nan
            combined.loc[i, 'tanimotoGM'] = np.nan
            combined.loc[i, 'tanimotoHM'] = np.nan
            combined.loc[i, 'occurance'] = 1
        
        if not HMDB_Scoring(combined, i) and not GNPS_Scoring(combined, i) and MB_Scoring(combined, i):
        
            combined.loc[i, 'annotation'] = 'MassBank'
            combined.loc[i, 'tanimotoHG'] = np.nan
            combined.loc[i, 'tanimotoGM'] = np.nan
            combined.loc[i, 'tanimotoHM'] = np.nan
            combined.loc[i, 'occurance'] = 1
        
        if not HMDB_Scoring(combined, i) and not GNPS_Scoring(combined, i) and not MB_Scoring(combined, i):
            combined.loc[i, 'annotation'] = 'none'
            combined.loc[i, 'tanimotoHG'] = np.nan
            combined.loc[i, 'tanimotoGM'] = np.nan
            combined.loc[i, 'tanimotoHM'] = np.nan
            combined.loc[i, 'occurance'] = 0
    combined.to_csv(input_dir + "MetabolomicsResults/combinedSpecDB.csv")
    return(combined)


# In[ ]:





# In[22]:


#SpectralDB_Results = scoring_spec(combined)


# In[ ]:





# ### Suspect List Screening

# In[24]:


#Suspect_list = pd.read_csv('/Users/mahnoorzulfiqar/OneDriveUNI/SuspectList/Use_This_CURATED_SUSPECT_LIST.csv')


# In[25]:


def suspectListScreening(Suspect_list, SpectralDB_Results):
    
    SpectralDB_Results['HLsmiles'] = np.nan
    SpectralDB_Results['HLname'] = np.nan
    SpectralDB_Results['GLsmiles'] = np.nan
    SpectralDB_Results['GLname'] = np.nan
    SpectralDB_Results['MLsmiles'] = np.nan
    SpectralDB_Results['MLname'] = np.nan
    
    for i, row in SpectralDB_Results.iterrows():
        if not isNaN(SpectralDB_Results['HMDBSMILES'][i]) and SpectralDB_Results['HMDBSMILES'][i] is not " ":
            for j, row in Suspect_list.iterrows():
                LHms2 = [Chem.MolFromSmiles(SpectralDB_Results['HMDBSMILES'][i]), Chem.MolFromSmiles(Suspect_list['SMILES'][j])]
                LHfps2 = [AllChem.GetMorganFingerprintAsBitVect(x2,2, nBits=1024) for x2 in LHms2]
                LHtn2 = DataStructs.FingerprintSimilarity(LHfps2[0],LHfps2[1])
                if LHtn2 == 1:
                    SpectralDB_Results.loc[i, 'HLsmiles'] = Suspect_list['SMILES'][j]
                    SpectralDB_Results.loc[i, 'HLname'] = Suspect_list['Name'][j]
        if not isNaN(SpectralDB_Results['GNPSSMILES'][i]) and SpectralDB_Results['GNPSSMILES'][i] is not " ":
            for k, row in Suspect_list.iterrows():
                LGms2 = [Chem.MolFromSmiles(SpectralDB_Results['GNPSSMILES'][i]), Chem.MolFromSmiles(Suspect_list['SMILES'][k])]
                LGfps2 = [AllChem.GetMorganFingerprintAsBitVect(x2,2, nBits=1024) for x2 in LGms2]
                LGtn2 = DataStructs.FingerprintSimilarity(LGfps2[0],LGfps2[1])
                if LGtn2 == 1:
                    SpectralDB_Results.loc[i, 'GLsmiles'] = Suspect_list['SMILES'][k]
                    SpectralDB_Results.loc[i, 'GLname'] = Suspect_list['Name'][k]
        if not isNaN(SpectralDB_Results['MBSMILES'][i]) and SpectralDB_Results['MBSMILES'][i] is not " ":
            for l, row in Suspect_list.iterrows():
                LMms2 = [Chem.MolFromSmiles(SpectralDB_Results['MBSMILES'][i]), Chem.MolFromSmiles(Suspect_list['SMILES'][l])]
                LMfps2 = [AllChem.GetMorganFingerprintAsBitVect(x2,2, nBits=1024) for x2 in LMms2]
                LMtn2 = DataStructs.FingerprintSimilarity(LMfps2[0],LMfps2[1])
                if LMtn2 == 1:
                    SpectralDB_Results.loc[i, 'MLsmiles'] = Suspect_list['SMILES'][l]
                    SpectralDB_Results.loc[i, 'MLname'] = Suspect_list['Name'][l]
                    #SpectralDB_Results['occurance'][i] = SpectralDB_Results['occurance'][i] + 1
                
    for i, row in SpectralDB_Results.iterrows():
        if not isNaN(SpectralDB_Results['HLname'][i]) or not isNaN(SpectralDB_Results['GLname'][i]) or not isNaN(SpectralDB_Results['MLname'][i]):
            SpectralDB_Results['occurance'][i] = SpectralDB_Results['occurance'][i] + 1
            SpectralDB_Results['annotation'][i] = SpectralDB_Results['annotation'][i] + ', Suspect_List'
            
    combined.to_csv(input_dir + "MetabolomicsResults/SpecDBvsSL.csv")
    return(SpectralDB_Results)


# In[26]:


#suspectListScreening(Suspect_list, SpectralDB_Results)


# In[ ]:





# In[ ]:





# In[ ]:




