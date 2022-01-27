#!/usr/bin/env python
# coding: utf-8

# In[147]:


import pandas as pd
import numpy as np
def isNaN(string):
    return string != string
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS


# In[148]:


#Import Libraries
import pubchempy as pcp
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import PandasTools
from rdkit.Chem import rdFMCS
import pandas.io.formats.style
import os
import re
import glob
import xlrd
import openpyxl
import statistics
import time


# In[149]:


from pybatchclassyfire import *
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


# In[150]:


#df_excel = pd.read_csv('checkCD_Results_QC.csv')


# ## MetFrag Curation

# In[151]:


metfrag = pd.read_csv("/Users/mahnoorzulfiqar/OneDriveUNI/MZML/MetabolomicsResults/MetFrag_combined.csv")


# In[152]:


metfrag.drop(["Unnamed: 0", "Unnamed: 0.1", "Unnamed: 0.1.1",
             "int", "col_eng", "pol", "ms2Peaks", "ms1Peaks", "KG_ID",
             'KG_Formula', 'KG_file','KG_MCSSstring', 'PC_ID',
             'PC_Formula','PC_file', 'PC_MCSSstring', 'KG_Name',
             'result_dir'], axis = 1)


# top candidates have a score of 1.0 in MetFrag. 
# So for selection of candidate, use: expl_peaks and whether the result from KEGG nad PubChem matches.
# 
# condition:
# 1. for compounds with only Pubchem, or only KEGG, will remain as it is.
# 2. if there is KEGG and PubChem, calculate a tanimoto score, if it is equal to 1, keep PubChem entry but add KEGG to source add well.
# 3. if there is less than 1 tanimoto, Check whether the entry has a suspcte list compound, add the entry eith Suspct list , if no suspectlist, then add PubChem 

# In[153]:


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
            
            #if there is an entry from suspect list WITH kegg
            elif not isNaN(metfrag['KG_SL_comp'][i]) and isNaN(metfrag['PC_SL_comp'][i]):
                metfrag.loc[i, 'Annotation_M'] = 'KEGG, SuspectList'
                #metfrag.loc[i, 'SMILES_final'] = metfrag['KG_SMILES'][i]
                
            #if there is an entry from suspect list WITH kegg
            elif not isNaN(metfrag['PC_SL_comp'][i]) and isNaN(metfrag['KG_SL_comp'][i]):
                metfrag.loc[i, 'Annotation_M'] = 'PubChem, SuspectList'
                #metfrag.loc[i, 'SMILES_final'] = metfrag['PC_SMILES'][i]


# In[154]:


metfrag.to_csv("metfrag_curated.csv")


# ## SIRIUS Curation

# In[155]:


sirius = pd.read_csv("/Users/mahnoorzulfiqar/OneDriveUNI/MZML/MetabolomicsResults/SIRIUS_combined.csv")


# In[156]:


sirius.columns


# In[157]:


sirius.drop(["Unnamed: 0", "Unnamed: 0.1", "Unnamed: 0.1.1", "Unnamed: 0.1.1.1",
             "int", "col_eng", "pol", "ms2Peaks", "ms1Peaks", "rtmed", "Adducts", "Formula", "all_classifications",
             "subclass", "level _5", "most_specific_class", "PubChemIDs", "premz", "X", "dir",
             'result_dir'], axis = 1)


# important scores to consider are
# 1. SIRIUSscore
# 2. CSIFingerIDscore
# 3. exp_int >= 0.70
# 4. SL_comp (if present, give priority and mention in the annotation)

# In[158]:


for i, row in sirius.iterrows():
    
    # If the explained intensity is greater than 0.70 and there is no suspect list entry
    if sirius['exp_int'][i] >= 0.70 and isNaN(sirius['SL_comp'][i]):
        sirius.loc[i, 'Annotation_S'] = 'SIRIUS'
        #sirius.loc[i, 'SMILES_final'] = sirius['SMILES'][i]
    
    #If the explained intensity is greater than 0.70 and there is an entry from suspect list
    elif sirius['exp_int'][i] >= 0.70 and not isNaN(sirius['SL_comp'][i]):
        sirius.loc[i, 'Annotation_S'] = 'SIRIUS, SuspectList'
        #sirius.loc[i, 'SMILES_final'] = sirius['SMILES'][i]
    
    # if the intensity is less thna 0.70 but it still is similar to an entry in Suspect list,
    elif sirius['exp_int'][i] < 0.70 and not isNaN(sirius['SL_comp'][i]):
        sirius.loc[i, 'Annotation_S'] = 'SIRIUS, SuspectList'
        #sirius.loc[i, 'SMILES_final'] = sirius['SMILES'][i]


# In[159]:


sirius.to_csv("sirius_curated.csv")


# In[ ]:





# In[160]:


S_M_CSV = pd.concat([sirius, metfrag], axis = 1, levels = ["id_X"])


# In[161]:


S_M_CSV.columns


# In[162]:


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
                    if PStn is 1:
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
                    if SKtn is 1:
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


# ## Spectral Database Dereplication Curation

# In[163]:


combined = pd.read_csv("/Users/mahnoorzulfiqar/OneDriveUNI/MZML/MetabolomicsResults/SDB_SDvsSLvsCD.csv")


# In[164]:


combined.columns


# In[165]:


combined.drop(["Unnamed: 0", 
             "Unnamed: 0.1", 
             "Unnamed: 0.1.1", 
             "Unnamed: 0.1.1.1", 
             'Unnamed: 0_x', 
             'Unnamed: 0.1_x',
             'rtmin_x', 
             'rtmax_x',
             'GNPSmzScore', 
             'GNPSintScore', 
             'GQMatchingPeaks', 
             'gQueryTotalPeaks', 
             'GNPSTotalPeaks', 
             'GNPSresult_dir', 
             'GNPSSMILES', 
             'GNPSspectrumID', 
             'GNPSresult_dir',
             'Results_x', 
             'corr_names', 
             'polarity_x', 
             'Unnamed: 0_y', 
             'polarity_y', 
             'Unnamed: 0.1.1.1.1', 
             'Unnamed: 0.1.1.1.1.1', 
             'index',
             'rtmin', 
             'rtmax',
             'Unnamed: 0.1_y', 
             'index_y', 
             'rtmin_y', 
             'rtmax_y', 
             'HMDBmzScore', 
             'HMDBintScore', 
             'HMDBQMatchingPeaks',
             'hmdbQueryTotalPeaks', 
             'HMDBTotalPeaks', 
             'HMDBcompound_id', 
             'HMDBresult_dir', 
             'Results_y', 
             'HMDBSMILES', 
             'MBmzScore', 
             'MBintScore',
             'MBQMatchingPeaks', 
             'mbQueryTotalPeaks', 
             'MBTotalPeaks', 
             'MBinchiKEY', 
             'MBspectrumID', 
             'MBresult_dir',
             'Results', 
             'MBSMILES', 
             'polarity', 
             'rt_med', 
             'annotation', 
             'occurance', 
             'HLsmiles', 
             'HLname',
             'GLsmiles', 
             'GLname', 
             'MLsmiles', 
             'MLname', 
             'HCname', 
             'HCsmiles',
             'HCannotation', 
             'GCname', 
             'GCsmiles', 
             'GCannotation', 
             'MCname',
             'MCsmiles', 
             'MCannotation', 
             'Occurence'], axis = 1)


# In[166]:


def HMDB_Scoring(db, i):
    if db['HMDBmax_similarity'][i] >= 0.75 and db['HMDBintScore'][i] >= 0.50 and db['HMDBmzScore'][i] >= 0.50 and db['HMDBQMatchingPeaks'][i]/db['hmdbQueryTotalPeaks'][i] >= 0.50:
        return True
    else:
        return False
def GNPS_Scoring(db, i):
    if db['GNPSmax_similarity'][i] >= 0.75 and db['GNPSintScore'][i] >= 0.50 and db['GNPSmzScore'][i] >= 0.50 and db['GQMatchingPeaks'][i]/db['gQueryTotalPeaks'][i] >= 0.50:
        return True
    else:
        return False
def MB_Scoring(db, i):
    if db['MBmax_similarity'][i] >= 0.75 and db['MBintScore'][i] >= 0.50 and db['MBmzScore'][i] >= 0.50 and db['MBQMatchingPeaks'][i]/db['mbQueryTotalPeaks'][i] >= 0.50:
        return True
    else:
        return False


# In[167]:


for i, row in combined.iterrows():
    if HMDB_Scoring(combined, i) or GNPS_Scoring(combined, i) or MB_Scoring(combined, i):
        print(i)
    else:
        combined['GNPSspectrumID'][i] = np.nan
        combined['MBspectrumID'][i] = np.nan
        combined['HMDBcompound_id'][i] = np.nan


# In[ ]:





# In[168]:


for i, row in combined.iterrows():
    
    ##### When there is an annotaion from all DBs #####
    
    #all entries with a high scoring annotation in all DBs,
    if not isNaN(combined['GNPSspectrumID'][i]) and not isNaN(combined['MBspectrumID'][i]) and not isNaN(combined['HMDBcompound_id'][i]):
        
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
                print(combined['GLsmiles'][i])
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
    if not isNaN(combined['GNPSspectrumID'][i]) and isNaN(combined['MBspectrumID'][i]) and not isNaN(combined['HMDBcompound_id'][i]):

        if not isNaN(combined['GLname'][i]):
            combined.loc[i, 'Annotation'] = 'GNPS, SuspectList'
        elif not isNaN(combined['HLname'][i]):
            combined.loc[i, 'Annotation'] = 'HMDB, SuspectList'
        elif not isNaN(combined['GNPSSMILES'][i]):
            combined.loc[i, 'Annotation'] = 'GNPS'
        else:
            combined.loc[i, 'Annotation'] = 'HMDB'
        
    
    # only GNPS and MassBank
    if not isNaN(combined['GNPSspectrumID'][i]) and not isNaN(combined['MBspectrumID'][i]) and isNaN(combined['HMDBcompound_id'][i]):

        if not isNaN(combined['GLname'][i]):
            combined.loc[i, 'Annotation'] = 'GNPS, SuspectList'
        elif not isNaN(combined['MLname'][i]):
            combined.loc[i, 'Annotation'] = 'MassBank, SuspectList'
        elif not isNaN(combined['GNPSSMILES'][i]):
            combined.loc[i, 'Annotation'] = 'GNPS'
        else:
            combined.loc[i, 'Annotation'] = 'MassBank'
    
    # only MassBank and HMDB
    
    if isNaN(combined['GNPSspectrumID'][i]) and not isNaN(combined['MBspectrumID'][i]) and not isNaN(combined['HMDBcompound_id'][i]):

        if not isNaN(combined['MLname'][i]):
            combined.loc[i, 'Annotation'] = 'MassBank, SuspectList'
        elif not isNaN(combined['HLname'][i]):
            combined.loc[i, 'Annotation'] = 'HMDB, SuspectList'
        else:
            combined.loc[i, 'Annotation'] = 'MassBank'
    
    
    
    ##### When there is an annotation from one DBs #####
    
    
    # only GNPS
    if not isNaN(combined['GNPSspectrumID'][i]) and isNaN(combined['MBspectrumID'][i]) and isNaN(combined['HMDBcompound_id'][i]):
        
        #If also SuspectList
        if not isNaN(combined['GLname'][i]):
            combined.loc[i, 'Annotation'] = 'GNPS, SuspectList'
        elif not isNaN(combined['GNPSSMILES'][i]):
            combined.loc[i, 'Annotation'] = 'GNPS'
        
    # only MassBank
    if isNaN(combined['GNPSspectrumID'][i]) and not isNaN(combined['MBspectrumID'][i]) and isNaN(combined['HMDBcompound_id'][i]):
        combined.loc[i, 'Annotation'] = 'MassBank'
        #If also SuspectList
        if not isNaN(combined['MLname'][i]):
            combined.loc[i, 'Annotation'] = 'MassBank, SuspectList'
    
    # only HMDB
    if isNaN(combined['GNPSspectrumID'][i]) and isNaN(combined['MBspectrumID'][i]) and not isNaN(combined['HMDBcompound_id'][i]):
        combined.loc[i, 'Annotation'] = 'HMDB'
        #If also SuspectList
        if not isNaN(combined['HLname'][i]):
            combined.loc[i, 'Annotation'] = 'HMDB, SuspectList'
            


# ### combine spectral dbs and compound dbs

# In[169]:


mega = pd.concat([S_M_CSV, combined], axis = 1, levels = ["id_X"])


# In[170]:


list(mega.columns)


# In[171]:


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
                    print(SKtn)
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
                    print(SKtn)
                    mega.loc[i, "Annotation_Source"] = mega['Annotation'][i] + ', ' + mega['Annotation_C'][i]
                else:
                    mega.loc[i, "Annotation_Source"] = mega['Annotation'][i]
            
            #for GNPS also check if there are any smiles given
            elif 'GNPS' in mega['Annotation'][i] and not isNaN(mega['GNPSSMILES'][i]):
                SKms = [Chem.MolFromSmiles(mega['GNPSSMILES'][i]), Chem.MolFromSmiles(mega['SMILES'][i])]
                SKfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in SKms]
                SKtn = DataStructs.FingerprintSimilarity(SKfps[0],SKfps[1])
                if SKtn is 1:
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
            # do something
            mega.loc[i, "Annotation_Source"] = mega['Annotation'][i]
        
        
        
        #######ONE SDB SOURCE#########
        
        #if only one sdb, but no SIRIUS
        #prioritize Spectral DB
        elif len(mega['Annotation'][i].split()) == 1 and 'SIRIUS' not in mega['Annotation_C'][i]:
            # do something
            mega.loc[i, "Annotation_Source"] = mega['Annotation'][i]
        
        #if only one sdb
        #prioritize SIRIUS
        elif len(mega['Annotation'][i].split()) == 1 and 'SIRIUS' in mega['Annotation_C'][i]:
            mega.loc[i, "Annotation_Source"] = mega['Annotation_C'][i]
            
            
    # if only spectral db results
    if isNaN(mega['Annotation'][i]) and isNaN(mega['Annotation_C'][i]) and not isNaN(mega['Formula'][i]):
        mega.loc[i, "Annotation_Source"] = 'SIRIUS_Formula'


# In[172]:


possibilities = list(np.unique(list(mega['Annotation_Source'])))
possibilities


# In[173]:


list(mega.columns)


# In[ ]:





# In[174]:


megaA = mega[['id_X', 'premz', 'rtmed', 'int', 'col_eng', 'pol','Adducts', 'SMILES', 'Formula', 'SIRIUSscore',
'CSIFingerIDscore', 'MCSS_SMILES', 'most_specific_class', 'level _5','subclass', 'class', 'superclass', 'all_classifications', 
             'KG_SMILES', 'KG_Formula', 'KG_MCSS_SMILES',  'HLsmiles',
 'HLname',
 'GLsmiles',
 'GLname',
 'MLsmiles',
 'MLname',
'PC_SMILES', 'PC_Formula', 'PC_MCSS_SMILES', 'GNPSSMILES', 'HMDBSMILES','MBSMILES', 'Annotation_Source', 'Classification_Source']]


# In[175]:


megaA[['id_X', 'premz', 'rtmed', 'int', 'col_eng', 'pol','Formula',
'subclass', 'class', 'superclass', 'Annotation_Source']]


# In[176]:


megaA


# In[177]:


bef_mega = megaA.loc[:,~megaA.columns.duplicated()]
bef_mega


# In[178]:


for i, row in bef_mega.iterrows():
    if not isNaN(bef_mega['Annotation_Source'][i]):
        
        if 'SIRIUS' in bef_mega['Annotation_Source'][i] and 'SIRIUS_Formula' not in bef_mega['Annotation_Source'][i]:
            bef_mega.loc[i, 'SMILES_final'] = bef_mega['SMILES'][i]
            bef_mega['PC_MCSS_SMILES'][i] = np.nan
            bef_mega['KG_MCSS_SMILES'][i] = np.nan
            
            
        elif 'MassBank' in bef_mega['Annotation_Source'][i]:
            bef_mega.loc[i, 'SMILES_final'] = bef_mega['MBSMILES'][i]
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
            bef_mega['Formula'][i] = 'MB_Formula'
            
        elif 'HMDB' in bef_mega['Annotation_Source'][i]:
            bef_mega.loc[i, 'SMILES_final'] = bef_mega['HMDBSMILES'][i]
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

        elif 'KEGG' is bef_mega['Annotation_Source'][i]:
            bef_mega.loc[i, 'SMILES_final'] = bef_mega['KG_SMILES'][i]
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


# In[179]:


bef_mega.columns


# In[180]:


frame = bef_mega[['id_X', 'premz', 'rtmed', 'int', 'col_eng', 'pol', 'Formula', 'MCSS_SMILES', 'most_specific_class', 'level _5', 'subclass', 'class', 'superclass',
       'all_classifications', 'KG_MCSS_SMILES', 'PC_MCSS_SMILES', 'Annotation_Source', 'SMILES_final', 'Classification_Source']]


# In[181]:


np.unique(list(bef_mega['Classification_Source']))


# In[182]:


def classification(frame):
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
            if df_merged["smiles_x"][p] is frame["SMILES_final"][q]:
                frame.loc[q, 'subclass'] = df_merged["subclass.name"][p]
                frame.loc[q, 'class'] = df_merged["class.name"][p]
                frame.loc[q, 'superclass'] = df_merged["superclass.name"][p]
                frame.loc[q, 'Classification_Source'] = "ClassyFire"
    #frame.to_csv(input_dir, '/SIRIUS_combined.csv')
    return(frame)


# In[183]:


clas = classification(frame)
clas


# In[184]:


frame


# In[185]:


frame.to_csv("Candidate_list.csv")


# In[186]:


for i, row in frame.iterrows():
    if not isNaN(frame['SMILES_final'][i]):
        sx = str(pcp.get_compounds(frame['SMILES_final'][i], 'smiles'))
        cp = [int(x) for x in re.findall(r'\b\d+\b', sx)]
        if len(cp) > 0: 
            comp = pcp.Compound.from_cid(cp)
            frame['Formula'][i] = comp.molecular_formula
            frame.loc[i, 'Molecular_mass'] = comp.molecular_weight
            if len(comp.synonyms) > 0:
                frame.loc[i, 'Name'] = comp.synonyms[0]
            else:
                frame.loc[i, 'Name'] = comp.iupac_name


# In[187]:


frame[['id_X', 
       'premz', 
       'rtmed', 
       'int', 
       'col_eng', 
       'pol', 
       'SMILES_final', 
       'Formula', 
       'Name',
       'Molecular_mass', 
       'MCSS_SMILES', 
       'PC_MCSS_SMILES', 
       'KG_MCSS_SMILES', 
       'subclass', 
       'class', 
       'superclass', 
       'Classification_Source', 
       'Annotation_Source'
      ]]


# In[188]:


for i, row in frame.iterrows():
    if not isNaN(frame['Annotation_Source'][i]) and 'SIRIUS_Formula' != frame['Annotation_Source'][i]:
        frame.loc[i, 'Occurence'] = len(frame['Annotation_Source'][i].split())


# In[ ]:





# In[ ]:





# In[ ]:





# In[189]:


frame.to_csv("Candidate_list.csv")


# In[190]:


metadata = ['id_x = index which depends on the filenumber, mz, rt and the id of the feature in that file, the filenumber are added to the new code, but not to these results yet',
           'premz = precursor m/z',
           'rtmed = retention time in median',
           'int = intensity',
           'col_eng = collision energy', 
            'pol = neg or pos mode',
            
           'SMILES_final = SMILES of the candidate from the annotation source',
            'Formula = formula for the candidate',
            'Name = name of the candidate from PubChem',
            'Molecular_mass = mass of the candidate from PubChem',
            
           'MCSS_SMILES = maximum common substructure from the SMILES of the top scoring candidates from SIRIUS',
            'PC_MCSS_SMILES = maximum common substructure from the SMILES of the top scoring candidates from MetFrag (PubChem)',
            'KG_MCSS_SMILES = maximum common substructure from the SMILES of the top scoring candidates from MetFrag (KEGG)',
            
            'most_specific_class',
            'level _5',
            'subclass',
            'class',
            'superclass',
            'all_classifications',
            
            'Classification_Source = either Canopus OR ClassyFire',
            'Annotation_Source = list of the database sources, can be either one or all of the following options: HMDB, GNPS, MassBank, SIRIUS, KEGG(MetFrag), PubChem(MetFrag), SIRIUS_Formula, SuspectList',
           'Occurence = no. of Annotation Sources'
            
           ]


# In[191]:


with open('metadata.txt', 'a') as f:
    f.writelines('\n'.join(metadata))


# ### Comparison with CD results

# In[192]:


results = pd.read_csv("Candidate_list.csv")


# In[193]:


results


# In[194]:


results = (results.drop_duplicates(subset = ['premz', 'rtmed', 'SMILES_final', 'Name']))


# In[195]:


results


# In[196]:


new = results.loc[results['Annotation_Source'] != 'SIRIUS_Formula']


# In[197]:


new.loc[isNaN(new['Annotation_Source'])]


# In[198]:


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


# In[199]:


cd = pd.read_csv('/Users/mahnoorzulfiqar/OneDriveUNI/MZML/CD/checkCD_Results_QC.csv')


# In[200]:


cd


# In[ ]:





# In[ ]:





# In[201]:


cd = cd.drop_duplicates(subset = ['SMILES'])


# In[202]:


for i, row in results.iterrows():
    if not isNaN(results['SMILES_final'][i]):
        if 'invalid_SMILES' not in results['SMILES_final'][i] and 'invalid_chemistry' not in results['SMILES_final'][i]:
            for j, row in cd.iterrows():
                if not isNaN(cd['SMILES'][j]):
                    CGms = [Chem.MolFromSmiles(results['SMILES_final'][i]), Chem.MolFromSmiles(cd['SMILES'][j])]
                    CGfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=1024) for x in CGms]
                    CGtn = DataStructs.FingerprintSimilarity(CGfps[0],CGfps[1])
                    if CGtn == 1 and 'CompoundDiscoverer' not in results['Annotation_Source'][i]:
                        results['Annotation_Source'][i] = results['Annotation_Source'][i] + ', ' + 'CompoundDiscoverer'
                        results['Occurence'][i] = results['Occurence'][i] + 1


# In[203]:


results


# In[204]:


results['Annotation_Source'][1850]


# In[205]:


results.to_csv("Candidate_list.csv")


# In[ ]:




