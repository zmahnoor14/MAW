#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python
# coding: utf-8

#!/usr/bin/env python
#make executable in bash chmod +x PyRun

# Libraries
import os
import glob
import re

import sys

import csv 
import time
import json

import pubchempy as pcp
import numpy as np
<<<<<<< HEAD
import pandas as pd
from pandas import json_normalize

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem import PandasTools

=======
from pandas import json_normalize

>>>>>>> 25c6491 (cleaned directory)
from pybatchclassyfire import *



<<<<<<< HEAD
=======

>>>>>>> 25c6491 (cleaned directory)
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
    
    def isNaN(string):
        return string != string
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
<<<<<<< HEAD
    

    frame.to_csv(input_dir + "MetabolomicsResults/final_curationList.csv")
    return(frame)
=======
    return(frame)

    frame.to_csv(input_dir + "MetabolomicsResults/final_curationList.csv")
    
>>>>>>> 25c6491 (cleaned directory)
classification(sys.argv[1], sys.argv[2])

