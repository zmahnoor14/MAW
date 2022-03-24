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

import csv 
import time
import json

import sys

import pandas as pd

import pubchempy as pcp
import numpy as np


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

combine_specdb(sys.argv[1])

