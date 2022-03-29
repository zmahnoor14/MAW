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

import pubchempy as pcp
import numpy as np

<<<<<<< HEAD
import pandas as pd

import sys

=======
>>>>>>> 25c6491 (cleaned directory)


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

<<<<<<< HEAD
    # create a new directory to store all results /MetabolomicsResults/
    path = os.path.join(input_dir, "MetabolomicsResults")
    if not os.path.isdir(path):
        os.mkdir(path)
=======
    
>>>>>>> 25c6491 (cleaned directory)
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
            
<<<<<<< HEAD
    
=======
>>>>>>> 25c6491 (cleaned directory)
    combined_csv.to_csv(input_dir + 'MetabolomicsResults/SD_post_processed_combined_results.csv')
    return(combined_csv)

combine_allspec(sys.argv[1])

