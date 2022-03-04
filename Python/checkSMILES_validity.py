#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
    def isNaN(string):
        return string != string
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

checkSMILES_validity(sys.argv[1], sys.argv[2])


# In[ ]:




