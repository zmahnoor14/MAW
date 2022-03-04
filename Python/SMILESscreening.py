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

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem import PandasTools

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
    def isNaN(string):
        return string != string
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
    
SMILESscreening(sys.argv[1], sys.argv[2], sys.argv[3])

