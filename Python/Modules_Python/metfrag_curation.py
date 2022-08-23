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
import sys
import pandas as pd
import numpy as np

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem import PandasTools
import pubchempy as pcp

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
    
metfrag_curation(sys.argv[1], sys.argv[2], sys.argv[3])

