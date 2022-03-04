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
import pubchempy as pcp
import numpy as np

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem import PandasTools
import pubchempy as pcp


def sirius_curation(input_dir, siriuscsv, sl = True):
    
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
    
    def isNaN(string):
        return string != string
    
    sirius = pd.read_csv(siriuscsv)
    for i, row in sirius.iterrows():
    
        # If the explained intensity is greater than 0.70 and there is no suspect list entry
        if sirius['exp_int'][i] >= 0.70 and isNaN(sirius['SL_comp'][i]):
            sirius.loc[i, 'Annotation_S'] = 'SIRIUS'
            #sirius.loc[i, 'SMILES_final'] = sirius['SMILES'][i]
        else:
            if sl:
                
                #If the explained intensity is greater than 0.70 and there is an entry from suspect list
                if sirius['exp_int'][i] >= 0.70 and not isNaN(sirius['SL_comp'][i]):
                    sirius.loc[i, 'Annotation_S'] = 'SIRIUS, SuspectList'
                    #sirius.loc[i, 'SMILES_final'] = sirius['SMILES'][i]
    
                # if the intensity is less thna 0.70 but it still is similar to an entry in Suspect list,
                elif sirius['exp_int'][i] < 0.70 and not isNaN(sirius['SL_comp'][i]):
                    sirius.loc[i, 'Annotation_S'] = 'SIRIUS, SuspectList'
                    #sirius.loc[i, 'SMILES_final'] = sirius['SMILES'][i]
        
    sirius.to_csv(input_dir + "MetabolomicsResults/sirius_curated.csv")
    return(sirius)

sirius_curation(sys.argv[1], sys.argv[2], sys.argv[3])

