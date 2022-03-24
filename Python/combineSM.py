#!/usr/bin/env python
# coding: utf-8

# In[ ]:


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

def combineSM(input_dir, metfragcsv, siriuscsv):
    
    """combineSM prioritizes Sirius and Suspect list over PubChem and
    KEGG

    Parameters:
    input_dir (str): This is the input directory where all the .mzML 
    files and their respective result directories are stored.
    sirius (dataframe): result of sirius_curation
    metfrag (dataframe): result of metfrag_curation
    
    Returns:
    dataframe: dataframe of combined curated sirius and metfrag results
    csv: "MetabolomicsResults/combinedSM.csv"
    
    Usage:
    combineSM(input_dir = "usr/project/", metfrag, sirius)

    """
    
    def isNaN(string):
        return string != string
    
    metfrag = pd.read_csv(metfragcsv)
    sirius = pd.read_csv(siriuscsv)
    S_M_CSV = pd.concat([sirius, metfrag], axis = 1, levels = ["id_X"])
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
                        if PStn == 1:
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
                        if SKtn == 1:
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
                    
    S_M_CSV.to_csv(input_dir + "MetabolomicsResults/combinedSM.csv")
    return(S_M_CSV)

combineSM(sys.argv[1], sys.argv[2], sys.argv[3])