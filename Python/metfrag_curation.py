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

def metfrag_curation(input_dir, metfragcsv, sl = True):
    
    
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
    
    def isNaN(string):
        return string != string
    
    metfrag = pd.read_csv(metfragcsv)
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
                else:
                    if sl:
                        
                        #if there is an entry from suspect list WITH kegg
                        if not isNaN(metfrag['KG_SL_comp'][i]) and isNaN(metfrag['PC_SL_comp'][i]):
                            metfrag.loc[i, 'Annotation_M'] = 'KEGG, SuspectList'
                            #metfrag.loc[i, 'SMILES_final'] = metfrag['KG_SMILES'][i]
                
                        #if there is an entry from suspect list WITH kegg
                        elif not isNaN(metfrag['PC_SL_comp'][i]) and isNaN(metfrag['KG_SL_comp'][i]):
                            metfrag.loc[i, 'Annotation_M'] = 'PubChem, SuspectList'
                            #metfrag.loc[i, 'SMILES_final'] = metfrag['PC_SMILES'][i]
                    
                                
    metfrag.to_csv(input_dir + "MetabolomicsResults/metfrag_curated.csv")  
    return(metfrag)
    
metfrag_curation(sys.argv[1], sys.argv[2], sys.argv[3])

