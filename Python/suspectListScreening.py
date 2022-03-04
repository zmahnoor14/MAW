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

def suspectListScreening(input_dir, slistcsv, SpectralDB_Results):
    
    """suspectListScreening runs tanoimoto similarity score to between
    compounds from the results from spectral DBs and suspect list

    Parameters:
    input_dir (str): This is the input directory where all the .mzML 
    files and their respective result directories are stored.
    slistcsv (str): path to suspect list
    SpectralDB_Results (dataframe): dataframe from scoring_spec
    
    Returns:
    dataframe: all features and specDB reults and suspect list screening 
    results
    csv: CSV reuslt file named MetabolomicsResults/SpecDBvsSL.csv
    which contains all the features and their Spec DB annotations
    and suspect list occurences if any
    
    Usage:
    suspectListScreening(input_dir = "usr/project/",
    slistcsv = "usr/project/suspect_list.csv", 
    SpectralDB_Results)

    """
    
    SpectralDB_Results = pd.read_csv(SpectralDB_Results)
    
    # add columns to the resulst from scoring_spec
    # these columns are for high similiarity canidtes between the databases and suspect list
    SpectralDB_Results['HLsmiles'] = np.nan
    SpectralDB_Results['HLname'] = np.nan
    SpectralDB_Results['GLsmiles'] = np.nan
    SpectralDB_Results['GLname'] = np.nan
    SpectralDB_Results['MLsmiles'] = np.nan
    SpectralDB_Results['MLname'] = np.nan
    
    Suspect_list = pd.read_csv(slistcsv)
    
    for i, row in SpectralDB_Results.iterrows():
        if not isNaN(SpectralDB_Results['HMDBSMILES'][i]) and SpectralDB_Results['HMDBSMILES'][i] != " ":
            for j, row in Suspect_list.iterrows():
                LHms2 = [Chem.MolFromSmiles(SpectralDB_Results['HMDBSMILES'][i]), Chem.MolFromSmiles(Suspect_list['SMILES'][j])]
                LHfps2 = [AllChem.GetMorganFingerprintAsBitVect(x2,2, nBits=2048) for x2 in LHms2]
                LHtn2 = DataStructs.FingerprintSimilarity(LHfps2[0],LHfps2[1])
                if LHtn2 >= 0.9:
                    SpectralDB_Results.loc[i, 'HLsmiles'] = Suspect_list['SMILES'][j]
                    SpectralDB_Results.loc[i, 'HLname'] = Suspect_list['Name'][j]
                    
                    
        if not isNaN(SpectralDB_Results['GNPSSMILES'][i]) and SpectralDB_Results['GNPSSMILES'][i] != " ":
            for k, row in Suspect_list.iterrows():
                LGms2 = [Chem.MolFromSmiles(SpectralDB_Results['GNPSSMILES'][i]), Chem.MolFromSmiles(Suspect_list['SMILES'][k])]
                LGfps2 = [AllChem.GetMorganFingerprintAsBitVect(x2,2, nBits=2048) for x2 in LGms2]
                LGtn2 = DataStructs.FingerprintSimilarity(LGfps2[0],LGfps2[1])
                if LGtn2 >= 0.9:
                    SpectralDB_Results.loc[i, 'GLsmiles'] = Suspect_list['SMILES'][k]
                    SpectralDB_Results.loc[i, 'GLname'] = Suspect_list['Name'][k]
                    
                    
        if not isNaN(SpectralDB_Results['MBSMILES'][i]) and SpectralDB_Results['MBSMILES'][i] != " ":
            for l, row in Suspect_list.iterrows():
                LMms2 = [Chem.MolFromSmiles(SpectralDB_Results['MBSMILES'][i]), Chem.MolFromSmiles(Suspect_list['SMILES'][l])]
                LMfps2 = [AllChem.GetMorganFingerprintAsBitVect(x2,2, nBits=2048) for x2 in LMms2]
                LMtn2 = DataStructs.FingerprintSimilarity(LMfps2[0],LMfps2[1])
                if LMtn2 >= 0.9:
                    SpectralDB_Results.loc[i, 'MLsmiles'] = Suspect_list['SMILES'][l]
                    SpectralDB_Results.loc[i, 'MLname'] = Suspect_list['Name'][l]
    
    # add annotations and occurences
    for i, row in SpectralDB_Results.iterrows():
        if not isNaN(SpectralDB_Results['HLname'][i]) or not isNaN(SpectralDB_Results['GLname'][i]) or not isNaN(SpectralDB_Results['MLname'][i]):
            SpectralDB_Results['occurence'][i] = SpectralDB_Results['occurence'][i] + 1
            SpectralDB_Results['annotation'][i] = SpectralDB_Results['annotation'][i] + ', Suspect_List'
            
    SpectralDB_Results.to_csv(input_dir + "MetabolomicsResults/SpecDBvsSL.csv")
    return(SpectralDB_Results)

suspectListScreening(sys.argv[1], sys.argv[2], sys.argv[3])

