#!/usr/bin/env python
# coding: utf-8

# In[24]:


import pandas as pd
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import pubchempy as pcp
import numpy as np
def isNaN(string):
    return string != string
import os
import glob
import re
from pybatchclassyfire import *
import csv 
import time
import json
from pandas import json_normalize
import pandas.io.formats.style
from rdkit.Chem import PandasTools
import xlrd
import openpyxl
import statistics


# In[29]:


from Workflow_Python_Functions import *


# In[19]:


input_dir = os.getcwd() +"/"


# In[20]:


input_dir


# In[9]:


input_table = pd.read_csv("/Users/mahnoorzulfiqar/Standards_CodeSet/ input_table.csv")
input_table


# In[21]:


hmdb_sdf = '/Users/mahnoorzulfiqar/OneDriveUNI/S_CResults/structures.sdf'


# In[22]:


dframe = PandasTools.LoadSDF(hmdb_sdf,
                                idName='HMDB_ID',
                                smilesName='SMILES',
                                molColName='Molecule',
                                includeFingerprints=False)


# In[ ]:


Suspect_list = pd.read_csv('/Users/mahnoorzulfiqar/OneDriveUNI/SuspectList/Use_This_CURATED_SUSPECT_LIST.csv')


# In[30]:


spec_df = spec_postproc(input_dir, dframe)


# In[31]:


Merged_Result_df = combine_specdb(spec_df)


# In[32]:


combined = combine_allspec(Merged_Result_df)


# In[33]:


SpectralDB_Results = scoring_spec(combined)


# In[34]:


suspectListScreening(Suspect_list, SpectralDB_Results)


# In[ ]:




