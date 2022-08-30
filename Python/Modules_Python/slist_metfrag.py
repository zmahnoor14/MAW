#!/usr/bin/env python
# coding: utf-8

# In[8]:


#!/usr/bin/env python
# make executable in bash chmod +x PyRun

# Libraries
import pandas as pd
import sys
import numpy as np
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem import PandasTools


def main():
    if len(sys.argv) != 3:
        print("Usage: python slist_metfrag.py input_dir, slist_csv, name")
    else:
        slist_metfrag(sys.argv[1], sys.argv[2], sys.argv[3])


# make sure your Smiles entries in the suspect list csv are in a column named "SMILES"
def slist_metfrag(input_dir, slist_csv, name):
    """slist_metfrag is used to create a txt file that contains a list of 
    InChIKeys. This list is later used by MetFrag to use these compounds 
    as a Suspect List.

    Parameters:
    input_dir (str): This is the input directory where all the .mzML 
    files and their respective result directories are stored. For this 
    function this directory must contain a csv file that has a column 
    named "SMILES".
    
    slist_csv (str): This is the csv file that contains a column of 
    "SMILES". Additionally this file can contain other information 
    about the compounds, but for this function, column of "SMILES", 
    named as "SMILES" is necessary.

    Returns:
    list: list of InChIKeys
    txt: a txt file of list of InChIKeys, is stored in input_dir
    
    Usage:
    slist_metfrag(input_dir = "/user/project/", slist_csv = 
    "suspectlist.csv")
    
    """

    def isNaN(string):
        return string != string

    sl = pd.read_csv(slist_csv)
    sl_mtfrag = []
    for i, rows in sl.iterrows():
        if i is not None:
            mols = Chem.MolFromSmiles(sl["SMILES"][i])
            try:
                sl.loc[i, "InChIKey"] = Chem.inchi.MolToInchiKey(mols)
                sl_mtfrag.append(sl["InChIKey"][i])
            except Exception as e:
                print(e)

    with open((input_dir + "/SL_" + name + ".txt"), "w") as filehandle:
        for listitem in sl_mtfrag:
            filehandle.write("%s\n" % listitem)

    return sl_mtfrag


if __name__ == "__main__":
    main()
