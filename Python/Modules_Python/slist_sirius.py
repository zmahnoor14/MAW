#!/usr/bin/env python
# coding: utf-8

# In[2]:


#!/usr/bin/env python
# make executable in bash chmod +x PyRun

# Libraries
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem import PandasTools
import sys
import pandas as pd
import numpy as np
import os


def main():
    if len(sys.argv) != 3:
        print("Usage: python slist_sirius.py input_dir, slist_csv, substring = None")
    else:
        slist_sirius(sys.argv[1], sys.argv[2], list1)


def slist_sirius(input_dir, slist_csv, substring=None):
    """slist_sirius is used to create a tsv file that contains a list of 
    SMILES. The function also runs the sirius command custom db to create
    fingerprints for each SMILES in a folder that we by default name as
    SL_Frag/. This fingerprints folder is later used by SIRIUS to use 
    these compounds as a another small list of compounds to match against
    the input spectra fingerprints.
    Since SIRIUS doesn't take disconnected structure, Multiply charged, 
    Incorrect syntax, wild card(*) in smiles; this function removes all
    such SMILES from the Suspect List.

    Parameters:
    input_dir (str): This is the input directory where all the .mzML 
    files and their respective result directories are stored. For this 
    function this directory must contain a csv file that has a column 
    named "SMILES".
    
    slist_csv (str): This is the csv file that contains a column of 
    "SMILES". Additionally this file can contain other information 
    about the compounds, but for this function, column of "SMILES", 
    named as "SMILES" is necessary.
    
    substring (list): provide a list of strings of SMILES that 
    shouldn't be considered, provide a list even if there is one string
    that shouldnt be considered. e.g: "[Fe+2]". 

    Returns:
    tsv: a tsv file of list of SMILES, named as SL_Sirius.tsv, is stored 
    in input_dir
    directory: a directory with compound fragmentations will be created 
    in a folder named SL_Frag/ within the same input_dir
    
    
    Usage:
    slist_sirius("/user/project/", "suspectlist.csv", 
    substring = None)

    """

    def isNaN(string):
        return string != string

    sl = pd.read_csv(slist_csv)

    # define function to neutralize the charged SMILES
    def neutralize_atoms(mol):

        pattern = Chem.MolFromSmarts(
            "[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]"
        )
        at_matches = mol.GetSubstructMatches(pattern)
        at_matches_list = [y[0] for y in at_matches]
        if len(at_matches_list) > 0:
            for at_idx in at_matches_list:
                atom = mol.GetAtomWithIdx(at_idx)
                chg = atom.GetFormalCharge()
                hcount = atom.GetTotalNumHs()
                atom.SetFormalCharge(0)
                atom.SetNumExplicitHs(hcount - chg)
                atom.UpdatePropertyCache()
        return mol

    for i, row in sl.iterrows():
        # remove SMILES with wild card
        if "*" in sl["SMILES"][i]:
            sl = sl.drop(labels=i, axis=0)
    for i, row in sl.iterrows():
        # remove SMILES with any string present in the substring
        if substring:
            if bool([ele for ele in substring if (ele in sl["SMILES"][i])]):
                sl = sl.drop(labels=i, axis=0)
    for i, row in sl.iterrows():
        if "." in sl["SMILES"][i]:
            sl.loc[i, "SMILES"] = sl["SMILES"][i].split(".")[0]
    # Neutralize the charged SMILES
    for i, row in sl.iterrows():
        if "+" in sl["SMILES"][i] or "-" in sl["SMILES"][i]:
            mol = Chem.MolFromSmiles(sl["SMILES"][i])
            neutralize_atoms(mol)
            sl.loc[i, "SMILES"] = Chem.MolToSmiles(mol)

            # Remove multiple charged SMILES
            if "+" in sl["SMILES"][i] or "-" in sl["SMILES"][i]:
                pos = sl["SMILES"][i].count("+")
                neg = sl["SMILES"][i].count("-")
                charge = pos + neg
                if charge > 1:
                    sl = sl.drop(labels=i, axis=0)

    slsirius = pd.DataFrame({"smiles": sl["SMILES"]})
    slsirius.to_csv(input_dir + "SL_Sirius.tsv", sep="\t", header=False, index=False)
    os.system(
        "sirius --input "
        + input_dir
        + "SL_Sirius.tsv custom-db --name=SL_Frag --output "
        + input_dir
    )


sys.argv[1]
sys.argv[2]
list1 = sys.argv[3].split(",")

if __name__ == "__main__":
    main()
