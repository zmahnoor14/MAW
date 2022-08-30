#!/usr/bin/env python
# coding: utf-8

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

import pubchempy as pcp
import numpy as np


import os
import glob
import re


def main():
    if len(sys.argv) != 2:
        print("Usage python3 sirius_postProc2.py input_directory input_tablecsv")
    else:
        sirius_postProc2(sys.argv[1], sys.argv[2])


def sirius_postProc2(input_dir, input_tablecsv):
    def isNaN(string):
        return string != string

    """sirius_postProc2 is the second part of the function 
    sirius_postProc defined in R part of the workflow. This function
    re-checks the Suspect list, if present or given as a parameter, 
    whether the candidates have a high similarity with compounds in
    Suspect List. It also calculates the Maximum Common Substructure
    (MCSS)

    Parameters:
    input_dir (str): This is the input directory where all the .mzML 
    files and their respective result directories are stored. For this 
    function this directory must contain a csv file that has a column 
    named "SMILES".
    
    input_tablecsv (str): This is the table in csv format (defined in R), 
    which stores a csv table containing columns "mzml_files", which 
    contains liat of all input files with their relative paths, second
    column is "ResultFileName" which is a list of the corresponding
    result relative directories to each mzml files. Lastly, "file_id", 
    contains a file directory. This table will be used to read the 
    SIRIUS json files
    

    Returns:
    csv: a result file with additional columns such as those for suspect
    list if one is used. It also adds columns on MCSS., named as 
    "input_dir/ResultFileName/insilico/SiriusResults.csv"
    
    
    Usage:
    sirius_postProc2(input_dir = "/user/project/", 
    input_table = "/user/project/suspectlist.csv")


    """

    # Describe the heavy atoms to be considered for MCSS
    heavy_atoms = ["C", "N", "P", "O", "S"]

    input_table = pd.read_csv(input_tablecsv)

    for m, row in input_table.iterrows():

        # Read the file result_dir/insilico/MS1DATAsirius.csv.
        # This file has been produced in R workflow and contains
        # SIRIUS results.

        file1 = pd.read_csv(
            input_dir
            + (
                input_table["ResultFileNames"][m] + "/insilico/MS1DATAsirius.csv"
            ).replace("./", "")
        )

        for i, row in file1.iterrows():

            # if the entry has SMILES extracted for MCSS calculation
            if not isNaN(file1["SMILESforMCSS"][i]):

                # split the SMILES using |
                top_smiles = file1["SMILESforMCSS"][i].split("|")

                # if there are more than 1 smiles in the top smiles,
                if len(top_smiles) > 1:
                    mol = []
                    for j in top_smiles:
                        n = Chem.MolFromSmiles(j)
                        mol.append(n)
                    # list of mol used to calaculate the MCSS
                    res = rdFMCS.FindMCS(mol)
                    sm_res = Chem.MolToSmiles(Chem.MolFromSmarts(res.smartsString))
                    # Check if the MCSS has one of the heavy atoms and whether they are
                    # more than 3
                    elem = [ele for ele in heavy_atoms if (ele in sm_res)]
                    if elem and len(sm_res) >= 3:
                        file1.loc[i, "MCSSstring"] = res.smartsString
                        file1.loc[i, "MCSS_SMILES"] = Chem.MolToSmiles(
                            Chem.MolFromSmarts(res.smartsString)
                        )

            if file1["FormulaRank"][i] == 1.0:
                sep = "json/"
                strpd = file1["dir"][i].split(sep, 1)[0] + "json/canopus_summary.tsv"
                if os.path.isfile(strpd):

                    canopus = pd.read_csv(strpd, sep="\t")
                    if len(canopus) > 0:
                        # file1.loc[i, 'most_specific_class'] = canopus["most specific class"][0]
                        # file1.loc[i, 'level _5'] = canopus["level 5"][0]
                        file1.loc[i, "subclass"] = canopus["subclass"][0]
                        file1.loc[i, "class"] = canopus["class"][0]
                        file1.loc[i, "superclass"] = canopus["superclass"][0]
                        # file1.loc[i, 'all_classifications'] = canopus["all classifications"][0]
                        file1.loc[i, "Classification_Source"] = "CANOPUS"

        file1.to_csv(
            input_dir
            + (
                input_table["ResultFileNames"][m] + "/insilico/SiriusResults.csv"
            ).replace("./", "")
        )
        return file1


if __name__ == "__main__":
    main()
