#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python
# coding: utf-8

#!/usr/bin/env python
# make executable in bash chmod +x PyRun

# Libraries
import os
import glob
import re

import sys

import csv
import time
import json

import pandas as pd
import pubchempy as pcp
import numpy as np

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem import PandasTools


def main():
    if len(sys.argv) != 4:
        print(
            "Usage python3 SMILESscreening.py input_directory resultcsv complist listname"
        )
    else:
        SMILESscreening(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])


def SMILESscreening(input_dir, resultcsv, complist, listname):
    def isNaN(string):
        return string != string

    """SMILESscreening takes a list of SMILES

    Parameters:
    input_dir (str): This is the input directory where all the .mzML 
    files and their respective result directories are stored.
    
    resultcsv: df from combine_CuratedR or checkSMILES_validity or classification
    complist: list of /n separated txt file conyaining smiles on each line
    listname: name of the list of compounds
    
    Returns:
    dataframe: comparison with another list of compounds
    csv: "MetabolomicsResults/final_curation_with_validSMILES.csv"
    
    Usage:
    checkSMILES_validity(input_dir = "usr/project/", results)

    """

    results = pd.read_csv(resultcsv)
    with open(complist, "r") as text_file:
        cd = text_file.read().split("\n")

    for i, row in results.iterrows():
        if not isNaN(results["SMILES"][i]):
            if (
                "invalid_SMILES" not in results["SMILES"][i]
                and "invalid_chemistry" not in results["SMILES"][i]
            ):
                for j in cd:
                    if not isNaN(j):
                        CGms = [
                            Chem.MolFromSmiles(results["SMILES"][i]),
                            Chem.MolFromSmiles(j),
                        ]
                        CGfps = [
                            AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=1024)
                            for x in CGms
                        ]
                        CGtn = DataStructs.FingerprintSimilarity(CGfps[0], CGfps[1])
                        if (
                            CGtn == 1
                            and listname not in results["Annotation_Source"][i]
                        ):
                            results["Annotation_Source"][i] = (
                                results["Annotation_Source"][i] + ", " + listname
                            )

    frame.to_csv(
        input_dir + "MetabolomicsResults/final_curationListVS" + listname + ".csv"
    )
    return frame


if __name__ == "__main__":
    main()
