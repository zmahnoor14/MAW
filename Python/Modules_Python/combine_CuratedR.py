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


def main():
    if len(sys.argv) != 4:
        print(
            "Usage python3 combine_CuratedR.py input_directory combinedSDBs combinedSMs data_type(optional)"
        )
    else:
        combine_CuratedR(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])


def combine_CuratedR(input_dir, combinedSDBs, combinedSMs, data_type="standards"):

    """combine_CuratedR prioritizes in the following manner: gnps>
    mbank>suspectlist>sirius>hmdb>metfrag

    Parameters:
    input_dir (str): This is the input directory where all the .mzML
    files and their respective result directories are stored.

    curatedSDB: df from specDB_Curation
    combinedSM: df from combineSM

    Returns:
    dataframe: with curated Spectral DB results and CDB (S+M) results
    csv: "MetabolomicsResults/final_curation_without_classes.csv"

    Usage:
    combine_CuratedR(input_dir = "usr/project/", curatedSDB, combinedSM)

    """

    def isNaN(string):
        return string != string

    combinedSDB = pd.read_csv(combinedSDBs)
    combinedSM = pd.read_csv(combinedSMs)
    mega = pd.concat([combinedSM, combinedSDB], axis=1, levels=["id_X"])

    for i, row in mega.iterrows():

        # if only compound database results
        if isNaN(mega["Annotation"][i]) and not isNaN(mega["Annotation_C"][i]):
            mega.loc[i, "Annotation_Source"] = mega["Annotation_C"][i]

        # if only spectral db results
        if not isNaN(mega["Annotation"][i]) and isNaN(mega["Annotation_C"][i]):
            mega.loc[i, "Annotation_Source"] = mega["Annotation"][i]

        # if both have results
        if not isNaN(mega["Annotation"][i]) and not isNaN(mega["Annotation_C"][i]):
            ########THREE OR FOUR SDB SOURCES########

            # if three sdb sources or more
            # prioritize Spectral DBs
            if (
                len(mega["Annotation"][i].split()) >= 3
                and "SIRIUS" in mega["Annotation_C"][i]
            ):
                if "MassBank" in mega["Annotation"][i]:
                    SKms = [
                        Chem.MolFromSmiles(mega["MBSMILES"][i]),
                        Chem.MolFromSmiles(mega["SMILES"][i]),
                    ]
                    SKfps = [
                        AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=2048)
                        for x in SKms
                    ]
                    SKtn = DataStructs.FingerprintSimilarity(SKfps[0], SKfps[1])
                    if SKtn == 1.0:
                        print(SKtn)
                        mega.loc[i, "Annotation_Source"] = (
                            mega["Annotation"][i] + ", " + mega["Annotation_C"][i]
                        )
                    else:
                        mega.loc[i, "Annotation_Source"] = mega["Annotation"][i]

                elif "HMDB" in mega["Annotation"][i]:
                    SKms = [
                        Chem.MolFromSmiles(mega["HMDBSMILES"][i]),
                        Chem.MolFromSmiles(mega["SMILES"][i]),
                    ]
                    SKfps = [
                        AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=2048)
                        for x in SKms
                    ]
                    SKtn = DataStructs.FingerprintSimilarity(SKfps[0], SKfps[1])
                    if SKtn == 1.0:
                        mega.loc[i, "Annotation_Source"] = (
                            mega["Annotation"][i] + ", " + mega["Annotation_C"][i]
                        )
                    else:
                        mega.loc[i, "Annotation_Source"] = mega["Annotation"][i]
            elif (
                len(mega["Annotation"][i].split()) >= 3
                and "SIRIUS" not in mega["Annotation_C"][i]
            ):
                if "KEGG" in mega["Annotation_C"][i]:
                    if "MassBank" in mega["Annotation"][i]:
                        SKms = [
                            Chem.MolFromSmiles(mega["MBSMILES"][i]),
                            Chem.MolFromSmiles(mega["KG_SMILES"][i]),
                        ]
                        SKfps = [
                            AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=2048)
                            for x in SKms
                        ]
                        SKtn = DataStructs.FingerprintSimilarity(SKfps[0], SKfps[1])
                        if SKtn == 1.0:
                            mega.loc[i, "Annotation_Source"] = (
                                mega["Annotation"][i] + ", " + mega["Annotation_C"][i]
                            )
                        else:
                            mega.loc[i, "Annotation_Source"] = mega["Annotation"][i]

                    elif "HMDB" in mega["Annotation"][i]:
                        SKms = [
                            Chem.MolFromSmiles(mega["HMDBSMILES"][i]),
                            Chem.MolFromSmiles(mega["KG_SMILES"][i]),
                        ]
                        SKfps = [
                            AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=2048)
                            for x in SKms
                        ]
                        SKtn = DataStructs.FingerprintSimilarity(SKfps[0], SKfps[1])
                        if SKtn == 1.0:
                            mega.loc[i, "Annotation_Source"] = (
                                mega["Annotation"][i] + ", " + mega["Annotation_C"][i]
                            )
                        else:
                            mega.loc[i, "Annotation_Source"] = mega["Annotation"][i]
                else:
                    mega.loc[i, "Annotation_Source"] = mega["Annotation"][i]

            #######TWO OR ONE SDB SOURCE#########

            # if both 2 SDBs and results from insilico tools
            elif len(mega["Annotation"][i].split()) <= 2:
                mega.loc[i, "Annotation_Source"] = mega["Annotation_C"][i]

        # if no results from any databases
        if (
            isNaN(mega["Annotation"][i])
            and isNaN(mega["Annotation_C"][i])
            and not isNaN(mega["Formula"][i])
        ):
            mega.loc[i, "Annotation_Source"] = "SIRIUS_Formula"

    bef_mega = mega.loc[:, ~mega.columns.duplicated()]
    for i, row in bef_mega.iterrows():
        if not isNaN(bef_mega["Annotation_Source"][i]):
            # check if SIRIUS is in the annotation source but keep in mind it shouldnt be SIRIUS_Formula
            if (
                "SIRIUS" in bef_mega["Annotation_Source"][i]
                and "SIRIUS_Formula" not in bef_mega["Annotation_Source"][i]
            ):
                bef_mega.loc[i, "SMILES_final"] = bef_mega["SMILES"][i]
                bef_mega.loc[i, "CompoundNames"] = bef_mega["name"][i]
                bef_mega["PC_MCSS_SMILES"][i] = np.nan
                bef_mega["KG_MCSS_SMILES"][i] = np.nan
            elif "KEGG" in bef_mega["Annotation_Source"][i]:
                bef_mega.loc[i, "SMILES_final"] = bef_mega["KG_SMILES"][i]
                bef_mega.loc[i, "CompoundNames"] = bef_mega["KG_Name"][i]
                # bef_mega['most_specific_class'][i] = np.nan
                # bef_mega['level _5'][i] = np.nan
                bef_mega["subclass"][i] = np.nan
                bef_mega["class"][i] = np.nan
                bef_mega["superclass"][i] = np.nan
                # bef_mega['all_classifications'][i] = np.nan
                bef_mega["Classification_Source"][i] = np.nan
                bef_mega["MCSS_SMILES"][i] = np.nan
                bef_mega["PC_MCSS_SMILES"][i] = np.nan
                bef_mega["Formula"][i] = np.nan

            elif "GNPS, SuspectList" in bef_mega["Annotation_Source"][i]:
                bef_mega.loc[i, "SMILES_final"] = bef_mega["GLsmiles"][i]
                bef_mega.loc[i, "CompoundNames"] = bef_mega["GLname"][i]
                bef_mega.loc[i, "CompoundNames"]
                # bef_mega['most_specific_class'][i] = np.nan
                # bef_mega['level _5'][i] = np.nan
                bef_mega["subclass"][i] = np.nan
                bef_mega["class"][i] = np.nan
                bef_mega["superclass"][i] = np.nan
                # bef_mega['all_classifications'][i] = np.nan
                bef_mega["MCSS_SMILES"][i] = np.nan
                bef_mega["Classification_Source"][i] = np.nan
                bef_mega["PC_MCSS_SMILES"][i] = np.nan
                bef_mega["KG_MCSS_SMILES"][i] = np.nan
                bef_mega["Formula"][i] = np.nan

            elif "GNPS" in bef_mega["Annotation_Source"][i]:
                bef_mega.loc[i, "SMILES_final"] = bef_mega["GNPSSMILES"][i]
                bef_mega.loc[i, "CompoundNames"] = bef_mega["GNPScompound_name"][i]
                # bef_mega['most_specific_class'][i] = np.nan
                # bef_mega['level _5'][i] = np.nan
                bef_mega["subclass"][i] = np.nan
                bef_mega["class"][i] = np.nan
                bef_mega["superclass"][i] = np.nan
                # bef_mega['all_classifications'][i] = np.nan
                bef_mega["MCSS_SMILES"][i] = np.nan
                bef_mega["Classification_Source"][i] = np.nan
                bef_mega["PC_MCSS_SMILES"][i] = np.nan
                bef_mega["KG_MCSS_SMILES"][i] = np.nan
                bef_mega["Formula"][i] = np.nan
            elif "MassBank" in bef_mega["Annotation_Source"][i]:
                bef_mega.loc[i, "SMILES_final"] = bef_mega["MBSMILES"][i]
                bef_mega.loc[i, "CompoundNames"] = bef_mega["MBcompound_name"][i]
                # bef_mega['most_specific_class'][i] = np.nan
                # bef_mega['level _5'][i] = np.nan
                bef_mega["subclass"][i] = np.nan
                bef_mega["class"][i] = np.nan
                bef_mega["superclass"][i] = np.nan
                # bef_mega['all_classifications'][i] = np.nan
                bef_mega["Classification_Source"][i] = np.nan
                bef_mega["MCSS_SMILES"][i] = np.nan
                bef_mega["PC_MCSS_SMILES"][i] = np.nan
                bef_mega["KG_MCSS_SMILES"][i] = np.nan
                bef_mega["Formula"][i] = np.nan

            elif "PubChem" in bef_mega["Annotation_Source"][i]:
                bef_mega.loc[i, "SMILES_final"] = bef_mega["PC_SMILES"][i]
                bef_mega.loc[i, "CompoundNames"] = bef_mega["PC_Name"][i]
                # bef_mega['most_specific_class'][i] = np.nan
                # bef_mega['level _5'][i] = np.nan
                bef_mega["subclass"][i] = np.nan
                bef_mega["class"][i] = np.nan
                bef_mega["superclass"][i] = np.nan
                # bef_mega['all_classifications'][i] = np.nan
                bef_mega["Classification_Source"][i] = np.nan
                bef_mega["MCSS_SMILES"][i] = np.nan
                bef_mega["KG_MCSS_SMILES"][i] = np.nan
                bef_mega["Formula"][i] = np.nan

            elif "HMDB" in bef_mega["Annotation_Source"][i]:
                bef_mega.loc[i, "SMILES_final"] = bef_mega["HMDBSMILES"][i]
                bef_mega.loc[i, "CompoundNames"] = bef_mega["HMDBcompound_name"][i]
                # bef_mega['most_specific_class'][i] = np.nan
                # bef_mega['level _5'][i] = np.nan
                bef_mega["subclass"][i] = np.nan
                bef_mega["class"][i] = np.nan
                bef_mega["superclass"][i] = np.nan
                # bef_mega['all_classifications'][i] = np.nan
                bef_mega["Classification_Source"][i] = np.nan
                bef_mega["MCSS_SMILES"][i] = np.nan
                bef_mega["PC_MCSS_SMILES"][i] = np.nan
                bef_mega["KG_MCSS_SMILES"][i] = np.nan
                bef_mega["Formula"][i] = np.nan

            elif "SIRIUS_Formula" in bef_mega["Annotation_Source"][i]:
                bef_mega["PC_MCSS_SMILES"][i] = np.nan
                bef_mega["KG_MCSS_SMILES"][i] = np.nan

    bef_megaA = bef_mega[
        [
            "id_X",
            "premz",
            "rtmed",
            "rtmean",
            "int",
            "col_eng",
            "pol",
            "SMILES_final",
            "CompoundNames",
            "MCSS_SMILES",
            "PC_MCSS_SMILES",
            "KG_MCSS_SMILES",
            "subclass",
            "class",
            "superclass",
            "Classification_Source",
            "Annotation_Source",
        ]
    ]

    bef_megaA.rename(columns={"SMILES_final": "SMILES"}, inplace=True)

    Standards = ["Experimental"]
    SpectralDB = ["GNPS", "HMDB", "MassBank"]
    CompoundDB = ["SuspectList", "SIRIUS", "KEGG", "PubChem"]
    Formula = ["SIRIUS_Formula"]

    # bef_megaA['MSI_Level'] = np.nan
    for i, rows in bef_megaA.iterrows():

        if not isNaN(bef_megaA["Annotation_Source"][i]):

            if data_type == "standards":
                bef_megaA.loc[i, "Annotation_Source"] = (
                    bef_megaA["Annotation_Source"][i] + ", Experimental"
                )

                if any(x in bef_megaA["Annotation_Source"][i] for x in SpectralDB):
                    bef_megaA.loc[i, "MSI_Level"] = "Level_1"

                elif any(
                    x in bef_megaA["Annotation_Source"][i] for x in CompoundDB
                ) and not any(x in bef_megaA["Annotation_Source"][i] for x in Formula):
                    bef_megaA.loc[i, "MSI_Level"] = "Level_2/Level_3"

                elif any(x in bef_megaA["Annotation_Source"][i] for x in Formula):
                    bef_megaA.loc[i, "MSI_Level"] = "Level_4"

            else:

                if any(x in bef_megaA["Annotation_Source"][i] for x in SpectralDB):
                    bef_megaA.loc[i, "MSI_Level"] = "Level_2"

                elif any(
                    x in bef_megaA["Annotation_Source"][i] for x in CompoundDB
                ) and not any(x in bef_megaA["Annotation_Source"][i] for x in Formula):
                    bef_megaA.loc[i, "MSI_Level"] = "Level_3"

                elif any(x in bef_megaA["Annotation_Source"][i] for x in Formula):
                    bef_megaA.loc[i, "MSI_Level"] = "Level_4"

        else:
            bef_megaA.loc[i, "MSI_Level"] = "Level_5"

    bef_megaA.to_csv(
        input_dir + "MetabolomicsResults/final_curation_without_classes.csv"
    )
    return bef_megaA


if __name__ == "__main__":
    main()
