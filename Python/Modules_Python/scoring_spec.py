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
import sys
import pubchempy as pcp
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem import PandasTools
import pubchempy as pcp


def main():
    if len(sys.argv) != 2:
        print("Usage python3 scoring_spec.py input_directory spec_file")
    else:
        scoring_spec(sys.argv[1], sys.argv[2])


def scoring_spec(input_dir, spec_file):

    """scoring_spec extracts the candidates with high scores from
    the results from combine_allspec function 

    Parameters:
    input_dir (str): This is the input directory where all the .mzML 
    files and their respective result directories are stored.
    combined (dataframe): dataframe from combine_allspec
    
    Returns:
    dataframe: of the all features and their results
    csv: CSV reuslt file named MetabolomicsResults/combinedSpecDB.csv
    which contains all the features and their Spec DB annotations
    
    Usage:
    scoring_spec(input_dir = "usr/project/", combined)

    """

    def isNaN(string):
        return string != string

    # the scoring highly depends on the following information:
    # similarity scores should be higher than 0.75
    # intScore >=0.50
    # mzScore >= 0.50
    # ratio of the matchingpeaks by the totalpeaks in the query >= 0.50

    combined = pd.read_csv(spec_file)

    def HMDB_Scoring(db, i):
        if (
            db["HMDBmax_similarity"][i] >= 0.75
            and db["HMDBintScore"][i] >= 0.50
            and db["HMDBmzScore"][i] >= 0.50
            and db["HQMatchingPeaks"][i] / db["hQueryTotalPeaks"][i] >= 0.50
        ):
            return True
        else:
            return False

    def GNPS_Scoring(db, i):
        if (
            db["GNPSmax_similarity"][i] >= 0.90
            and db["GNPSintScore"][i] >= 0.50
            and db["GNPSmzScore"][i] >= 0.50
            and db["GQMatchingPeaks"][i] / db["gQueryTotalPeaks"][i] >= 0.50
        ):
            return True
        else:
            return False

    def MB_Scoring(db, i):
        if (
            db["MBmax_similarity"][i] >= 0.50
            and db["MBintScore"][i] >= 0.50
            and db["MBmzScore"][i] >= 0.50
            and db["MQMatchingPeaks"][i] / db["mQueryTotalPeaks"][i] >= 0.50
        ):
            return True
        else:
            return False

    for i, row in combined.iterrows():

        if (
            "HMDBSMILES" in combined.columns
            and "MBSMILES" in combined.columns
            and "GNPSSMILES" in combined.columns
        ):

            # if all DBs show good candidates accorindg to the scoring
            if (
                HMDB_Scoring(combined, i)
                and GNPS_Scoring(combined, i)
                and MB_Scoring(combined, i)
                and not isNaN(combined["GNPSSMILES"][i])
                and not isNaN(combined["MBSMILES"][i])
                and not isNaN(combined["HMDBSMILES"][i])
            ):

                # calulate the tanimoto similarity between the candidates from three DBs

                # hmdb and gnps
                HGms = [
                    Chem.MolFromSmiles(combined["HMDBSMILES"][i]),
                    Chem.MolFromSmiles(combined["GNPSSMILES"][i]),
                ]
                HGfps = [
                    AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=2048)
                    for x in HGms
                ]
                HGtn = DataStructs.FingerprintSimilarity(HGfps[0], HGfps[1])

                # gnps and mbank
                GMms = [
                    Chem.MolFromSmiles(combined["GNPSSMILES"][i]),
                    Chem.MolFromSmiles(combined["MBSMILES"][i]),
                ]
                GMfps = [
                    AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=2048)
                    for x in GMms
                ]
                GMtn = DataStructs.FingerprintSimilarity(GMfps[0], GMfps[1])

                # mbank and hmdb
                HMms = [
                    Chem.MolFromSmiles(combined["HMDBSMILES"][i]),
                    Chem.MolFromSmiles(combined["MBSMILES"][i]),
                ]
                HMfps = [
                    AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=2048)
                    for x in HMms
                ]
                HMtn = DataStructs.FingerprintSimilarity(HMfps[0], HMfps[1])

                # add the following columns
                combined.loc[i, "annotation"] = "HMDB, GNPS, MassBank"
                combined.loc[i, "tanimotoHG"] = HGtn
                combined.loc[i, "tanimotoGM"] = GMtn
                combined.loc[i, "tanimotoHM"] = HMtn
                combined.loc[i, "occurence"] = 3

            # if HMDB and GNPS show good candidates accorindg to the scoring
            if (
                HMDB_Scoring(combined, i)
                and GNPS_Scoring(combined, i)
                and not MB_Scoring(combined, i)
                and not isNaN(combined["GNPSSMILES"][i])
                and not isNaN(combined["HMDBSMILES"][i])
            ):
                HGms = [
                    Chem.MolFromSmiles(combined["HMDBSMILES"][i]),
                    Chem.MolFromSmiles(combined["GNPSSMILES"][i]),
                ]
                HGfps = [
                    AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=2048)
                    for x in HGms
                ]
                HGtn = DataStructs.FingerprintSimilarity(HGfps[0], HGfps[1])

                combined.loc[i, "annotation"] = "HMDB, GNPS"
                combined.loc[i, "tanimotoHG"] = HGtn
                combined.loc[i, "tanimotoGM"] = np.nan
                combined.loc[i, "tanimotoHM"] = np.nan
                combined.loc[i, "occurence"] = 2

            # if MassBank and GNPS show good candidates accorindg to the scoring
            if (
                not HMDB_Scoring(combined, i)
                and GNPS_Scoring(combined, i)
                and MB_Scoring(combined, i)
                and not isNaN(combined["MBSMILES"][i])
                and not isNaN(combined["GNPSSMILES"][i])
            ):
                GMms = [
                    Chem.MolFromSmiles(combined["GNPSSMILES"][i]),
                    Chem.MolFromSmiles(combined["MBSMILES"][i]),
                ]
                GMfps = [
                    AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=2048)
                    for x in GMms
                ]
                GMtn = DataStructs.FingerprintSimilarity(GMfps[0], GMfps[1])

                combined.loc[i, "annotation"] = "GNPS, MassBank"
                combined.loc[i, "tanimotoHG"] = np.nan
                combined.loc[i, "tanimotoGM"] = GMtn
                combined.loc[i, "tanimotoHM"] = np.nan
                combined.loc[i, "occurence"] = 2

            # if MassBank and HMDB show good candidates accorindg to the scoring
            if (
                HMDB_Scoring(combined, i)
                and not GNPS_Scoring(combined, i)
                and MB_Scoring(combined, i)
                and not isNaN(combined["MBSMILES"][i])
                and not isNaN(combined["HMDBSMILES"][i])
            ):
                HMms = [
                    Chem.MolFromSmiles(combined["HMDBSMILES"][i]),
                    Chem.MolFromSmiles(combined["MBSMILES"][i]),
                ]
                HMfps = [
                    AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=2048)
                    for x in HMms
                ]
                HMtn = DataStructs.FingerprintSimilarity(HMfps[0], HMfps[1])

                combined.loc[i, "annotation"] = "HMDB, MassBank"
                combined.loc[i, "tanimotoHG"] = np.nan
                combined.loc[i, "tanimotoGM"] = np.nan
                combined.loc[i, "tanimotoHM"] = HMtn
                combined.loc[i, "occurence"] = 2

            # only HMDB
            if (
                HMDB_Scoring(combined, i)
                and not GNPS_Scoring(combined, i)
                and not MB_Scoring(combined, i)
            ):

                combined.loc[i, "annotation"] = "HMDB"
                combined.loc[i, "tanimotoHG"] = np.nan
                combined.loc[i, "tanimotoGM"] = np.nan
                combined.loc[i, "tanimotoHM"] = np.nan
                combined.loc[i, "occurence"] = 1

            # only GNPS
            if (
                not HMDB_Scoring(combined, i)
                and GNPS_Scoring(combined, i)
                and not MB_Scoring(combined, i)
            ):

                combined.loc[i, "annotation"] = "GNPS"
                combined.loc[i, "tanimotoHG"] = np.nan
                combined.loc[i, "tanimotoGM"] = np.nan
                combined.loc[i, "tanimotoHM"] = np.nan
                combined.loc[i, "occurence"] = 1

            # only MassBank
            if (
                not HMDB_Scoring(combined, i)
                and not GNPS_Scoring(combined, i)
                and MB_Scoring(combined, i)
            ):

                combined.loc[i, "annotation"] = "MassBank"
                combined.loc[i, "tanimotoHG"] = np.nan
                combined.loc[i, "tanimotoGM"] = np.nan
                combined.loc[i, "tanimotoHM"] = np.nan
                combined.loc[i, "occurence"] = 1

            # none
            if (
                not HMDB_Scoring(combined, i)
                and not GNPS_Scoring(combined, i)
                and not MB_Scoring(combined, i)
            ):
                combined.loc[i, "annotation"] = "none"
                combined.loc[i, "tanimotoHG"] = np.nan
                combined.loc[i, "tanimotoGM"] = np.nan
                combined.loc[i, "tanimotoHM"] = np.nan
                combined.loc[i, "occurence"] = 0

        if (
            "HMDBSMILES" not in combined.columns
            and "MBSMILES" in combined.columns
            and "GNPSSMILES" in combined.columns
        ):

            # if MassBank and GNPS show good candidates accorindg to the scoring
            if (
                GNPS_Scoring(combined, i)
                and MB_Scoring(combined, i)
                and not isNaN(combined["MBSMILES"][i])
                and not isNaN(combined["GNPSSMILES"][i])
            ):
                GMms = [
                    Chem.MolFromSmiles(combined["GNPSSMILES"][i]),
                    Chem.MolFromSmiles(combined["MBSMILES"][i]),
                ]
                GMfps = [
                    AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=2048)
                    for x in GMms
                ]
                GMtn = DataStructs.FingerprintSimilarity(GMfps[0], GMfps[1])

                combined.loc[i, "annotation"] = "GNPS, MassBank"
                combined.loc[i, "tanimotoGM"] = GMtn
                combined.loc[i, "occurence"] = 2
            # only GNPS
            if GNPS_Scoring(combined, i) and not MB_Scoring(combined, i):

                combined.loc[i, "annotation"] = "GNPS"
                combined.loc[i, "tanimotoGM"] = np.nan
                combined.loc[i, "occurence"] = 1

            # only MassBank
            if not GNPS_Scoring(combined, i) and MB_Scoring(combined, i):

                combined.loc[i, "annotation"] = "MassBank"
                combined.loc[i, "tanimotoGM"] = np.nan
                combined.loc[i, "occurence"] = 1

            # none
            if not GNPS_Scoring(combined, i) and not MB_Scoring(combined, i):
                combined.loc[i, "annotation"] = "none"
                combined.loc[i, "tanimotoGM"] = np.nan
                combined.loc[i, "occurence"] = 0

        if (
            "HMDBSMILES" in combined.columns
            and "MBSMILES" not in combined.columns
            and "GNPSSMILES" in combined.columns
        ):
            # if HMDB and GNPS show good candidates accorindg to the scoring
            if (
                HMDB_Scoring(combined, i)
                and GNPS_Scoring(combined, i)
                and not isNaN(combined["GNPSSMILES"][i])
                and not isNaN(combined["HMDBSMILES"][i])
            ):
                HGms = [
                    Chem.MolFromSmiles(combined["HMDBSMILES"][i]),
                    Chem.MolFromSmiles(combined["GNPSSMILES"][i]),
                ]
                HGfps = [
                    AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=2048)
                    for x in HGms
                ]
                HGtn = DataStructs.FingerprintSimilarity(HGfps[0], HGfps[1])

                combined.loc[i, "annotation"] = "HMDB, GNPS"
                combined.loc[i, "tanimotoHG"] = HGtn
                combined.loc[i, "occurence"] = 2

            # only HMDB
            if HMDB_Scoring(combined, i) and not GNPS_Scoring(combined, i):

                combined.loc[i, "annotation"] = "HMDB"
                combined.loc[i, "tanimotoHG"] = np.nan
                combined.loc[i, "occurence"] = 1

            # only GNPS
            if not HMDB_Scoring(combined, i) and GNPS_Scoring(combined, i):

                combined.loc[i, "annotation"] = "GNPS"
                combined.loc[i, "tanimotoHG"] = np.nan
                combined.loc[i, "occurence"] = 1
            # none
            if not HMDB_Scoring(combined, i) and not GNPS_Scoring(combined, i):
                combined.loc[i, "annotation"] = "none"
                combined.loc[i, "tanimotoHG"] = np.nan
                combined.loc[i, "occurence"] = 0

        if (
            "HMDBSMILES" in combined.columns
            and "MBSMILES" in combined.columns
            and "GNPSSMILES" not in combined.columns
        ):

            # if MassBank and HMDB show good candidates accorindg to the scoring
            if (
                HMDB_Scoring(combined, i)
                and MB_Scoring(combined, i)
                and not isNaN(combined["MBSMILES"][i])
                and not isNaN(combined["HMDBSMILES"][i])
            ):
                HMms = [
                    Chem.MolFromSmiles(combined["HMDBSMILES"][i]),
                    Chem.MolFromSmiles(combined["MBSMILES"][i]),
                ]
                HMfps = [
                    AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=2048)
                    for x in HMms
                ]
                HMtn = DataStructs.FingerprintSimilarity(HMfps[0], HMfps[1])

                combined.loc[i, "annotation"] = "HMDB, MassBank"
                combined.loc[i, "tanimotoHM"] = HMtn
                combined.loc[i, "occurence"] = 2

            # only HMDB
            if HMDB_Scoring(combined, i) and not MB_Scoring(combined, i):

                combined.loc[i, "annotation"] = "HMDB"
                combined.loc[i, "tanimotoHM"] = np.nan
                combined.loc[i, "occurence"] = 1

            # only MassBank
            if not HMDB_Scoring(combined, i) and MB_Scoring(combined, i):

                combined.loc[i, "annotation"] = "MassBank"
                combined.loc[i, "tanimotoHM"] = np.nan
                combined.loc[i, "occurence"] = 1

            # none
            if not HMDB_Scoring(combined, i) and not MB_Scoring(combined, i):
                combined.loc[i, "annotation"] = "none"
                combined.loc[i, "tanimotoHM"] = np.nan
                combined.loc[i, "occurence"] = 0

        # If only HMDB was used

        if (
            "HMDBSMILES" in combined.columns
            and "MBSMILES" not in combined.columns
            and "GNPSSMILES" not in combined.columns
        ):
            # only HMDB
            if HMDB_Scoring(combined, i):

                combined.loc[i, "annotation"] = "HMDB"
                combined.loc[i, "occurence"] = 1

            # none
            if not HMDB_Scoring(combined, i):
                combined.loc[i, "annotation"] = "none"
                combined.loc[i, "occurence"] = 0

        # If only MassBank was used

        if (
            "HMDBSMILES" not in combined.columns
            and "MBSMILES" in combined.columns
            and "GNPSSMILES" not in combined.columns
        ):
            # only MassBank
            if MB_Scoring(combined, i):

                combined.loc[i, "annotation"] = "MassBank"
                combined.loc[i, "occurence"] = 1

            # none
            if not MB_Scoring(combined, i):
                combined.loc[i, "annotation"] = "none"
                combined.loc[i, "occurence"] = 0

        # If only GNPS was used

        if (
            "HMDBSMILES" not in combined.columns
            and "MBSMILES" not in combined.columns
            and "GNPSSMILES" in combined.columns
        ):
            # only GNPS
            if GNPS_Scoring(combined, i):

                combined.loc[i, "annotation"] = "GNPS"
                combined.loc[i, "occurence"] = 1

            # none
            if not GNPS_Scoring(combined, i):
                combined.loc[i, "annotation"] = "none"
                combined.loc[i, "occurence"] = 0

    combined.to_csv(input_dir + "MetabolomicsResults/scoredSpecDB.csv")
    return combined


if __name__ == "__main__":
    main()
