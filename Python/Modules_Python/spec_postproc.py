#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python
# coding: utf-8

#!/usr/bin/env python
# make executable in bash chmod +x PyRun

# Libraries

import sys

import pandas as pd
import numpy as np

import pubchempy as pcp

import os
import glob
import re
from pandas import json_normalize
from rdkit.Chem import PandasTools


def main():
    if len(sys.argv) != 2:
        print("Usage: python spec_postproc.py input_dir, Source = all(default)")
    else:
        spec_postproc(sys.argv[1], sys.argv[2])


def spec_postproc(input_dir, Source="all"):

    """spec_postproc function processes the resulst from dereplication 
    using different spectral DBs. 

    Parameters:
    input_dir (str): This is the input directory where all the .mzML 
    files and their respective result directories are stored.
    
    Source (str): either "mbank" or "hmdb" or "gnps", or "all"

    Returns:
    
    dataframe: of the paths of the processed DB results
    
    
    Usage:
    spec_postproc(input_dir = "/user/project/", Source = "all")

    """

    def isNaN(string):
        return string != string

    # empty lists of csv files paths for each database
    GNPScsvfiles = []
    HMDBcsvfiles = []
    MassBankcsvfiles = []

    # list all files and directories
    for entry in os.listdir(input_dir):
        if os.path.isdir(os.path.join(input_dir, entry)):

            # enter the directory with /spectral_dereplication/ results
            sub_dir = input_dir + entry + "/spectral_dereplication"
            if os.path.exists(sub_dir):
                files = glob.glob(sub_dir + "/*.csv")

                for f in files:
                    if "gnps." in f:
                        GNPScsvfiles.append(f)
                    if "hmdb." in f:
                        HMDBcsvfiles.append(f)
                    if "mbank." in f:
                        MassBankcsvfiles.append(f)

    if Source == "hmdb" or Source == "all":

        if not os.path.exists(input_dir + "structures.sdf"):
            # download SDF structures
            os.system(
                "wget -P "
                + input_dir
                + " https://hmdb.ca/system/downloads/current/structures.zip"
            )
            os.system("unzip " + input_dir + "structures.zip" + " -d " + input_dir)

        # Load the sdf
        dframe = PandasTools.LoadSDF(
            (input_dir + "structures.sdf"),
            idName="HMDB_ID",
            smilesName="SMILES",
            molColName="Molecule",
            includeFingerprints=False,
        )

        #### read sdf file from HMDB to collect names and smiles ####

        # HMDB CSV Result file pre_processing

        # open another csv path holding empty list, which will be filled
        # with post processed csv results
        HMDBcsvfiles2 = []

        for k in HMDBcsvfiles:

            # read the csv files
            hmdb_df = pd.read_csv(k)

            # merge on basis of id, frame and hmdb result files
            SmilesHM = pd.merge(
                hmdb_df,
                dframe,
                left_on=hmdb_df.HMDBcompoundID,
                right_on=dframe.DATABASE_ID,
            )

            for i, row in hmdb_df.iterrows():

                for j, row in SmilesHM.iterrows():

                    # where index for both match, add the name and SMILES
                    if hmdb_df["id_X"][i] == SmilesHM["id_X"][j]:
                        hmdb_df.loc[i, "HMDBSMILES"] = SmilesHM["SMILES"][
                            j
                        ]  # add SMILES
                        hmdb_df.loc[i, "HMDBcompound_name"] = SmilesHM["GENERIC_NAME"][
                            j
                        ]  # add name
                        hmdb_df.loc[i, "HMDBformula"] = SmilesHM["FORMULA"][
                            j
                        ]  # add formula

            csvname = (
                (os.path.splitext(k)[0]) + "proc" + ".csv"
            )  # name for writing it in a new file
            hmdb_df.to_csv(csvname)  # write
            HMDBcsvfiles2.append(csvname)  # add to a list
            dict1 = {"HMDBr": HMDBcsvfiles2}
            df = pd.DataFrame(dict1)

    # MassBank CSV Result file pre_processing

    if Source == "mbank" or Source == "all":

        # open another csv path holding empty list, which will be filled
        # with post processed csv results
        MassBankcsvfiles2 = []

        for l in MassBankcsvfiles:

            # read mbank csv file
            mbank_df = pd.read_csv(l)

            for i, row in mbank_df.iterrows():

                inchiK = str(mbank_df["MBinchiKEY"][i])

                # extract inchikeys
                y = pcp.get_compounds(inchiK, "inchikey")  # compound based on inchikey

                for compound in y:

                    # add smiles
                    smles = compound.isomeric_smiles
                    mbank_df.loc[i, "MBSMILES"] = smles

            csvname = (os.path.splitext(l)[0]) + "proc" + ".csv"
            mbank_df.to_csv(csvname)
            MassBankcsvfiles2.append(csvname)

            dict1 = {"MBr": MassBankcsvfiles2}
            df = pd.DataFrame(dict1)

    # GNPS CSV Result file pre_processing
    if Source == "gnps" or Source == "all":
        # open another csv path holding empty list, which will be filled
        # with post processed csv results
        GNPScsvfiles2 = []
        # currently only these subsets are removed from the names from GNPS
        matches = [
            "M+",
            "[M",
            "M-",
            "2M",
            "M*" "20.0",
            "50.0",
            "30.0",
            "40.0",
            "60.0",
            "70.0",
            "eV",
            "Massbank",
            "Spectral",
            "Match",
            "to",
            "from",
            "NIST14",
            "MoNA",
            "[IIN-based:",
            "[IIN-based",
            "on:",
            "CCMSLIB00003136269]",
        ]

        for l in GNPScsvfiles:
            gnps_df = pd.read_csv(l)

            for i, row in gnps_df.iterrows():
                # if compound name is present
                if not isNaN(gnps_df["GNPScompound_name"][i]):
                    # split if there is a gap in the names
                    string_chng = gnps_df["GNPScompound_name"][i].split(" ")

                    # create an empty list
                    newstr = []

                    # for each part of the string in the names
                    chng = []

                    for j in range(len(string_chng)):

                        # check if the substrings are present in the matches and no - is present
                        if not any(
                            x in string_chng[j] for x in matches
                        ):  # and not '-' == string_chng[j]:

                            # IF | and ! not in the substring
                            if "|" not in string_chng[j] or "!" not in string_chng[j]:
                                newstr.append(string_chng[j])

                            # if | present in the substring
                            elif "|" in string_chng[j]:

                                # split the string
                                jlen = string_chng[j].split("|")
                                # how many substrings are left now
                                lst = len(jlen) - 1
                                # append this to chng
                                chng.append(jlen[lst])
                                break

                    # now append chng to newstr
                    chng.append(" ".join(newstr))
                    # save this as the correct name
                    gnps_df.loc[i, "corr_names"] = chng[0]
                    if not isNaN(gnps_df["GNPSSMILES"][i]):
                        if chng == "":
                            break
                        elif gnps_df["GNPSSMILES"][i].isalpha():
                            s = pcp.get_compounds(chng[0], "name")
                            if s:
                                for comp in s:
                                    gnps_df["GNPSSMILES"][i] = comp.isomeric_smiles
                            else:
                                gnps_df["GNPSSMILES"][i] = ""
                else:
                    gnps_df["GNPSSMILES"][i] = ""

            for i, row in gnps_df.iterrows():
                if isNaN(gnps_df["GNPSSMILES"][i]):
                    if "[" in gnps_df["GNPScompound_name"][i].split(" ")[-1]:
                        string_chng = gnps_df["GNPScompound_name"][i].split("[")
                        # print(gnps_df['GNPScompound_name'][i])
                        keep_names = []
                        for j in range(len(string_chng) - 1):
                            gnps_df.loc[i, "corr_names"] == string_chng[j]
                            s = pcp.get_compounds(string_chng[j], "name")

                            if s:
                                for comp in s:
                                    gnps_df["GNPSSMILES"][i] = comp.isomeric_smiles
                            else:
                                gnps_df["GNPSSMILES"][i] = ""
                if not isNaN(gnps_df["GNPSSMILES"][i]):
                    try:
                        sx = pcp.get_compounds(gnps_df["GNPSSMILES"][i], "smiles")
                        if sx:
                            sx = str(sx)
                            comp = pcp.Compound.from_cid(
                                [int(x) for x in re.findall(r"\b\d+\b", sx)]
                            )
                            gnps_df.loc[i, "GNPSformula"] = comp.molecular_formula
                    except:
                        gnps_df.loc[i, "GNPSformula"] = ""

            csvname = (os.path.splitext(l)[0]) + "proc" + ".csv"
            gnps_df.to_csv(csvname)
            GNPScsvfiles2.append(csvname)
            dict1 = {"GNPSr": GNPScsvfiles2}
            df = pd.DataFrame(dict1)

    if Source == "all":

        dict1 = {
            "GNPSr": GNPScsvfiles2,
            "HMDBr": HMDBcsvfiles2,
            "MBr": MassBankcsvfiles2,
        }
        df = pd.DataFrame(dict1)

        return df


if __name__ == "__main__":
    main()
