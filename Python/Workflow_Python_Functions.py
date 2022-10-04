#!/usr/bin/env python
# coding: utf-8

# In[6]:


#!/usr/bin/env python
# coding: utf-8

import glob
import json
import os
import re
import time
import wget
import urllib.parse


# In[7]:


import numpy as np
import pandas as pd
import pubchempy as pcp


# In[8]:


from pybatchclassyfire import *
from pandas import json_normalize
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem import PandasTools


# In[9]:


def isNaN(string):
    return string != string


# In[10]:


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


# In[ ]:





# In[7]:


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

    # Define scoring for all DBs
    def HMDB_Scoring(db, i):
        if (
            db["HMDBintScore"][i] >= 0.50
            and db["HMDBmzScore"][i] >= 0.50
            and db["HQMatchingPeaks"][i] / db["hQueryTotalPeaks"][i] >= 0.50
        ):
            return True
        else:
            return False

    def GNPS_Scoring(db, i):
        if (
            db["GNPSintScore"][i] >= 0.50
            and db["GNPSmzScore"][i] >= 0.50
            and db["GQMatchingPeaks"][i] / db["gQueryTotalPeaks"][i] >= 0.50
        ):
            return True
        else:
            return False

    def MB_Scoring(db, i):
        if (
            db["MBintScore"][i] >= 0.50
            and db["MBmzScore"][i] >= 0.50
            and db["MQMatchingPeaks"][i] / db["mQueryTotalPeaks"][i] >= 0.50
        ):
            return True
        else:
            return False

    # list all files and directories
    for entry in os.listdir(input_dir):

        if os.path.isdir(os.path.join(input_dir, entry)):
            # print(entry)

            msp_file = glob.glob(
                input_dir + "/" + entry + "/spectral_dereplication" + "/*.csv"
            )
            # print(msp_file)

            if len(msp_file) > 0:

                if os.path.exists(msp_file[0]):

                    msp = pd.read_csv(msp_file[0])
                    # enter the directory with /spectral_dereplication/ results

                    # enter the directory with /spectral_dereplication/ results
                    # GNPS Results
                    if Source == "gnps" or Source == "all":
                        msp["gnps_results_csv"] = np.nan

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

                        # open another csv path holding empty list, which will be filled
                        # with post processed csv results
                        # GNPScsvfiles2 = []

                        # print(entry)
                        # enter the directory with /spectral_dereplication/ results
                        sub_dir = (
                            input_dir + "/" + entry + "/spectral_dereplication/GNPS/"
                        )

                        if os.path.exists(sub_dir):
                            files = glob.glob(sub_dir + "/*.csv")
                            # print(files)

                            for mz, row in msp.iterrows():
                                # print(msp["id_X"][mz])

                                for fls_g in files:

                                    if msp["id_X"][mz] in fls_g:

                                        gnps_df = pd.read_csv(fls_g)
                                        gnps_df = gnps_df.drop_duplicates(
                                            subset=["GNPSSMILES"]
                                        )

                                        if len(gnps_df) > 0:

                                            for i, row in gnps_df.iterrows():
                                                # if compound name is present

                                                if GNPS_Scoring(gnps_df, i):

                                                    if not isNaN(
                                                        gnps_df["GNPScompound_name"][i]
                                                    ):
                                                        # split if there is a gap in the names

                                                        string_chng = gnps_df[
                                                            "GNPScompound_name"
                                                        ][i].split(" ")

                                                        # create an empty list
                                                        newstr = []

                                                        # for each part of the string in the names
                                                        chng = []

                                                        for j in range(
                                                            len(string_chng)
                                                        ):
                                                            # check if the substrings are present in the matches and no - is present

                                                            if not any(
                                                                x in string_chng[j]
                                                                for x in matches
                                                            ):  # and not '-' == string_chng[j]:

                                                                # IF | and ! not in the substring
                                                                if (
                                                                    "|"
                                                                    not in string_chng[
                                                                        j
                                                                    ]
                                                                    or "!"
                                                                    not in string_chng[
                                                                        j
                                                                    ]
                                                                ):

                                                                    newstr.append(
                                                                        string_chng[j]
                                                                    )
                                                                # if | present in the substring
                                                                elif (
                                                                    "|"
                                                                    in string_chng[j]
                                                                ):

                                                                    # split the string
                                                                    jlen = string_chng[
                                                                        j
                                                                    ].split("|")
                                                                    # how many substrings are left now
                                                                    lst = len(jlen) - 1
                                                                    # append this to chng
                                                                    chng.append(
                                                                        jlen[lst]
                                                                    )
                                                                    break

                                                                    # now append chng to newstr
                                                        chng.append(" ".join(newstr))

                                                        # save this as the correct name
                                                        gnps_df.loc[
                                                            i, "corr_names"
                                                        ] = chng[0]

                                                        if not isNaN(
                                                            gnps_df["GNPSSMILES"][i]
                                                        ):
                                                            if chng == "":
                                                                break
                                                            elif gnps_df["GNPSSMILES"][
                                                                i
                                                            ].isalpha():
                                                                s = pcp.get_compounds(
                                                                    chng[0], "name"
                                                                )
                                                                if s:
                                                                    for comp in s:
                                                                        gnps_df[
                                                                            "GNPSSMILES"
                                                                        ][
                                                                            i
                                                                        ] = (
                                                                            comp.isomeric_smiles
                                                                        )
                                                                else:
                                                                    gnps_df[
                                                                        "GNPSSMILES"
                                                                    ][i] = ""
                                                    else:
                                                        gnps_df["GNPSSMILES"][i] = ""
                                                else:
                                                    gnps_df.drop(
                                                        [i], axis=0, inplace=True
                                                    )

                                            for k, row in gnps_df.iterrows():

                                                if isNaN(gnps_df["GNPSSMILES"][k]):

                                                    if (
                                                        "["
                                                        in gnps_df["GNPScompound_name"][
                                                            k
                                                        ].split(" ")[-1]
                                                    ):
                                                        string_chng = gnps_df[
                                                            "GNPScompound_name"
                                                        ][k].split("[")
                                                        # print(gnps_df['GNPScompound_name'][i])

                                                        # keep_names = []
                                                        for j in range(
                                                            len(string_chng) - 1
                                                        ):
                                                            gnps_df.loc[
                                                                k, "corr_names"
                                                            ] == string_chng[j]
                                                            s = pcp.get_compounds(
                                                                string_chng[j], "name"
                                                            )

                                                            if s:
                                                                for comp in s:
                                                                    gnps_df[
                                                                        "GNPSSMILES"
                                                                    ][
                                                                        k
                                                                    ] = (
                                                                        comp.isomeric_smiles
                                                                    )
                                                                    gnps_df.loc[
                                                                        k, "GNPSformula"
                                                                    ] = (
                                                                        comp.molecular_formula
                                                                    )
                                                                    gnps_df.loc[
                                                                        k, "GNPSinchi"
                                                                    ] = Chem.MolToInchi(
                                                                        Chem.MolFromSmiles(
                                                                            comp.isomeric_smiles
                                                                        )
                                                                    )

                                                            else:
                                                                gnps_df["GNPSSMILES"][
                                                                    k
                                                                ] = ""
                                                                gnps_df.loc[
                                                                    k, "GNPSformula"
                                                                ] = ""
                                                                gnps_df.loc[
                                                                    k, "GNPSinchi"
                                                                ] = ""
                                                if not isNaN(gnps_df["GNPSSMILES"][k]):
                                                    try:
                                                        sx = pcp.get_compounds(
                                                            gnps_df["GNPSSMILES"][k],
                                                            "smiles",
                                                        )
                                                        gnps_df.loc[
                                                            k, "GNPSinchi"
                                                        ] = Chem.MolToInchi(
                                                            Chem.MolFromSmiles(
                                                                comp.isomeric_smiles
                                                            )
                                                        )
                                                        if sx:
                                                            sx = str(sx)
                                                            comp = pcp.Compound.from_cid(
                                                                [
                                                                    int(x)
                                                                    for x in re.findall(
                                                                        r"\b\d+\b", sx
                                                                    )
                                                                ]
                                                            )
                                                            gnps_df.loc[
                                                                k, "GNPSformula"
                                                            ] = comp.molecular_formula

                                                    except Exception:
                                                        gnps_df.loc[
                                                            k, "GNPSformula"
                                                        ] = ""
                                                        gnps_df.loc[k, "GNPSinchi"] = ""

                                        gnps_df = gnps_df.dropna(axis=0, how="all")
                                        csvname = (
                                            (os.path.splitext(fls_g)[0])
                                            + "proc"
                                            + ".csv"
                                        )
                                        gnps_results_csv = csvname.replace(
                                            input_dir, "."
                                        )
                                        msp.loc[
                                            mz, "gnps_results_csv"
                                        ] = gnps_results_csv
                                        gnps_df.to_csv(csvname)
                                        # GNPScsvfiles2.append(csvname)
                                    # dict1 = {'GNPSr': GNPScsvfiles2}
                                    # df = pd.DataFrame(dict1)
                                    # return(df)

                    msp.to_csv(msp_file[0])

                    # HMDB Results
                    if Source == "hmdb" or Source == "all":

                        if not os.path.exists(input_dir + "/hmdb_dframe_str.csv"):

                            # download SDF structures
                            os.system(
                                "wget -P "
                                + input_dir
                                + " https://hmdb.ca/system/downloads/current/structures.zip"
                            )
                            os.system(
                                "unzip "
                                + input_dir
                                + "/structures.zip"
                                + " -d "
                                + input_dir
                            )

                            # Load the sdf
                            dframe = PandasTools.LoadSDF(
                                (input_dir + "/structures.sdf"),
                                idName="HMDB_ID",
                                smilesName="SMILES",
                                molColName="Molecule",
                                includeFingerprints=False,
                            )

                            dframe = dframe[
                                [
                                    "DATABASE_ID",
                                    "SMILES",
                                    "INCHI_IDENTIFIER",
                                    "INCHI_KEY",
                                    "FORMULA",
                                    "MOLECULAR_WEIGHT",
                                    "EXACT_MASS",
                                    "GENERIC_NAME",
                                    "SYNONYMS",
                                ]
                            ]

                        elif os.path.exists(input_dir + "/hmdb_dframe_str.csv"):

                            dframe = pd.read_csv(
                                input_dir + "/hmdb_dframe_str.csv", low_memory=False
                            )

                        # HMDBcsvfiles2 = []
                        # print(entry)
                        # enter the directory with /spectral_dereplication/ results
                        sub_dir = (
                            input_dir + "/" + entry + "/spectral_dereplication/HMDB/"
                        )

                        if os.path.exists(sub_dir):

                            # print(sub_dir)
                            files = glob.glob(sub_dir + "/*.csv")
                            # print(files)
                            for mz, row in msp.iterrows():
                                # print(msp["id_X"][mz])
                                for fls_h in files:
                                    if msp["id_X"][mz] in fls_h:
                                        hmdb_df = pd.read_csv(fls_h)
                                        hmdb_df = hmdb_df.drop_duplicates(
                                            subset=["HMDBcompoundID"]
                                        )

                                        if len(hmdb_df) > 0:
                                            print(entry)
                                            # merge on basis of id, frame and hmdb result files
                                            SmilesHM = pd.merge(
                                                hmdb_df,
                                                dframe,
                                                left_on=hmdb_df.HMDBcompoundID,
                                                right_on=dframe.DATABASE_ID,
                                            )

                                            for i, row in hmdb_df.iterrows():
                                                if HMDB_Scoring(hmdb_df, i):

                                                    for j, row in SmilesHM.iterrows():

                                                        # where index for both match, add the name and SMILES
                                                        if (
                                                            hmdb_df["HMDBcompoundID"][i]
                                                            == SmilesHM[
                                                                "HMDBcompoundID"
                                                            ][j]
                                                        ):
                                                            hmdb_df.loc[
                                                                i, "HMDBSMILES"
                                                            ] = SmilesHM["SMILES"][
                                                                j
                                                            ]  # add SMILES
                                                            hmdb_df.loc[
                                                                i, "HMDBcompound_name"
                                                            ] = SmilesHM[
                                                                "GENERIC_NAME"
                                                            ][
                                                                j
                                                            ]  # add name
                                                            hmdb_df.loc[
                                                                i, "HMDBformula"
                                                            ] = SmilesHM["FORMULA"][
                                                                j
                                                            ]  # add formula
                                                            # hmdb_df.loc[i, 'HMDBinchi'] = Chem.MolToInchi(Chem.MolFromSmiles(SmilesHM['SMILES'][j]))
                                                else:
                                                    hmdb_df.drop(
                                                        [i], axis=0, inplace=True
                                                    )

                                        csvname = (
                                            (os.path.splitext(fls_h)[0])
                                            + "proc"
                                            + ".csv"
                                        )  # name for writing it in a new file
                                        hmdb_results_csv = csvname.replace(
                                            input_dir, "."
                                        )
                                        msp.loc[
                                            mz, "hmdb_results_csv"
                                        ] = hmdb_results_csv
                                        hmdb_df.to_csv(csvname)  # write
                                        # HMDBcsvfiles2.append(csvname)# add to a list
                                    # dict1 = {'HMDBr': HMDBcsvfiles2}
                                    # df = pd.DataFrame(dict1)
                                    # return(df)

                    msp.to_csv(msp_file[0])

                    # MASSBANK Results

                    # enter the directory with /spectral_dereplication/ results
                    if Source == "mbank" or Source == "all":
                        # open another csv path holding empty list, which will be filled
                        # with post processed csv results
                        # MassBankcsvfiles2 = []
                        # print(entry)
                        # enter the directory with /spectral_dereplication/ results
                        sub_dir = (
                            input_dir
                            + "/"
                            + entry
                            + "/spectral_dereplication/MassBank/"
                        )
                        if os.path.exists(sub_dir):
                            files = glob.glob(sub_dir + "/*.csv")
                            # print(files)
                            for mz, row in msp.iterrows():
                                # print(msp["id_X"][mz])
                                for fls_m in files:
                                    if msp["id_X"][mz] in fls_m:
                                        print(fls_m)
                                        mbank_df = pd.read_csv(fls_m)
                                        mbank_df = mbank_df.drop_duplicates(
                                            subset=["MBSMILES"]
                                        )
#                                         if len(mbank_df) > 0:

#                                             for i, row in mbank_df.iterrows():
#                                                 if MB_Scoring(mbank_df, i):

#                                                     inchiK = str(
#                                                         mbank_df["MBinchiKEY"][i]
#                                                     )

#                                                     # extract inchikeys
#                                                     y = pcp.get_compounds(
#                                                         inchiK, "inchikey"
#                                                     )  # compound based on inchikey

#                                                     for compound in y:

#                                                         # add smiles
#                                                         smles = compound.isomeric_smiles
#                                                         mbank_df.loc[
#                                                             i, "MBSMILES"
#                                                         ] = smles
#                                                         # mbank_df.loc[i, 'MBinchi'] =Chem.MolToInchi(Chem.MolFromSmiles(smles))
#                                                 else:
#                                                     mbank_df.drop(
#                                                         [i], axis=0, inplace=True
#                                                     )

                                        csvname = (
                                            (os.path.splitext(fls_m)[0])
                                            + "proc"
                                            + ".csv"
                                        )
                                        mbank_results_csv = csvname.replace(
                                            input_dir, "."
                                        )
                                        msp.loc[
                                            mz, "mbank_results_csv"
                                        ] = mbank_results_csv
                                        mbank_df.to_csv(csvname)
                                        # MassBankcsvfiles2.append(csvname)

                                    # dict1 = {'MBr': MassBankcsvfiles2}
                                    # df = pd.DataFrame(dict1)
                                    # return(df)

                    msp.to_csv(msp_file[0])


# # SIRIUS Post Processing

# In[12]:


def sirius_postproc(input_dir):
#exp_int=0.90, csi_score=-150:

#     def str_can_score(db, i):
#         if (
#             db["explainedIntensity"][i] >= exp_int
#             and db["CSI:FingerIDScore"][i] >= csi_score
#         ):
#             return True
#         else:
#             return False

    # entry is all files and folders in input_dir
    for entry in os.listdir(input_dir):
        # if the entry is also a directory
        if os.path.isdir(os.path.join(input_dir, entry)):
            sub_dir = input_dir + "/" + entry + "/insilico/SIRIUS/"
            msp_csv = input_dir + "/" + entry + "/insilico/MS1DATA.csv"
            if os.path.exists(msp_csv) and os.path.exists(sub_dir):
                # output json files from SIRIUS
                files_S = glob.glob(sub_dir + "/*.json")
                # list of precursor m/z
                msp = pd.read_csv(msp_csv)

                # for each mz
                for mz, row in msp.iterrows():
                    # make a list of files with this mz
                    files_for_mz = []

                    for file in files_S:
                        if str(msp["premz"][mz]) in file:
                            files_for_mz.append(file)

                    # extract the formula and structure files
                    json_dirALL = next(os.walk(files_for_mz[0]))[1]
                    if len(json_dirALL) > 0:
                        sub_sub_dirALL_structure_can = (
                            files_for_mz[0]
                            + "/"
                            + json_dirALL[0]
                            + "/structure_candidates.tsv"
                        )
                        sub_sub_dirALL_formula_can = (
                            files_for_mz[0]
                            + "/"
                            + json_dirALL[0]
                            + "/formula_candidates.tsv"
                        )
                        ALL_Canopus_csv = files_for_mz[0] + "/canopus_summary.tsv"

                        # if both structure files exist
                        if (
                            os.path.exists(sub_sub_dirALL_structure_can)
                            and len(pd.read_csv(sub_sub_dirALL_structure_can, sep="\t")) > 0
                        ):
                            if (
                                os.path.exists(sub_sub_dirALL_formula_can)
                                and len(pd.read_csv(sub_sub_dirALL_formula_can, sep="\t"))
                                > 0
                            ):
                                ALL_structure_csv = pd.read_csv(
                                    sub_sub_dirALL_structure_can, sep="\t"
                                )
                                ALL_formula_csv = pd.read_csv(
                                    sub_sub_dirALL_formula_can, sep="\t"
                                )
                                ALL_Canopus = pd.read_csv(ALL_Canopus_csv, sep="\t")
                                # Add the structure and formula files together
                                for structure, rows in ALL_structure_csv.iterrows():
                                    for formula, rows in ALL_formula_csv.iterrows():
                                        if (
                                            ALL_structure_csv["formulaRank"][structure]
                                            == ALL_formula_csv["rank"][formula]
                                        ):
                                            ALL_structure_csv.loc[
                                                structure, "SiriusScore"
                                            ] = ALL_formula_csv["SiriusScore"][formula]
                                            ALL_structure_csv.loc[
                                                structure, "numExplainedPeaks"
                                            ] = ALL_formula_csv["numExplainedPeaks"][
                                                formula
                                            ]
                                            ALL_structure_csv.loc[
                                                structure, "explainedIntensity"
                                            ] = ALL_formula_csv["explainedIntensity"][
                                                formula
                                            ]
                                            # ALL_structure_csv.loc[structure, "SuspectListEntry"] = "FALSE"
                                            if len(ALL_Canopus) > 0:
                                                if (
                                                    ALL_formula_csv["molecularFormula"][
                                                        formula
                                                    ]
                                                    == ALL_Canopus["molecularFormula"][0]
                                                ):
                                                    ALL_structure_csv.loc[
                                                        structure, "superclass"
                                                    ] = ALL_Canopus["superclass"][0]
                                                    ALL_structure_csv.loc[
                                                        structure, "class"
                                                    ] = ALL_Canopus["class"][0]
                                                    ALL_structure_csv.loc[
                                                        structure, "subclass"
                                                    ] = ALL_Canopus["subclass"][0]

#                                 for str_siriusA, row in ALL_structure_csv.iterrows():
#                                     if not str_can_score(ALL_structure_csv, str_siriusA):
#                                         ALL_structure_csv = ALL_structure_csv.drop(
#                                             str_siriusA, inplace=False
#                                         )

                            result_sirius_name = (
                                sub_dir
                                + "results_for_"
                                + json_dirALL[0].split("_")[-1]
                                + "_"
                                + "structure.csv"
                            )
                            msp.loc[mz, "sirius_result_dir"] = result_sirius_name.replace(
                                input_dir, "."
                            )

                            ALL_structure_csv.to_csv(result_sirius_name)

                        elif not (
                            os.path.exists(sub_sub_dirALL_structure_can)
                            and len(pd.read_csv(sub_sub_dirALL_structure_can, sep="\t")) == 0
                        ):
                            if (
                                os.path.exists(sub_sub_dirALL_formula_can)
                                and len(pd.read_csv(sub_sub_dirALL_formula_can, sep="\t"))
                                > 0
                            ):
                                ALL_formula_csv = pd.read_csv(
                                    sub_sub_dirALL_formula_can, sep="\t"
                                )
                                ALL_Canopus = pd.read_csv(ALL_Canopus_csv, sep="\t")
                                if len(ALL_Canopus) > 0:
                                    for formula, rows in ALL_formula_csv.iterrows():
                                        ALL_formula_csv.loc[
                                            formula, "superclass"
                                        ] = ALL_Canopus["superclass"][0]
                                        ALL_formula_csv.loc[formula, "class"] = ALL_Canopus[
                                            "class"
                                        ][0]
                                        ALL_formula_csv.loc[
                                            formula, "subclass"
                                        ] = ALL_Canopus["subclass"][0]

                                for for_siriusA, row in ALL_formula_csv.iterrows():
                                    if (
                                        not ALL_formula_csv["explainedIntensity"][
                                            for_siriusA
                                        ]
                                        >= exp_int
                                    ):
                                        ALL_formula_csv = ALL_formula_csv.drop(
                                            for_siriusA, inplace=False
                                        )

                                result_sirius_name = (
                                    sub_dir
                                    + "results_for_"
                                    + json_dirALL[0].split("_")[-1]
                                    + "_"
                                    + "formula.csv"
                                )
                                msp.loc[
                                    mz, "sirius_result_dir"
                                ] = result_sirius_name.replace(input_dir, ".")

                                ALL_formula_csv.to_csv(
                                    sub_dir
                                    + "results_for_"
                                    + json_dirALL[0].split("_")[-1]
                                    + "_"
                                    + "formula.csv"
                                )

                            else:
                                print("no file for formula")
                        else:
                            print("no file for structure or formula")
            msp.to_csv(msp_csv)



# In[ ]:


# def str_can_score(db, i):
#     if (
#         db["explainedIntensity"][i] >= exp_int
#         #and db["CSI:FingerIDScore"][i] >= csi_score
#     ):
#         return True
#     else:
#         return False

# # entry is all files and folders in input_dir
# for entry in os.listdir(input_dir):
#     # if the entry is also a directory
#     if os.path.isdir(os.path.join(input_dir, entry)):
#         sub_dir = input_dir + "/" + entry + "/insilico/SIRIUS/"
#         msp_csv = input_dir + "/" + entry + "/insilico/MS1DATA.csv"
#         if os.path.exists(msp_csv) and os.path.exists(sub_dir):
#             # output json files from SIRIUS
#             files_S = glob.glob(sub_dir + "/*.json")
#             # list of precursor m/z
#             msp = pd.read_csv(msp_csv)

#             # for each mz
#             for mz, row in msp.iterrows():
#                 # make a list of files with this mz
#                 #print(mz)
#                 files_for_mz = []

#                 for file in files_S:
#                     if str(msp["premz"][mz]) in file:
#                         files_for_mz.append(file)

#                 # extract the formula and structure files
#                 json_dirALL = next(os.walk(files_for_mz[0]))[1]
#                 if len(json_dirALL)>0:
#                     sub_sub_dirALL_structure_can = (
#                         files_for_mz[0]
#                         + "/"
#                         + json_dirALL[0]
#                         + "/structure_candidates.tsv"
#                     )
#                     sub_sub_dirALL_formula_can = (
#                         files_for_mz[0]
#                         + "/"
#                         + json_dirALL[0]
#                         + "/formula_candidates.tsv"
#                     )
#                     ALL_Canopus_csv = files_for_mz[0] + "/canopus_summary.tsv"

#                     # if both structure files exist
#                     if (
#                         os.path.exists(sub_sub_dirALL_structure_can)
#                         and len(pd.read_csv(sub_sub_dirALL_structure_can, sep="\t")) > 0
#                     ):
#                         if (
#                             os.path.exists(sub_sub_dirALL_formula_can)
#                             and len(pd.read_csv(sub_sub_dirALL_formula_can, sep="\t"))
#                             > 0
#                         ):
#                             ALL_structure_csv = pd.read_csv(
#                                 sub_sub_dirALL_structure_can, sep="\t"
#                             )
#                             ALL_formula_csv = pd.read_csv(
#                                 sub_sub_dirALL_formula_can, sep="\t"
#                             )
#                             ALL_Canopus = pd.read_csv(ALL_Canopus_csv, sep="\t")
#                             # Add the structure and formula files together
#                             for structure, rows in ALL_structure_csv.iterrows():
#                                 for formula, rows in ALL_formula_csv.iterrows():
#                                     if (
#                                         ALL_structure_csv["formulaRank"][structure]
#                                         == ALL_formula_csv["rank"][formula]
#                                     ):
#                                         ALL_structure_csv.loc[
#                                             structure, "SiriusScore"
#                                         ] = ALL_formula_csv["SiriusScore"][formula]
#                                         ALL_structure_csv.loc[
#                                             structure, "numExplainedPeaks"
#                                         ] = ALL_formula_csv["numExplainedPeaks"][
#                                             formula
#                                         ]
#                                         ALL_structure_csv.loc[
#                                             structure, "explainedIntensity"
#                                         ] = ALL_formula_csv["explainedIntensity"][
#                                             formula
#                                         ]
#                                         # ALL_structure_csv.loc[structure, "SuspectListEntry"] = "FALSE"
#                                         if len(ALL_Canopus) > 0:
#                                             if (
#                                                 ALL_formula_csv["molecularFormula"][
#                                                     formula
#                                                 ]
#                                                 == ALL_Canopus["molecularFormula"][0]
#                                             ):
#                                                 ALL_structure_csv.loc[
#                                                     structure, "superclass"
#                                                 ] = ALL_Canopus["superclass"][0]
#                                                 ALL_structure_csv.loc[
#                                                     structure, "class"
#                                                 ] = ALL_Canopus["class"][0]
#                                                 ALL_structure_csv.loc[
#                                                     structure, "subclass"
#                                                 ] = ALL_Canopus["subclass"][0]

#                             for str_siriusA, row in ALL_structure_csv.iterrows():
#                                 if not str_can_score(ALL_structure_csv, str_siriusA):
#                                     ALL_structure_csv = ALL_structure_csv.drop(
#                                         str_siriusA, inplace=False
#                                     )

#                         result_sirius_name = (
#                             sub_dir
#                             + "results_for_"
#                             + json_dirALL[0].split("_")[-1]
#                             + "_"
#                             + "structure.csv"
#                         )
#                         msp.loc[mz, "sirius_result_dir"] = result_sirius_name.replace(
#                             input_dir, "."
#                         )

#                         ALL_structure_csv.to_csv(result_sirius_name)

#                     elif not (
#                         os.path.exists(sub_sub_dirALL_structure_can)
#                         and len(pd.read_csv(sub_sub_dirALL_structure_can, sep="\t")) > 0
#                     ):
#                         if (
#                             os.path.exists(sub_sub_dirALL_formula_can)
#                             and len(pd.read_csv(sub_sub_dirALL_formula_can, sep="\t"))
#                             > 0
#                         ):
#                             ALL_formula_csv = pd.read_csv(
#                                 sub_sub_dirALL_formula_can, sep="\t"
#                             )
#                             ALL_Canopus = pd.read_csv(ALL_Canopus_csv, sep="\t")
#                             if len(ALL_Canopus) > 0:
#                                 for formula, rows in ALL_formula_csv.iterrows():
#                                     ALL_formula_csv.loc[
#                                         formula, "superclass"
#                                     ] = ALL_Canopus["superclass"][0]
#                                     ALL_formula_csv.loc[formula, "class"] = ALL_Canopus[
#                                         "class"
#                                     ][0]
#                                     ALL_formula_csv.loc[
#                                         formula, "subclass"
#                                     ] = ALL_Canopus["subclass"][0]

#                             for for_siriusA, row in ALL_formula_csv.iterrows():
#                                 if (
#                                     not ALL_formula_csv["explainedIntensity"][
#                                         for_siriusA
#                                     ]
#                                     >= exp_int
#                                 ):
#                                     ALL_formula_csv = ALL_formula_csv.drop(
#                                         for_siriusA, inplace=False
#                                     )

#                             result_sirius_name = (
#                                 sub_dir
#                                 + "results_for_"
#                                 + json_dirALL[0].split("_")[-1]
#                                 + "_"
#                                 + "formula.csv"
#                             )
#                             msp.loc[
#                                 mz, "sirius_result_dir"
#                             ] = result_sirius_name.replace(input_dir, ".")

#                             ALL_formula_csv.to_csv(
#                                 sub_dir
#                                 + "results_for_"
#                                 + json_dirALL[0].split("_")[-1]
#                                 + "_"
#                                 + "formula.csv"
#                             )

#                         else:
#                             print("no file for formula")
#                     else:
#                         print("no file for structure or formula")
#         msp.to_csv(msp_csv)



# In[8]:


def MCSS_for_SpecDB(input_dir, Source):
    # Describe the heavy atoms to be considered for MCSS
    heavy_atoms = ["C", "N", "P", "O", "S"]
    # list all files and directories
    for entry in os.listdir(input_dir):

        if os.path.isdir(os.path.join(input_dir, entry)):

            # for specdb
            specdb_msp_file = glob.glob(
                input_dir + "/" + entry + "/spectral_dereplication" + "/*.csv"
            )

            if len(specdb_msp_file) > 0:

                if os.path.exists(specdb_msp_file[0]):

                    spec_msp = pd.read_csv(specdb_msp_file[0])

                    for mz, row in spec_msp.iterrows():

                        if Source == "gnps" or Source == "specdb" or Source == "all":

                            sub_dir = (
                                input_dir
                                + "/"
                                + entry
                                + "/spectral_dereplication/GNPS/"
                            )
                            if os.path.exists(sub_dir):
                                gnps_files = glob.glob(sub_dir + "/*proc.csv")

                                for files in gnps_files:
                                    if spec_msp["id_X"][mz] in files:
                                        gnpsproc = pd.read_csv(files)

                                        if len(gnpsproc) > 0:
                                            G_Smiles = gnpsproc["GNPSSMILES"]
                                            G_Smiles = list(filter(None, G_Smiles))
                                            # print(G_Smiles)
                                            # create empty list of GNPS top smiles
                                            GNPS_Mol = []
                                            # extract only the InChI of the top 5
                                            for j in list(G_Smiles):
                                                if not isNaN(j):
                                                    print(type(j))
                                                    mol2 = Chem.MolFromSmiles(j)
                                                    GNPS_Mol.append(mol2)

                                            if len(GNPS_Mol) >= 2:
                                                res = rdFMCS.FindMCS(GNPS_Mol, timeout=60)
                                                sm_res = Chem.MolToSmiles(
                                                    Chem.MolFromSmarts(res.smartsString)
                                                )
                                                # if there are atleast 3 heavy atoms in the MCSS, then add it to the result file
                                                elem = [
                                                    ele
                                                    for ele in heavy_atoms
                                                    if (ele in sm_res)
                                                ]
                                                if elem and len(sm_res) >= 3:
                                                    spec_msp.loc[
                                                        mz, "GNPS_MCSSstring"
                                                    ] = res.smartsString
                                                    spec_msp.loc[
                                                        mz, "GNPS_MCSS_SMILES"
                                                    ] = Chem.MolToSmiles(
                                                        Chem.MolFromSmarts(
                                                            res.smartsString
                                                        )
                                                    )
                        if Source == "hmdb" or Source == "specdb" or Source == "all":
                            sub_dir = (
                                input_dir
                                + "/"
                                + entry
                                + "/spectral_dereplication/HMDB/"
                            )
                            if os.path.exists(sub_dir):
                                hmdb_files = glob.glob(sub_dir + "/*proc.csv")

                                for files in hmdb_files:
                                    if spec_msp["id_X"][mz] in files:
                                        hmdbproc = pd.read_csv(files)

                                        if len(hmdbproc) > 0:
                                            H_Smiles = hmdbproc["HMDBSMILES"]
                                            H_Smiles = list(filter(None, H_Smiles))

                                            HMDB_Mol = []
                                            # extract only the InChI of the top 5
                                            for j in list(H_Smiles):
                                                if not isNaN(j):

                                                    mol2 = Chem.MolFromSmiles(j)
                                                    HMDB_Mol.append(mol2)

                                            if len(HMDB_Mol) >= 2:
                                                res = rdFMCS.FindMCS(HMDB_Mol, timeout=60)
                                                sm_res = Chem.MolToSmiles(
                                                    Chem.MolFromSmarts(res.smartsString)
                                                )
                                                # if there are atleast 3 heavy atoms in the MCSS, then add it to the result file
                                                elem = [
                                                    ele
                                                    for ele in heavy_atoms
                                                    if (ele in sm_res)
                                                ]
                                                if elem and len(sm_res) >= 3:
                                                    spec_msp.loc[
                                                        mz, "HMDB_MCSSstring"
                                                    ] = res.smartsString
                                                    spec_msp.loc[
                                                        mz, "HMDB_MCSS_SMILES"
                                                    ] = Chem.MolToSmiles(
                                                        Chem.MolFromSmarts(
                                                            res.smartsString
                                                        )
                                                    )
                        if Source == "mbank" or Source == "specdb" or Source == "all":
                            sub_dir = (
                                input_dir
                                + "/"
                                + entry
                                + "/spectral_dereplication/MassBank/"
                            )
                            if os.path.exists(sub_dir):
                                mbank_files = glob.glob(sub_dir + "/*proc.csv")

                                for files in mbank_files:
                                    if spec_msp["id_X"][mz] in files:
                                        mbankproc = pd.read_csv(files)

                                        if len(mbankproc) > 0:
                                            M_Smiles = mbankproc["MBSMILES"]
                                            M_Smiles = list(filter(None, M_Smiles))

                                            MB_Mol = []
                                            # extract only the InChI of the top 5
                                            for j in list(M_Smiles):
                                                if not isNaN(j):

                                                    mol2 = Chem.MolFromSmiles(j)
                                                    MB_Mol.append(mol2)

                                            if len(MB_Mol) >= 2:
                                                res = rdFMCS.FindMCS(MB_Mol, timeout=60)
                                                sm_res = Chem.MolToSmiles(
                                                    Chem.MolFromSmarts(res.smartsString)
                                                )
                                                # if there are atleast 3 heavy atoms in the MCSS, then add it to the result file
                                                elem = [
                                                    ele
                                                    for ele in heavy_atoms
                                                    if (ele in sm_res)
                                                ]
                                                if elem and len(sm_res) >= 3:
                                                    spec_msp.loc[
                                                        mz, "MB_MCSSstring"
                                                    ] = res.smartsString
                                                    spec_msp.loc[
                                                        mz, "MB_MCSS_SMILES"
                                                    ] = Chem.MolToSmiles(
                                                        Chem.MolFromSmarts(
                                                            res.smartsString
                                                        )
                                                    )
                    spec_msp.to_csv(specdb_msp_file[0])


# In[9]:


def MCSS_for_SIRIUS(input_dir):


    # Describe the heavy atoms to be considered for MCSS
    heavy_atoms = ["C", "N", "P", "O", "S"]
    # list all files and directories
    for entry in os.listdir(input_dir):
        if os.path.isdir(os.path.join(input_dir, entry)):
            # for sirius
            sirius_msp_csv = input_dir + "/" + entry + "/insilico/MS1DATA.csv"
            sub_dir = input_dir + "/" + entry + "/insilico/SIRIUS/"
            if os.path.exists(sirius_msp_csv) and os.path.exists(sub_dir):
                sirius_msp = pd.read_csv(sirius_msp_csv)
                sirius_files = glob.glob(sub_dir)
                for sir_file in sirius_files:
                    r = [s for s in os.listdir(sir_file) if "structure" in s]
                    for filenames in r:
                        sirius_f = sir_file + filenames
                        for mz, row in sirius_msp.iterrows():
                            if str(sirius_msp["id_X"][mz].split("_")[1]) in sirius_f:
                                s_f = pd.read_csv(str(sirius_f))
                                if (
                                    len(s_f) > 0
                                    and "smiles" in s_f.columns.values.tolist()
                                ):
                                    S_Smiles = s_f["smiles"]
                                    # create empty list of MB top smiles
                                    SIRIUS_Mol = []

                                    # extract only the InChI of the top 5
                                    for j in list(S_Smiles):
                                        mol2 = Chem.MolFromSmiles(j)
                                        SIRIUS_Mol.append(mol2)
                                    if len(SIRIUS_Mol) >= 2:
                                        res = rdFMCS.FindMCS(SIRIUS_Mol, timeout=60)
                                        sm_res = Chem.MolToSmiles(
                                            Chem.MolFromSmarts(res.smartsString)
                                        )
                                        # if there are atleast 3 heavy atoms in the MCSS, then add it to the result file
                                        elem = [
                                            ele
                                            for ele in heavy_atoms
                                            if (ele in sm_res)
                                        ]
                                        if elem and len(sm_res) >= 3:
                                            sirius_msp.loc[
                                                mz, "SIRIUS_MCSSstring"
                                            ] = res.smartsString
                                            sirius_msp.loc[
                                                mz, "SIRIUS_MCSS_SMILES"
                                            ] = Chem.MolToSmiles(
                                                Chem.MolFromSmarts(res.smartsString)
                                            )
                sirius_msp.to_csv(sirius_msp_csv)


# In[10]:


def SuspectListScreening(input_dir, SuspectListPath, tanimoto, Source):
    def isNaN(string):
        return string != string

    SuspectList = pd.read_csv(SuspectListPath)

    for entry in os.listdir(input_dir):

        if os.path.isdir(os.path.join(input_dir, entry)):

            if Source == "gnps" or Source == "specdb" or Source == "all":

                sub_dir = input_dir + "/" + entry + "/spectral_dereplication/GNPS/"
                if os.path.exists(sub_dir):
                    gnps_files = glob.glob(sub_dir + "/*proc.csv")
                    for file in gnps_files:
                        gnpsproc = pd.read_csv(file)
                        if len(gnpsproc) > 0:
                            for g, row in gnpsproc.iterrows():
                                for s, row in SuspectList.iterrows():
                                    if (
                                        not isNaN(gnpsproc["GNPSSMILES"][g])
                                        and gnpsproc["GNPSSMILES"][g] != " "
                                    ):
                                        if (
                                            not isNaN(SuspectList["SMILES"][s])
                                            and SuspectList["SMILES"][s] != " "
                                        ):
                                            LHms2 = [
                                                Chem.MolFromSmiles(
                                                    gnpsproc["GNPSSMILES"][g]
                                                ),
                                                Chem.MolFromSmiles(
                                                    SuspectList["SMILES"][s]
                                                ),
                                            ]
                                            LHfps2 = [
                                                AllChem.GetMorganFingerprintAsBitVect(
                                                    x2, 2, nBits=2048
                                                )
                                                for x2 in LHms2
                                            ]
                                            LHtn2 = DataStructs.FingerprintSimilarity(
                                                LHfps2[0], LHfps2[1]
                                            )
                                            if LHtn2 >= tanimoto:
                                                gnpsproc.loc[
                                                    g, "SLGsmiles"
                                                ] = SuspectList["SMILES"][s]
                                                gnpsproc.loc[
                                                    g, "SLGname"
                                                ] = SuspectList["Name"][s]
                                                gnpsproc.loc[g, "SLGtanimoto"] = LHtn2
                        gnpsproc.to_csv(file)
                        return gnpsproc
            if Source == "hmdb" or Source == "specdb" or Source == "all":

                sub_dir = input_dir + "/" + entry + "/spectral_dereplication/HMDB/"
                if os.path.exists(sub_dir):
                    hmdb_files = glob.glob(sub_dir + "/*proc.csv")
                    for file in hmdb_files:

                        hmdbproc = pd.read_csv(file)
                        if len(hmdbproc) > 0:
                            for h, row in hmdbproc.iterrows():
                                for s, row in SuspectList.iterrows():
                                    if (
                                        not isNaN(hmdbproc["HMDBSMILES"][h])
                                        and hmdbproc["HMDBSMILES"][h] != " "
                                    ):
                                        if (
                                            not isNaN(SuspectList["SMILES"][s])
                                            and SuspectList["SMILES"][s] != " "
                                        ):
                                            LHms2 = [
                                                Chem.MolFromSmiles(
                                                    hmdbproc["HMDBSMILES"][h]
                                                ),
                                                Chem.MolFromSmiles(
                                                    SuspectList["SMILES"][s]
                                                ),
                                            ]
                                            LHfps2 = [
                                                AllChem.GetMorganFingerprintAsBitVect(
                                                    x2, 2, nBits=2048
                                                )
                                                for x2 in LHms2
                                            ]
                                            LHtn2 = DataStructs.FingerprintSimilarity(
                                                LHfps2[0], LHfps2[1]
                                            )
                                            if LHtn2 >= tanimoto:
                                                hmdbproc.loc[
                                                    h, "SLHsmiles"
                                                ] = SuspectList["SMILES"][s]
                                                hmdbproc.loc[
                                                    h, "SLHname"
                                                ] = SuspectList["Name"][s]
                                                hmdbproc.loc[h, "SLHtanimoto"] = LHtn2

                        hmdbproc.to_csv(file)
                        return hmdbproc

            if Source == "mbank" or Source == "specdb" or Source == "all":

                sub_dir = input_dir + "/" + entry + "/spectral_dereplication/MassBank/"
                if os.path.exists(sub_dir):
                    mbank_files = glob.glob(sub_dir + "/*proc.csv")
                    for file in mbank_files:

                        mbankproc = pd.read_csv(file)

                        if len(mbankproc) > 0:
                            for m, row in mbankproc.iterrows():
                                for s, row in SuspectList.iterrows():
                                    if (
                                        not isNaN(mbankproc["MBSMILES"][m])
                                        and mbankproc["MBSMILES"][m] != " "
                                    ):
                                        if (
                                            not isNaN(SuspectList["SMILES"][s])
                                            and SuspectList["SMILES"][s] != " "
                                        ):
                                            LHms2 = [
                                                Chem.MolFromSmiles(
                                                    mbankproc["MBSMILES"][m]
                                                ),
                                                Chem.MolFromSmiles(
                                                    SuspectList["SMILES"][s]
                                                ),
                                            ]
                                            LHfps2 = [
                                                AllChem.GetMorganFingerprintAsBitVect(
                                                    x2, 2, nBits=2048
                                                )
                                                for x2 in LHms2
                                            ]
                                            LHtn2 = DataStructs.FingerprintSimilarity(
                                                LHfps2[0], LHfps2[1]
                                            )
                                            if LHtn2 >= tanimoto:
                                                mbankproc.loc[
                                                    m, "SLMsmiles"
                                                ] = SuspectList["SMILES"][s]
                                                mbankproc.loc[
                                                    m, "SLMname"
                                                ] = SuspectList["Name"][s]
                                                mbankproc.loc[m, "SLMtanimoto"] = LHtn2

                        mbankproc.to_csv(file)
                        return mbankproc


# In[11]:


# Candidate Selection
def chemMN_CandidateSelection(df, tn_sim=0.85):

    """chemMN_CandidateSelection function is used to generate a Cytoscape readable tsv file.
    This file contains start(starting SMILES) and end(target SMILES) nodes and the tanimoto
    similarity scores between the nodes. User can visualize the structural similarity
    between the given SMILES. It provides an "ALL against ALL" network.

    Parameters:
    df: dataframe that contains "SMILES", "ranks", "Source". This function is specifically for
    candidate selection and so these columns are necessary.


    Returns:
    dataframe: it returns a df with following columns to be loaded into Cytoscape.
    1. Start, starting node/SMILES
    2. End, ending node/SMILES
    3. Tanimoto, Tanimoto between Start and End node
    4. Start_SMILES
    5. End_SMILES
    6. Start_Source
    7. End_Source
    8. MCSS, Maximum Common Substructure between start and end node/SMILES
    9. sorted_row, contains ids of the start and end nodes as a list


    Usage:
    chemMN_CandidateSelection(df)

    """

    # define an empty variable
    # one_df = []
    # define empty variable to save the edges
    dbn = []
    # for each entry in the df
    for i, row in df.iterrows():
        # to compare each element with each other element of the df
        for j, row in df.iterrows():
            try:
                # calcultae tanimoto
                ms = [
                    Chem.MolFromSmiles(df["SMILES"][i]),
                    Chem.MolFromSmiles(df["SMILES"][j]),
                ]
                fps = [
                    AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=2048) for x in ms
                ]
                tn = DataStructs.FingerprintSimilarity(fps[0], fps[1])
                if tn>= tn_sim:
                    if df["SMILES"][i]!= df["SMILES"][j]:
                        # save all entries to a matrix
                        dbn.append(
                            {
                                "Name_i": df["ranks"][i],
                                "Name_j": df["ranks"][j],
                                "i": df["SMILES"][i],
                                "j": df["SMILES"][j],
                                "Source_i": df["Source"][i],
                                "Source_j": df["Source"][j],
                                "Tanimoto": tn,
                            }
                        )

            except Exception:
                # print(e.string)
                pass

    # save chemical similarities
    db_edgenode = pd.DataFrame(dbn)
    
    if not db_edgenode.empty:

        # another empty variable to store the results for final tsv file
        dfe = []

        # heavy atoms for MCSS Calculation
        heavy_atoms = ["C", "N", "P", "O", "S"]

        # for the previous dataframe
        for i, row in db_edgenode.iterrows():
            # if the tanimoto > 0.85 for high similarity
            if db_edgenode["Tanimoto"][i] >= tn_sim:

                # calculate MCSS
                n = [
                    Chem.MolFromSmiles(db_edgenode["i"][i]),
                    Chem.MolFromSmiles(db_edgenode["j"][i]),
                ]
                res = rdFMCS.FindMCS(n, timeout=60)
                sm_res = Chem.MolToSmiles(Chem.MolFromSmarts(res.smartsString))

                # Check if the MCSS has one of the heavy atoms and whether they are
                # more than 3
                elem = [ele for ele in heavy_atoms if (ele in sm_res)]
                if elem and len(sm_res) >= 3:
                    MCSS_SMILES = Chem.MolToSmiles(Chem.MolFromSmarts(res.smartsString))

                # save everything into a dataframe
                dfe.append(
                    {
                        "Start": db_edgenode["Name_i"][i],
                        "End": db_edgenode["Name_j"][i],
                        "Tanimoto": db_edgenode["Tanimoto"][i],
                        "Start_SMILES": db_edgenode["i"][i],
                        "End_SMILES": db_edgenode["j"][i],
                        "Start_Source": db_edgenode["Source_i"][i],
                        "End_Source": db_edgenode["Source_j"][i],
                        "MCSS": MCSS_SMILES,
                    }
                )
        df_edge = pd.DataFrame(dfe)
        # generate a column called sorted_row which contains ids of the start and end nodes as a list
        df_edge["Start"] = df_edge["Start"].astype(str)
        df_edge["End"] = df_edge["End"].astype(str)
        df_edge["sorted_row"] = [sorted([a, b]) for a, b in zip(df_edge.Start, df_edge.End)]
        df_edge["sorted_row"] = df_edge["sorted_row"].astype(str)
        df_edge.drop_duplicates(subset=["sorted_row"], inplace=True)

        return df_edge


# In[12]:


def one_candidate_selection(
    df,
    Source="SGHM",
    tn_ident=0.99,
    sirius_df=None,
    mbank_df=None,
    gnps_df=None,
    hmdb_df=None,
):

    """one_candidate_selection function is used to generate a dataframe that tells,
    for each candidate SMILES, what was the source or how many sources had the same
    candidate. The idea is to merge all candidate SMILES into one list, preserving
    the rank and source, and then checking whether these SMILES come from SIRIUS or
    any spectral DB. If a SMILE is repeated in more sources, its confidence score
    increases and is considered the most likely candidate structure. This function
    is not stand-alone and is part of the function CandidateSelection_SimilarityandIdentity


    Parameters:
    df: dataframe that contains "SMILES", "ranks", "Source". This function is
    specifically for candidate selection and so these columns are necessary.
    Source: this depends on how many sources were used. Possiblilities are:
    1. SGHM (all)
    2. SGM (SIRIUS, GNPS, MassBank)
    3. SHM (SIRIUS, HMDB, MassBank)
    4. SGH (SIRIUS, GNPS, HMDB)
    5. GHM (GNPS, HMDB, MassBank)
    6. SG (SIRIUS, GNPS)
    7. SH (SIRIUS, HMDB)
    8. SM (SIRIUS, MassBank)
    9. GM (GNPS, MassBank)
    10. GH (GNPS, HMDB)
    11. HM (HMDB, MassBank)
    12. S
    13. G
    14. H
    15. M

    Returns:
    dataframe: it returns a df with follwoing columns which can be used to
    prioritize a database for the final candidate selection.
    1. Source, contains name of the source (SIRIUS, GNPS, HMDB or MassBank)
    2. ranks, contains first letter of the source and a rank number seperated
    by _ e.g: G_1(GNPS, 1st rank)
    3. SMILES
    4. SIRIUS, the rank again but only when the corresponding row SMILES is
    also part of SIRIUS results
    5. GNPS , same as SIRIUS but for GNPS
    5. MassBank
    6. HMDB


    Usage:
    chemMN_CandidateSelection(df, Source = "SGHM")

    """

    # define empty columns for each Source to only fill if the corresponding
    # SMILES is also present in the source

    df["SIRIUS"] = np.nan
    df["GNPS"] = np.nan
    df["MassBank"] = np.nan
    df["HMDB"] = np.nan

    # for each SMILES in df
    for smiles, rows in df.iterrows():

        # If the source contains SIRIUS
        if (
            Source == "SGHM"
            or Source == "SGM"
            or Source == "SGH"
            or Source == "SHM"
            or Source == "SG"
            or Source == "SM"
            or Source == "SH"
            or Source == "S"
        ):
            # sirius_df comes from within the function CandidateSelection_SimilarityandIdentity
            for sirius_i, row in sirius_df.iterrows():
                # calculate tanimoto
                try:
                    ms = [
                        Chem.MolFromSmiles(df["SMILES"][smiles]),
                        Chem.MolFromSmiles(sirius_df["smiles"][sirius_i]),
                    ]
                    fps = [
                        AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=2048)
                        for x in ms
                    ]
                    tn = DataStructs.FingerprintSimilarity(fps[0], fps[1])
                    # since we are dealing with idenity here so tanimoto of 0.99 is appropriate
                    if tn >= tn_ident:

                        # if SIRIUS is blank, add the SIRIUS id
                        if isNaN(df["SIRIUS"][smiles]):

                            df.loc[smiles, "SIRIUS"] = sirius_df["rank_ids"][sirius_i]
                        # if not empty, add SIRIUS id, with a comma
                        else:
                            df.loc[smiles, "SIRIUS"] = (
                                str(df["SIRIUS"][smiles])
                                + ", "
                                + sirius_df["rank_ids"][sirius_i]
                            )

                except Exception:
                    # print(e.string)
                    pass

        # If the Source contains GNPS
        if (
            Source == "SGHM"
            or Source == "SGM"
            or Source == "SGH"
            or Source == "GHM"
            or Source == "SG"
            or Source == "GM"
            or Source == "GH"
            or Source == "G"
        ):

            # gnps_df comes from within the function CandidateSelection_SimilarityandIdentity
            for gnps_i, row in gnps_df.iterrows():
                try:
                    # calculate tanimoto
                    ms = [
                        Chem.MolFromSmiles(df["SMILES"][smiles]),
                        Chem.MolFromSmiles(gnps_df["GNPSSMILES"][gnps_i]),
                    ]
                    fps = [
                        AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=2048)
                        for x in ms
                    ]
                    tn = DataStructs.FingerprintSimilarity(fps[0], fps[1])

                    # since we are dealing with idenity here so tanimoto of 0.99 is appropriate
                    if tn >= tn_ident:

                        # if GNPS is blank, add the GNPS id
                        if isNaN(df["GNPS"][smiles]):

                            df.loc[smiles, "GNPS"] = gnps_df["rank_ids"][gnps_i]
                        # if not empty, add GNPS id, with a comma
                        else:
                            df.loc[smiles, "GNPS"] = (
                                str(df["GNPS"][smiles])
                                + ", "
                                + gnps_df["rank_ids"][gnps_i]
                            )

                except Exception:
                    # print(e.string)
                    pass

        # If the source contains HMDB
        if (
            Source == "SGHM"
            or Source == "SGH"
            or Source == "SHM"
            or Source == "GHM"
            or Source == "SH"
            or Source == "GH"
            or Source == "HM"
            or Source == "H"
        ):
            for hmdb_i, row in hmdb_df.iterrows():
                try:
                    # calculate tanimoto
                    ms = [
                        Chem.MolFromSmiles(df["SMILES"][smiles]),
                        Chem.MolFromSmiles(hmdb_df["HMDBSMILES"][hmdb_i]),
                    ]
                    fps = [
                        AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=2048)
                        for x in ms
                    ]
                    tn = DataStructs.FingerprintSimilarity(fps[0], fps[1])

                    # since we are dealing with idenity here so tanimoto of 0.99 is appropriate
                    if tn >= tn_ident:

                        # if HMDB is blank, add the HMDB id
                        if isNaN(df["HMDB"][smiles]):

                            df.loc[smiles, "HMDB"] = hmdb_df["rank_ids"][hmdb_i]
                        # if not empty, add HMDB id, with a comma
                        else:
                            df.loc[smiles, "HMDB"] = (
                                str(df["HMDB"][smiles])
                                + ", "
                                + hmdb_df["rank_ids"][hmdb_i]
                            )

                except Exception:
                    # print(e.string)
                    pass


        # If the source contains MassBank
        if (
            Source == "SGHM"
            or Source == "SGM"
            or Source == "SHM"
            or Source == "GHM"
            or Source == "SM"
            or Source == "GM"
            or Source == "HM"
            or Source == "M"
        ):
            # mbank_df comes from within the function CandidateSelection_SimilarityandIdentity
            for mbank_i, row in mbank_df.iterrows():
                try:
                    # calculate tanimoto
                    ms = [
                        Chem.MolFromSmiles(df["SMILES"][smiles]),
                        Chem.MolFromSmiles(mbank_df["MBSMILES"][mbank_i]),
                    ]
                    fps = [
                        AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=2048)
                        for x in ms
                    ]
                    tn = DataStructs.FingerprintSimilarity(fps[0], fps[1])

                    # since we are dealing with idenity here so tanimoto of 0.99 is appropriate
                    if tn >= tn_ident:

                        # if MassBank is blank, add the MassBank id
                        if isNaN(df["MassBank"][smiles]):

                            df.loc[smiles, "MassBank"] = mbank_df["rank_ids"][mbank_i]
                        # if not empty, add MassBank id, with a comma
                        else:
                            df.loc[smiles, "MassBank"] = (
                                str(df["MassBank"][smiles])
                                + ", "
                                + mbank_df["rank_ids"][mbank_i]
                            )

                except Exception:
                    # print(e.string)
                    pass
    return df


# In[13]:


def add_count_column(df_one_candidate):
    df_one_candidate = df_one_candidate.dropna(axis=0, how="all", subset = ['SIRIUS', 'GNPS', 'MassBank', 'HMDB'])
    # create new df only with the Sources column
    df = pd.DataFrame(
        {
            "SIRIUS": df_one_candidate["SIRIUS"],
            "GNPS": df_one_candidate["GNPS"],
            "MassBank": df_one_candidate["MassBank"],
            "HMDB": df_one_candidate["HMDB"],
        }
    )

    # df_one_candidate = df_one_candidate.dropna(subset=["SIRIUS", "GNPS", "HMDB", "MassBank"], how='all', inplace=True)

    index_SIRIUS = [x for x, row in df.iterrows() if not isNaN(df["SIRIUS"][x])]
    index_GNPS = [x for x, row in df.iterrows() if not isNaN(df["GNPS"][x])]
    index_MassBank = [x for x, row in df.iterrows() if not isNaN(df["MassBank"][x])]
    index_HMDB = [x for x, row in df.iterrows() if not isNaN(df["HMDB"][x])]

    # make a list of the rows
    list_of_indices = [index_SIRIUS] + [index_GNPS] + [index_MassBank] + [index_HMDB]
    length_of_list = len([idx for idx, x in enumerate(list_of_indices)  if x])
    if length_of_list == 1:
        #print("1 is correct")
        # change number of counts, since its 1 source, so all counts should be 1
        df_one_candidate["Count"] = 1
        # extract ranks numbers from ranks
        df_one_candidate["rank_num"] = [counts.split("_")[1] for counts in df_one_candidate["ranks"]]
        # convert any str to int
        df_one_candidate["rank_num"] = [int(x) for x in df_one_candidate["rank_num"]]
        #Sort by ranknum
        df_one_candidate = df_one_candidate.sort_values(
            by="rank_num", ascending=False
        )
        for r, rows in df_one_candidate.iterrows():
            if df_one_candidate["Source"][r]=="GNPS":
                df_one_candidate.loc[r, "rank_db"] = 1
            if df_one_candidate["Source"][r]=="SIRIUS":
                df_one_candidate.loc[r, "rank_db"] = 2
            if df_one_candidate["Source"][r]=="MassBank":
                df_one_candidate.loc[r, "rank_db"] = 3
            if df_one_candidate["Source"][r]=="HMDB":
                df_one_candidate.loc[r, "rank_db"] = 4
        df_one_candidate.sort_values(
            by=["rank_num", "rank_db"], ascending=[True, True],
            inplace=True)

        return df_one_candidate
    elif length_of_list > 1:
        # make a list of the rows
        list_of_indices = index_SIRIUS + index_GNPS + index_MassBank + index_HMDB
        # count how many times one of the rows is appearing and add count
        count_list = [[x, list_of_indices.count(x)] for x in set(list_of_indices)]
        # add this info to one_can
        df_one_candidate["Count"] = [count_list[x][1] for x in range(len(count_list))]
        # sort the list by count in descending order
        sorted_count_one_candidate = df_one_candidate.sort_values(
            by="Count", ascending=False
        )
        sorted_count_one_candidate["rank_num"] = [counts.split("_")[1] for counts in sorted_count_one_candidate["ranks"]]
        sorted_count_one_candidate["rank_num"] = [int(x) for x in sorted_count_one_candidate["rank_num"]]
        for r, rows in sorted_count_one_candidate.iterrows():
            if sorted_count_one_candidate["Source"][r]=="GNPS":
                sorted_count_one_candidate.loc[r, "rank_db"] = 1
            if sorted_count_one_candidate["Source"][r]=="SIRIUS":
                sorted_count_one_candidate.loc[r, "rank_db"] = 2
            if sorted_count_one_candidate["Source"][r]=="MassBank":
                sorted_count_one_candidate.loc[r, "rank_db"] = 3
            if sorted_count_one_candidate["Source"][r]=="HMDB":
                sorted_count_one_candidate.loc[r, "rank_db"] = 4

        sorted_count_one_candidate.sort_values(
            by=["Count", "rank_num", "rank_db"], ascending=[False, True, True],
            inplace=True)
        return sorted_count_one_candidate


# In[2]:


def sources_1(candidates_with_counts, merged_df, mer, sirius_df):
    """if only 1 source has confirmed the presence of a certain SMILES.
    This holds true when each candidate SMILES has only one source. The
    function selects the best candidate

    Parameters:
    candidates_with_counts: this is the result from the function add_count_column
    and contains a ordered dataframe, with the most sourced SMILES at top.
    merged_df: dataframe that contains all features from the input mzML file

    Returns:
    merged_df: with added top SMILES, Annotation Sources, Annotation Count, and
    MSI-Level

    Usage:
    sources_1(candidates_with_counts, merged_df)

    """

    df_count_1 = candidates_with_counts[candidates_with_counts["Count"] == 1]
    # rank number e.g:1
    #df_count_1["rank_num"] = [counts.split("_")[1] for counts in df_count_1["ranks"]]
    #df_count_1["rank_num"] = [int(x) for x in df_count_1["rank_num"]]

    # if only one unique source e.g: only SIRIUS gave out candidates
    #if len(np.unique(df_count_1["Source"])) == 1:

    #df_count_1 = df_count_1.sort_values(by="rank_num")

    df_count_1 = df_count_1[df_count_1["rank_num"] == min(df_count_1["rank_num"])]

    df_count_1["count_min"] = [
        str(df_count_1["SIRIUS"][x])
        + str(df_count_1["GNPS"][x])
        + str(df_count_1["MassBank"][x])
        + str(df_count_1["HMDB"][x])
        for x, row in df_count_1.iterrows()
    ]

    df_count_1["count_max"] = [x.count("_") for x in df_count_1["count_min"]]

    df_count_1 = df_count_1.sort_values(by="count_max", ascending=False)

    df_count_1.reset_index(drop=True, inplace=True)



#         #  choose the index/row from df_count_1 that contains priority wise databases
#         if "GNPS" in list(df_count_1["Source"]):
#             can1 = df_count_1.index[df_count_1["Source"].str.contains("GNPS")].tolist()
#             merged_df.loc[mer, "SMILES"] = list(df_count_1["SMILES"][can1])[0]
#             comp = pcp.get_compounds(merged_df["SMILES"][mer], 'smiles')
#             try:
#                 if comp:
#                     for c in comp:
#                         merged_df["synonyms"][mer] = c.synonyms
#                         merged_df.loc[mer, "IUPAC"] = c.iupac_name
#                         merged_df.loc[mer, "Formula"] = c.molecular_formula
#             except Exception:
#                 pass
#             merged_df["superclass"][mer] = np.nan
#             merged_df["class"][mer] = np.nan
#             merged_df["subclass"][mer] = np.nan
#             merged_df["ClassificationSource"][mer] = np.nan
#         elif "SIRIUS" in list(df_count_1["Source"]):
#             can1 = df_count_1.index[
#                 df_count_1["Source"].str.contains("SIRIUS")
#             ].tolist()
#             merged_df.loc[mer, "SMILES"] = list(df_count_1["SMILES"][can1])[0]
#             merged_df["Formula"] = sirius_df["molecularFormula"][0]
#             if "class" in sirius_df.columns:
#                 merged_df["superclass"] = sirius_df["superclass"][0]
#                 merged_df["class"] = sirius_df["class"][0]
#                 merged_df["subclass"] = sirius_df["subclass"][0]
#                 merged_df["ClassificationSource"] = "CANOPUS"

#         elif "MassBank" in list(df_count_1["Source"]):
#             can1 = df_count_1.index[
#                 df_count_1["Source"].str.contains("MassBank")
#             ].tolist()
#             merged_df.loc[mer, "SMILES"] = list(df_count_1["SMILES"][can1])[0]
#             comp = pcp.get_compounds(merged_df["SMILES"][mer], 'smiles')
#             try:
#                 if comp:
#                     for c in comp:
#                         merged_df["synonyms"][mer] = c.synonyms
#                         merged_df.loc[mer, "IUPAC"] = c.iupac_name
#                         merged_df.loc[mer, "Formula"] = c.molecular_formula
#             except Exception:
#                 pass
#             merged_df["superclass"][mer] = np.nan
#             merged_df["class"][mer] = np.nan
#             merged_df["subclass"][mer] = np.nan
#             merged_df["ClassificationSource"][mer] = np.nan
#         elif "HMDB" in list(df_count_1["Source"]):
#             can1 = df_count_1.index[df_count_1["Source"].str.contains("HMDB")].tolist()
#             merged_df.loc[mer, "SMILES"] = list(df_count_1["SMILES"][can1])[0]
#             comp = pcp.get_compounds(merged_df["SMILES"][mer], 'smiles')
#             try:
#                 if comp:
#                     for c in comp:
#                         merged_df["synonyms"][mer] = c.synonyms
#                         merged_df.loc[mer, "IUPAC"] = c.iupac_name
#                         merged_df.loc[mer, "Formula"] = c.molecular_formula
#             except Exception:
#                 pass
#             merged_df["superclass"][mer] = np.nan
#             merged_df["class"][mer] = np.nan
#             merged_df["subclass"][mer] = np.nan
#             merged_df["ClassificationSource"][mer] = np.nan





    merged_df.loc[mer, "AnnotationCount"] = df_count_1["Count"][0]

    gnps_indices = list(df_count_1[(df_count_1["GNPS"].notnull())].index)
    mbank_indices = list(df_count_1[(df_count_1["MassBank"].notnull())].index)
    hmdb_indices = list(df_count_1[(df_count_1["HMDB"].notnull())].index)
    sirius_indices = list(df_count_1[(df_count_1["SIRIUS"].notnull())].index)

    if 0 in sirius_indices:
        merged_df.loc[mer, "AnnotationSources"] = (
            str(merged_df["AnnotationSources"][mer]) + "|SIRIUS"
        )
    if 0 in mbank_indices:
        merged_df.loc[mer, "AnnotationSources"] = (
            str(merged_df["AnnotationSources"][mer]) + "|MassBank"
        )
        # print("mbank")
    if 0 in hmdb_indices:
        merged_df.loc[mer, "AnnotationSources"] = (
            str(merged_df["AnnotationSources"][mer]) + "|HMDB"
        )
        # print("hmdb")
    if 0 in gnps_indices:
        merged_df.loc[mer, "AnnotationSources"] = (
            str(merged_df["AnnotationSources"][mer]) + "|GNPS"
        )
        # print("gnps")

    

    merged_df["AnnotationSources"][mer] = merged_df["AnnotationSources"][
        mer
    ].replace("nan|", "")
    
    merged_df.loc[mer, "SMILES"] = df_count_1["SMILES"][0]
    comp = pcp.get_compounds(df_count_1["SMILES"][0], 'smiles')
    try:
        if comp:
            for c in comp:
                if c.cid:
                    merged_df["synonyms"][mer] = c.synonyms
                    merged_df.loc[mer, "IUPAC"] = c.iupac_name
                    merged_df.loc[mer, "Formula"] = c.molecular_formula
                    merged_df.loc[mer, "PubChemID"] = c.cid
                else:
                    merged_df.loc[mer, "Formula"] = sirius_df["molecularFormula"][0]
                    merged_df.loc[mer, "PubChemID"] = sirius_df["pubchemids"][0]
        
                
    except Exception:
        if "SIRIUS" == merged_df["AnnotationSources"][mer]:
            
            comp = pcp.get_compounds(sirius_df["pubchemids"][0], 'cid')
            try:
                if comp:
                    for c in comp:
                        if c.cid:
                            merged_df["synonyms"][mer] = c.synonyms
                            merged_df.loc[mer, "IUPAC"] = c.iupac_name
                            merged_df.loc[mer, "Formula"] = sirius_df["molecularFormula"][0]
                            merged_df.loc[mer, "PubChemID"] = sirius_df["pubchemids"][0]
                        else:
                        
                            merged_df.loc[mer, "Formula"] = sirius_df["molecularFormula"][0]
                            merged_df.loc[mer, "PubChemID"] = sirius_df["pubchemids"][0]
                else:
                    merged_df.loc[mer, "Formula"] = sirius_df["molecularFormula"][0]
                    merged_df["synonyms"][mer] = sirius_df["name"][0]
                    merged_df.loc[mer, "PubChemID"] = sirius_df["pubchemids"][0]
            except Exception:
                merged_df.loc[mer, "Formula"] = sirius_df["molecularFormula"][0]
                merged_df["synonyms"][mer] = sirius_df["name"][0]
                merged_df.loc[mer, "PubChemID"] = sirius_df["pubchemids"][0]
        else:
            pass
    if "SIRIUS" not in merged_df["AnnotationSources"][mer]:
        merged_df["superclass"][mer] = np.nan
        merged_df["class"][mer] = np.nan
        merged_df["subclass"][mer] = np.nan
        merged_df["ClassificationSource"][mer] = np.nan

    if (
        "HMDB" in merged_df["AnnotationSources"][mer]
        or "GNPS" in merged_df["AnnotationSources"][mer]
        or "MassBank" in merged_df["AnnotationSources"][mer]
    ):
        merged_df.loc[mer, "MSILevel"] = 2
    elif "SIRIUS" == merged_df["AnnotationSources"][mer]:
        merged_df.loc[mer, "MSILevel"] = 3

    

        
    # if only two or more unique source
#     else:

#         df_count_1 = candidates_with_counts[candidates_with_counts["Count"] == 1]
#         df_count_1["rank_num"] = [
#             counts.split("_")[1] for counts in df_count_1["ranks"]
#         ]
#         df_count_1["rank_num"] = [int(x) for x in df_count_1["rank_num"]]
#         df_count_1["count_min"] = [
#             str(df_count_1["SIRIUS"][x])
#             + str(df_count_1["GNPS"][x])
#             + str(df_count_1["MassBank"][x])
#             + str(df_count_1["HMDB"][x])
#             for x, row in df_count_1.iterrows()
#         ]
#         lengths = df_count_1["count_min"].str.len()
        
#         if len(np.unique(lengths)) > 1:
#             argmax = np.where(lengths == lengths.max())[0]
#             merged_df.loc[mer, "SMILES"] = df_count_1["SMILES"][argmax[0]]
#             merged_df.loc[mer, "AnnotationCount"] = df_count_1["Count"][argmax[0]]
#             df_count_na = df_count_1.notna()
            
#             list_sources = [*filter(df_count_na.loc[7].get, df_count_na.loc[7].index)]
            
#             if "GNPS" in list_sources:
#                 merged_df.loc[mer, "AnnotationSources"] = (
#                     str(merged_df["AnnotationSources"][mer]) + "|GNPS"
#                 )
#                 comp = pcp.get_compounds(merged_df["SMILES"][mer], 'smiles')
#                 try:
#                     if comp:
#                         for c in comp:
#                             merged_df["synonyms"][mer] = c.synonyms
#                             merged_df.loc[mer, "IUPAC"] = c.iupac_name
#                             merged_df.loc[mer, "Formula"] = c.molecular_formula
#                 except Exception:
#                     pass
#                 merged_df["superclass"][mer] = np.nan
#                 merged_df["class"][mer] = np.nan
#                 merged_df["subclass"][mer] = np.nan
#                 merged_df["ClassificationSource"][mer] = np.nan

            
#             elif "SIRIUS" in list_sources:
#                 merged_df.loc[mer, "AnnotationSources"] = (
#                     str(merged_df["AnnotationSources"][mer]) + "|SIRIUS"
#                 )
#                 merged_df["Formula"] = sirius_df["molecularFormula"][0]
#                 if "class" in sirius_df.columns:
#                     merged_df["superclass"] = sirius_df["superclass"][0]
#                     merged_df["class"] = sirius_df["class"][0]
#                     merged_df["subclass"] = sirius_df["subclass"][0]
#                     merged_df["ClassificationSource"] = "CANOPUS"
                
#             elif "MassBank" in list_sources:
#                 merged_df.loc[mer, "AnnotationSources"] = (
#                     str(merged_df["AnnotationSources"][mer]) + "|Massbank"
#                 )
#                 comp = pcp.get_compounds(merged_df["SMILES"][mer], 'smiles')
#                 try:
#                     if comp:
#                         for c in comp:
#                             merged_df["synonyms"][mer] = c.synonyms
#                             merged_df.loc[mer, "IUPAC"] = c.iupac_name
#                             merged_df.loc[mer, "Formula"] = c.molecular_formula
#                 except Exception:
#                     pass
#                 merged_df["superclass"][mer] = np.nan
#                 merged_df["class"][mer] = np.nan
#                 merged_df["subclass"][mer] = np.nan
#                 merged_df["ClassificationSource"][mer] = np.nan

#             elif "HMDB" in list_sources:
#                 merged_df.loc[mer, "AnnotationSources"] = (
#                     str(merged_df["AnnotationSources"][mer]) + "|HMDB"
#                 )
#                 comp = pcp.get_compounds(merged_df["SMILES"][mer], 'smiles')
#                 try:
#                     if comp:
#                         for c in comp:
#                             merged_df["synonyms"][mer] = c.synonyms
#                             merged_df.loc[mer, "IUPAC"] = c.iupac_name
#                             merged_df.loc[mer, "Formula"] = c.molecular_formula
#                 except Exception:
#                     pass
#                 merged_df["superclass"][mer] = np.nan
#                 merged_df["class"][mer] = np.nan
#                 merged_df["subclass"][mer] = np.nan
#                 merged_df["ClassificationSource"][mer] = np.nan

#             merged_df["AnnotationSources"][mer] = merged_df["AnnotationSources"][
#                 mer
#             ].replace("nan|", "")
            
            
#             if (
#                 "HMDB" in merged_df["AnnotationSources"][mer]
#                 or "GNPS" in merged_df["AnnotationSources"][mer]
#                 or "MassBank" in merged_df["AnnotationSources"][mer]
#             ):
#                 merged_df.loc[mer, "MSILevel"] = 2
                
#             if "SIRIUS" == merged_df["AnnotationSources"][mer]:
#                 merged_df.loc[mer, "MSILevel"] = 3
#         else:
            
#             if "GNPS" in np.unique(df_count_1["Source"]):
#                 index = np.where(df_count_1["Source"] == "GNPS")[0]
#                 merged_df.loc[mer, "SMILES"] = df_count_1["SMILES"][index[0]]
#                 merged_df.loc[mer, "AnnotationCount"] = df_count_1["Count"][index[0]]
#                 comp = pcp.get_compounds(df_count_1["SMILES"][index[0]], 'smiles')
#                 try:
#                     if comp:
#                         for c in comp:
#                             merged_df["synonyms"][mer] = c.synonyms
#                             merged_df.loc[mer, "IUPAC"] = c.iupac_name
#                             merged_df.loc[mer, "Formula"] = c.molecular_formula
#                 except Exception:
#                     pass
#                 merged_df["superclass"][mer] = np.nan
#                 merged_df["class"][mer] = np.nan
#                 merged_df["subclass"][mer] = np.nan
#                 merged_df["ClassificationSource"][mer] = np.nan

#                 merged_df.loc[mer, "AnnotationSources"] = (
#                     str(merged_df["AnnotationSources"][mer]) + "|GNPS"
#                 )
                
#             elif "SIRIUS" in np.unique(df_count_1["Source"]):
#                 index = np.where(df_count_1["Source"] == "SIRIUS")[0]
#                 merged_df.loc[mer, "SMILES"] = df_count_1["SMILES"][index[0]]
#                 merged_df.loc[mer, "AnnotationCount"] = df_count_1["Count"][index[0]]
#                 merged_df["Formula"] = sirius_df["molecularFormula"][0]
#                 if "class" in sirius_df.columns:
#                     merged_df["superclass"] = sirius_df["superclass"][0]
#                     merged_df["class"] = sirius_df["class"][0]
#                     merged_df["subclass"] = sirius_df["subclass"][0]
#                     merged_df["ClassificationSource"] = "CANOPUS"
#                 merged_df.loc[mer, "AnnotationSources"] = (
#                     str(merged_df["AnnotationSources"][mer]) + "|SIRIUS"
#                 )
#             elif "MassBank" in np.unique(df_count_1["Source"]):
#                 index = np.where(df_count_1["Source"] == "MassBank")[0]
#                 merged_df.loc[mer, "SMILES"] = df_count_1["SMILES"][index[0]]
#                 merged_df.loc[mer, "AnnotationCount"] = df_count_1["Count"][index[0]]
                
#                 comp = pcp.get_compounds(df_count_1["SMILES"][index[0]], 'smiles')
#                 try:
#                     if comp:
#                         for c in comp:
#                             merged_df["synonyms"][mer] = c.synonyms
#                             merged_df.loc[mer, "IUPAC"] = c.iupac_name
#                             merged_df.loc[mer, "Formula"] = c.molecular_formula
#                 except Exception:
#                     pass
#                 merged_df["superclass"][mer] = np.nan
#                 merged_df["class"][mer] = np.nan
#                 merged_df["subclass"][mer] = np.nan
#                 merged_df["ClassificationSource"][mer] = np.nan
#                 merged_df.loc[mer, "AnnotationSources"] = (
#                     str(merged_df["AnnotationSources"][mer]) + "|MassBank"
#                 )

#             elif "HMDB" in np.unique(df_count_1["Source"]):
#                 index = np.where(df_count_1["Source"] == "HMDB")[0]
#                 merged_df.loc[mer, "SMILES"] = df_count_1["SMILES"][index[0]]
#                 merged_df.loc[mer, "AnnotationCount"] = df_count_1["Count"][index[0]]
                
#                 comp = pcp.get_compounds(df_count_1["SMILES"][index[0]], 'smiles')
#                 try:
#                     if comp:
#                         for c in comp:
#                             merged_df["synonyms"][mer] = c.synonyms
#                             merged_df.loc[mer, "IUPAC"] = c.iupac_name
#                             merged_df.loc[mer, "Formula"] = c.molecular_formula
#                 except Exception:
#                     pass
#                 merged_df["superclass"][mer] = np.nan
#                 merged_df["class"][mer] = np.nan
#                 merged_df["subclass"][mer] = np.nan
#                 merged_df["ClassificationSource"][mer] = np.nan
#                 merged_df.loc[mer, "AnnotationSources"] = (
#                     str(merged_df["AnnotationSources"][mer]) + "|HMDB"
#                 )

            
#             merged_df["AnnotationSources"][mer] = merged_df["AnnotationSources"][
#                 mer
#             ].replace("nan|", "")
#             if "SIRIUS" == merged_df["AnnotationSources"][mer]:
#                 merged_df.loc[mer, "MSILevel"] = 3

#             if (
#                 "HMDB" in merged_df["AnnotationSources"][mer]
#                 or "GNPS" in merged_df["AnnotationSources"][mer]
#                 or "MassBank" in merged_df["AnnotationSources"][mer]
#             ):
#                 merged_df.loc[mer, "MSILevel"] = 2
#                 merged_df["superclass"][mer] = np.nan
#                 merged_df["class"][mer] = np.nan
#                 merged_df["subclass"][mer] = np.nan
#                 merged_df["ClassificationSource"][mer] = np.nan
            
    return merged_df


# In[1]:


def sources_2(candidates_with_counts, merged_df, mer, sirius_df):

    """if only 2 sources have confirmed the presence of a certain SMILES.
    This holds true when each candidate SMILES has only two sources. The
    function selects the best candidate and adds the two sources as
    annotation sources

    Parameters:
    candidates_with_counts: this is the result from the function add_count_column
    and contains a ordered dataframe, with the most sourced SMILES at top.
    merged_df: dataframe that contains all features from the input mzML file

    Returns:
    merged_df: with added top SMILES, Annotation Sources, Annotation Count, and
    MSI-Level

    Usage:
    sources_2(candidates_with_counts, merged_df, mer)

    """

    df_count_2 = candidates_with_counts[candidates_with_counts["Count"] == 2]
    #df_count_2["rank_num"] = [counts.split("_")[1] for counts in df_count_2["ranks"]]
    #df_count_2["rank_num"] = [int(x) for x in df_count_2["rank_num"]]
    #df_count_2 = df_count_2.sort_values(by="rank_num")
    df_countnew = df_count_2[df_count_2["rank_num"] == min(df_count_2["rank_num"])]
    df_countnew["count_min"] = [
        str(df_countnew["SIRIUS"][x])
        + str(df_countnew["GNPS"][x])
        + str(df_countnew["MassBank"][x])
        + str(df_countnew["HMDB"][x])
        for x, row in df_countnew.iterrows()
    ]
    df_countnew["count_max"] = [x.count("_") for x in df_countnew["count_min"]]
    df_countnew = df_countnew.sort_values(by="count_max", ascending=False)

    df_countnew.reset_index(drop=True, inplace=True)

    gnps_indices = list(df_countnew[(df_countnew["GNPS"].notnull())].index)
    mbank_indices = list(df_countnew[(df_countnew["MassBank"].notnull())].index)
    hmdb_indices = list(df_countnew[(df_countnew["HMDB"].notnull())].index)
    sirius_indices = list(df_countnew[(df_countnew["SIRIUS"].notnull())].index)

    if 0 in sirius_indices:
        # print("sirius")
        merged_df.loc[mer, "AnnotationSources"] = (
            str(merged_df["AnnotationSources"][mer]) + "|SIRIUS"
        )
    if 0 in mbank_indices:
        merged_df.loc[mer, "AnnotationSources"] = (
            str(merged_df["AnnotationSources"][mer]) + "|MassBank"
        )
        # print("mbank")
    if 0 in hmdb_indices:
        merged_df.loc[mer, "AnnotationSources"] = (
            str(merged_df["AnnotationSources"][mer]) + "|HMDB"
        )
        # print("hmdb")
    if 0 in gnps_indices:
        merged_df.loc[mer, "AnnotationSources"] = (
            str(merged_df["AnnotationSources"][mer]) + "|GNPS"
        )
        # print("gnps")
#     if (
#         "HMDB" in merged_df["AnnotationSources"][mer]
#         or "GNPS" in merged_df["AnnotationSources"][mer]
#         or "MassBank" in merged_df["AnnotationSources"][mer]
#     ):
    merged_df.loc[mer, "MSILevel"] = 2
#     if "nan|SIRIUS" == merged_df["AnnotationSources"][mer]:
#         merged_df.loc[mer, "MSILevel"] = 3
        

    merged_df["AnnotationSources"][mer] = merged_df["AnnotationSources"][mer].replace(
        "nan|", ""
    )
    merged_df.loc[mer, "AnnotationCount"] = df_countnew["Count"][0]
    
    merged_df.loc[mer, "SMILES"] = df_countnew["SMILES"][0]
    comp = pcp.get_compounds(df_countnew["SMILES"][0], 'smiles')
    try:
        if comp:
            for c in comp:
                merged_df["synonyms"][mer] = c.synonyms
                merged_df.loc[mer, "IUPAC"] = c.iupac_name
                merged_df.loc[mer, "Formula"] = c.molecular_formula
                merged_df.loc[mer, "PubChemID"] = c.cid
    except Exception:
        pass
    if "SIRIUS" not in merged_df["AnnotationSources"][mer]:
        merged_df["superclass"][mer] = np.nan
        merged_df["class"][mer] = np.nan
        merged_df["subclass"][mer] = np.nan
        merged_df["ClassificationSource"][mer] = np.nan

#     if "GNPS" in merged_df["AnnotationSources"][mer]:
#         if len(df_count_2.loc[df_count_2["Source"] == "GNPS"]) == 1:
#             new = df_count_2.loc[df_count_2["Source"] == "GNPS"]
#             new.reset_index(drop=True, inplace=True)
#             merged_df.loc[mer, "SMILES"] = new["SMILES"][0]
#             comp = pcp.get_compounds(new["SMILES"][0], 'smiles')
#             try:
#                 if comp:
#                     for c in comp:
#                         merged_df["synonyms"][mer] = c.synonyms
#                         merged_df.loc[mer, "IUPAC"] = c.iupac_name
#                         merged_df.loc[mer, "Formula"] = c.molecular_formula
#             except Exception:
#                 pass
#             merged_df["superclass"][mer] = np.nan
#             merged_df["class"][mer] = np.nan
#             merged_df["subclass"][mer] = np.nan
#             merged_df["ClassificationSource"][mer] = np.nan
#         elif len(df_count_2.loc[df_count_2["Source"] == "GNPS"]) > 1:
#             new = df_count_2.loc[df_count_2["Source"] == "GNPS"]
#             new.reset_index(drop=True, inplace=True)
#             merged_df.loc[mer, "SMILES"] = new["SMILES"][0]
#             comp = pcp.get_compounds(new["SMILES"][0], 'smiles')
#             try:
#                 if comp:
#                     for c in comp:
#                         merged_df["synonyms"][mer] = c.synonyms
#                         merged_df.loc[mer, "IUPAC"] = c.iupac_name
#                         merged_df.loc[mer, "Formula"] = c.molecular_formula
#             except Exception:
#                 pass
#             merged_df["superclass"][mer] = np.nan
#             merged_df["class"][mer] = np.nan
#             merged_df["subclass"][mer] = np.nan
#             merged_df["ClassificationSource"][mer] = np.nan
#         else:
#             pass
        
        
        
#     elif "SIRIUS" in merged_df["AnnotationSources"][mer]:
#         if len(df_count_2.loc[df_count_2["Source"] == "SIRIUS"]) == 1:
#             new = df_count_2.loc[df_count_2["Source"] == "SIRIUS"]
#             new.reset_index(drop=True, inplace=True)
#             merged_df.loc[mer, "SMILES"] = new["SMILES"][0]
#             merged_df["Formula"] = sirius_df["molecularFormula"][0]
#             if "class" in sirius_df.columns:
#                 merged_df["superclass"] = sirius_df["superclass"][0]
#                 merged_df["class"] = sirius_df["class"][0]
#                 merged_df["subclass"] = sirius_df["subclass"][0]
#                 merged_df["ClassificationSource"] = "CANOPUS"
#         elif len(df_count_2.loc[df_count_2["Source"] == "SIRIUS"]) > 1:
#             new = df_count_2.loc[df_count_2["Source"] == "SIRIUS"]
#             new.reset_index(drop=True, inplace=True)
#             merged_df.loc[mer, "SMILES"] = new["SMILES"][0]
#             merged_df["Formula"] = sirius_df["molecularFormula"][0]
#             if "class" in sirius_df.columns:
#                 merged_df["superclass"] = sirius_df["superclass"][0]
#                 merged_df["class"] = sirius_df["class"][0]
#                 merged_df["subclass"] = sirius_df["subclass"][0]
#                 merged_df["ClassificationSource"] = "CANOPUS"
            
#         else:
#             pass
#     elif "MassBank" in merged_df["AnnotationSources"][mer]:
#         if len(df_count_2.loc[df_count_2["Source"] == "MassBank"]) == 1:
#             new = df_count_2.loc[df_count_2["Source"] == "MassBank"]
#             new.reset_index(drop=True, inplace=True)
#             merged_df.loc[mer, "SMILES"] = new["SMILES"][0]
#             comp = pcp.get_compounds(merged_df["SMILES"][mer], 'smiles')
#             try:
#                 if comp:
#                     for c in comp:
#                         merged_df["synonyms"][mer] = c.synonyms
#                         merged_df.loc[mer, "IUPAC"] = c.iupac_name
#                         merged_df.loc[mer, "Formula"] = c.molecular_formula
#             except Exception:
#                 pass
#             merged_df["superclass"][mer] = np.nan
#             merged_df["class"][mer] = np.nan
#             merged_df["subclass"][mer] = np.nan
#             merged_df["ClassificationSource"][mer] = np.nan
#         elif len(df_count_2.loc[df_count_2["Source"] == "MassBank"]) > 1:
#             new = df_count_2.loc[df_count_2["Source"] == "MassBank"]
#             new.reset_index(drop=True, inplace=True)
#             merged_df.loc[mer, "SMILES"] = new["SMILES"][0]
#             comp = pcp.get_compounds(new["SMILES"][0], 'smiles')
#             try:
#                 if comp:
#                     for c in comp:
#                         merged_df["synonyms"][mer] = c.synonyms
#                         merged_df.loc[mer, "IUPAC"] = c.iupac_name
#                         merged_df.loc[mer, "Formula"] = c.molecular_formula
#             except Exception:
#                 pass
#             merged_df["superclass"][mer] = np.nan
#             merged_df["class"][mer] = np.nan
#             merged_df["subclass"][mer] = np.nan
#             merged_df["ClassificationSource"][mer] = np.nan
#         else:
#             pass
#     elif "HMDB" in merged_df["AnnotationSources"][mer]:
#         if len(df_count_2.loc[df_count_2["Source"] == "HMDB"]) == 1:
#             new = df_count_2.loc[df_count_2["Source"] == "HMDB"]
#             new.reset_index(drop=True, inplace=True)
#             merged_df.loc[mer, "SMILES"] = new["SMILES"][0]
#             comp = pcp.get_compounds(new["SMILES"][0], 'smiles')
#             try:
#                 if comp:
#                     for c in comp:
#                         merged_df["synonyms"][mer] = c.synonyms
#                         merged_df.loc[mer, "IUPAC"] = c.iupac_name
#                         merged_df.loc[mer, "Formula"] = c.molecular_formula
#             except Exception:
#                 pass
#             merged_df["superclass"] = np.nan
#             merged_df["class"] = np.nan
#             merged_df["subclass"] = np.nan
#             merged_df["ClassificationSource"] = np.nan
#         elif len(df_count_2.loc[df_count_2["Source"] == "HMDB"]) > 1:
#             new = df_count_2.loc[df_count_2["Source"] == "HMDB"]
#             new.reset_index(drop=True, inplace=True)
#             merged_df.loc[mer, "SMILES"] = new["SMILES"][0]
#             comp = pcp.get_compounds(new["SMILES"][0], 'smiles')
#             try:
#                 if comp:
#                     for c in comp:
#                         merged_df["synonyms"][mer] = c.synonyms
#                         merged_df.loc[mer, "IUPAC"] = c.iupac_name
#                         merged_df.loc[mer, "Formula"] = c.molecular_formula
#             except Exception:
#                 pass
#             merged_df["superclass"][mer] = np.nan
#             merged_df["class"][mer] = np.nan
#             merged_df["subclass"][mer] = np.nan
#             merged_df["ClassificationSource"][mer] = np.nan
#         else:
#             pass
    
    
    return merged_df


# In[29]:


def sources_3(candidates_with_counts, merged_df, mer, sirius_df):

    """if only 3 sources have confirmed the presence of a certain SMILES.
    This holds true when each candidate SMILES has only 3 sources. The
    function selects the best candidate and adds the 3 sources as
    annotation sources

    Parameters:
    candidates_with_counts: this is the result from the function add_count_column
    and contains a ordered dataframe, with the most sourced SMILES at top.
    merged_df: dataframe that contains all features from the input mzML file

    Returns:
    merged_df: with added top SMILES, Annotation Sources, Annotation Count, and
    MSI-Level

    Usage:
    sources_2(candidates_with_counts, merged_df, mer)

    """
    # if the count is 3
    df_count_3 = candidates_with_counts[candidates_with_counts["Count"] == 3]
    # extracts the ranks again
    #df_count_3["rank_num"] = [counts.split("_")[1] for counts in df_count_3["ranks"]]
    #df_count_3["rank_num"] = [int(x) for x in df_count_3["rank_num"]]

    #df_count_3 = df_count_3.sort_values(by="rank_num")
    df_count_3 = df_count_3[df_count_3["rank_num"] == min(df_count_3["rank_num"])]
    df_count_3["count_min"] = [
        str(df_count_3["SIRIUS"][x])
        + str(df_count_3["GNPS"][x])
        + str(df_count_3["MassBank"][x])
        + str(df_count_3["HMDB"][x])
        for x, row in df_count_3.iterrows()
    ]

    df_count_3["count_max"] = [x.count("_") for x in df_count_3["count_min"]]
    df_count_3 = df_count_3.sort_values(by="count_max", ascending=False)

    df_count_3.reset_index(drop=True, inplace=True)

    merged_df.loc[mer, "AnnotationCount"] = df_count_3["Count"][0]
    
    heavy_atoms = ["C", "N", "P", "O", "S"]

    Mol = []

    for j in list(df_count_3["SMILES"]):
        if not isNaN(j):
            # print(type(j))
            mol2 = Chem.MolFromSmiles(j)
            Mol.append(mol2)

    if len(Mol) >= 2:
        res = rdFMCS.FindMCS(Mol, timeout=60)
        sm_res = Chem.MolToSmiles(Chem.MolFromSmarts(res.smartsString))
        # if there are atleast 3 heavy atoms in the MCSS, then add it to the result file
        elem = [ele for ele in heavy_atoms if (ele in sm_res)]
        if elem and len(sm_res) >= 3:
            merged_df.loc[mer, "MCSS"] = Chem.MolToSmiles(
                Chem.MolFromSmarts(res.smartsString)
            )

    gnps_indices = list(df_count_3[(df_count_3["GNPS"].notnull())].index)
    mbank_indices = list(df_count_3[(df_count_3["MassBank"].notnull())].index)
    hmdb_indices = list(df_count_3[(df_count_3["HMDB"].notnull())].index)
    sirius_indices = list(df_count_3[(df_count_3["SIRIUS"].notnull())].index)

    if 0 in sirius_indices:
        # print("sirius")
        merged_df.loc[mer, "AnnotationSources"] = (
            str(merged_df["AnnotationSources"][mer]) + "|SIRIUS"
        )
    if 0 in mbank_indices:
        merged_df.loc[mer, "AnnotationSources"] = (
            str(merged_df["AnnotationSources"][mer]) + "|MassBank"
        )
        # print("mbank")
    if 0 in hmdb_indices:
        merged_df.loc[mer, "AnnotationSources"] = (
            str(merged_df["AnnotationSources"][mer]) + "|HMDB"
        )
        # print("hmdb")
    if 0 in gnps_indices:
        merged_df.loc[mer, "AnnotationSources"] = (
            str(merged_df["AnnotationSources"][mer]) + "|GNPS"
        )
        # print("gnps")
    if "nan|SIRIUS" == merged_df["AnnotationSources"][mer]:
        merged_df.loc[mer, "MSILevel"] = 3

    merged_df["AnnotationSources"][mer] = merged_df["AnnotationSources"][mer].replace(
        "nan|", ""
    )
#     if (
#         "HMDB" in merged_df["AnnotationSources"][mer]
#         or "GNPS" in merged_df["AnnotationSources"][mer]
#         or "MassBank" in merged_df["AnnotationSources"][mer]
#     ):
    merged_df.loc[mer, "MSILevel"] = 2
    merged_df.loc[mer, "SMILES"] = df_count_3["SMILES"][0]
    comp = pcp.get_compounds(df_count_3["SMILES"][0], 'smiles')
    try:
        if comp:
            for c in comp:
                merged_df["synonyms"][mer] = c.synonyms
                merged_df.loc[mer, "IUPAC"] = c.iupac_name
                merged_df.loc[mer, "Formula"] = c.molecular_formula
                merged_df.loc[mer, "PubChemID"] = c.cid
    except Exception:
        pass
    if "SIRIUS" not in merged_df["AnnotationSources"][mer]:
        merged_df["superclass"][mer] = np.nan
        merged_df["class"][mer] = np.nan
        merged_df["subclass"][mer] = np.nan
        merged_df["ClassificationSource"][mer] = np.nan
    
    

#     if "GNPS" in merged_df["AnnotationSources"][mer]:
#         if len(df_count_3.loc[df_count_3["Source"] == "GNPS"]) == 1:
#             new = df_count_3.loc[df_count_3["Source"] == "GNPS"]
#             new.reset_index(drop=True, inplace=True)
#             merged_df.loc[mer, "SMILES"] = new["SMILES"][0]
#             comp = pcp.get_compounds(new["SMILES"][0], 'smiles')
#             try:
#                 if comp:
#                     for c in comp:
#                         merged_df["synonyms"][mer] = c.synonyms
#                         merged_df.loc[mer, "IUPAC"] = c.iupac_name
#                         merged_df.loc[mer, "Formula"] = c.molecular_formula
#             except Exception:
#                 pass
#             merged_df["superclass"][mer] = np.nan
#             merged_df["class"][mer] = np.nan
#             merged_df["subclass"][mer] = np.nan
#             merged_df["ClassificationSource"][mer] = np.nan
#         elif len(df_count_3.loc[df_count_3["Source"] == "GNPS"]) > 1:
#             new = df_count_3.loc[df_count_3["Source"] == "GNPS"]
#             new.reset_index(drop=True, inplace=True)
#             merged_df.loc[mer, "SMILES"] = new["SMILES"][0]
#             comp = pcp.get_compounds(new["SMILES"][0], 'smiles')
#             try:
#                 if comp:
#                     for c in comp:
#                         merged_df["synonyms"][mer] = c.synonyms
#                         merged_df.loc[mer, "IUPAC"] = c.iupac_name
#                         merged_df.loc[mer, "Formula"] = c.molecular_formula
#             except Exception:
#                 pass
#             merged_df["superclass"][mer] = np.nan
#             merged_df["class"][mer] = np.nan
#             merged_df["subclass"][mer] = np.nan
#             merged_df["ClassificationSource"][mer] = np.nan
#         else:
#             pass
    
    
#     elif "SIRIUS" in merged_df["AnnotationSources"][mer]:
#         if len(df_count_3.loc[df_count_3["Source"] == "SIRIUS"]) == 1:
#             new = df_count_3.loc[df_count_3["Source"] == "SIRIUS"]
#             new.reset_index(drop=True, inplace=True)
#             merged_df.loc[mer, "SMILES"] = new["SMILES"][0]
#             merged_df["Formula"] = sirius_df["molecularFormula"][0]
#             if "class" in sirius_df.columns:
#                 merged_df["superclass"] = sirius_df["superclass"][0]
#                 merged_df["class"] = sirius_df["class"][0]
#                 merged_df["subclass"] = sirius_df["subclass"][0]
#                 merged_df["ClassificationSource"] = "CANOPUS"
#         elif len(df_count_3.loc[df_count_3["Source"] == "SIRIUS"]) > 1:
#             new = df_count_3.loc[df_count_3["Source"] == "SIRIUS"]
#             new.reset_index(drop=True, inplace=True)
#             merged_df.loc[mer, "SMILES"] = new["SMILES"][0]
#             merged_df["Formula"] = sirius_df["molecularFormula"][0]
#             if "class" in sirius_df.columns:
#                 merged_df["superclass"] = sirius_df["superclass"][0]
#                 merged_df["class"] = sirius_df["class"][0]
#                 merged_df["subclass"] = sirius_df["subclass"][0]
#                 merged_df["ClassificationSource"] = "CANOPUS"
#         else:
#             pass
    
#     elif "MassBank" in merged_df["AnnotationSources"][mer]:
#         if len(df_count_3.loc[df_count_3["Source"] == "MassBank"]) == 1:
#             new = df_count_3.loc[df_count_3["Source"] == "MassBank"]
#             new.reset_index(drop=True, inplace=True)
#             merged_df.loc[mer, "SMILES"] = new["SMILES"][0]
#             comp = pcp.get_compounds(new["SMILES"][0], 'smiles')
#             try:
#                 if comp:
#                     for c in comp:
#                         merged_df["synonyms"][mer] = c.synonyms
#                         merged_df.loc[mer, "IUPAC"] = c.iupac_name
#                         merged_df.loc[mer, "Formula"] = c.molecular_formula
#             except Exception:
#                 pass
#             merged_df["superclass"][mer] = np.nan
#             merged_df["class"][mer] = np.nan
#             merged_df["subclass"][mer] = np.nan
#             merged_df["ClassificationSource"][mer] = np.nan
#         elif len(df_count_3.loc[df_count_3["Source"] == "MassBank"]) > 1:
#             new = df_count_3.loc[df_count_3["Source"] == "MassBank"]
#             new.reset_index(drop=True, inplace=True)
#             merged_df.loc[mer, "SMILES"] = new["SMILES"][0]
#             comp = pcp.get_compounds(new["SMILES"][0], 'smiles')
#             try:
#                 if comp:
#                     for c in comp:
#                         merged_df["synonyms"][mer] = c.synonyms
#                         merged_df.loc[mer, "IUPAC"] = c.iupac_name
#                         merged_df.loc[mer, "Formula"] = c.molecular_formula
#             except Exception:
#                 pass
#             merged_df["superclass"][mer] = np.nan
#             merged_df["class"][mer] = np.nan
#             merged_df["subclass"][mer] = np.nan
#             merged_df["ClassificationSource"][mer] = np.nan
#         else:
#             pass
#     elif "HMDB" in merged_df["AnnotationSources"][mer]:
#         if len(df_count_3.loc[df_count_3["Source"] == "HMDB"]) == 1:
#             new = df_count_3.loc[df_count_3["Source"] == "HMDB"]
#             new.reset_index(drop=True, inplace=True)
#             merged_df.loc[mer, "SMILES"] = new["SMILES"][0]
#             comp = pcp.get_compounds(new["SMILES"][0], 'smiles')
#             try:
#                 if comp:
#                     for c in comp:
#                         merged_df["synonyms"][mer] = c.synonyms
#                         merged_df.loc[mer, "IUPAC"] = c.iupac_name
#                         merged_df.loc[mer, "Formula"] = c.molecular_formula
#             except Exception:
#                 pass
#             merged_df["superclass"][mer] = np.nan
#             merged_df["class"][mer] = np.nan
#             merged_df["subclass"][mer] = np.nan
#             merged_df["ClassificationSource"][mer] = np.nan
#         elif len(df_count_3.loc[df_count_3["Source"] == "HMDB"]) > 1:
#             new = df_count_3.loc[df_count_3["Source"] == "HMDB"]
#             new.reset_index(drop=True, inplace=True)
#             merged_df.loc[mer, "SMILES"] = new["SMILES"][0]
#             comp = pcp.get_compounds(new["SMILES"][0], 'smiles')
#             try:
#                 if comp:
#                     for c in comp:
#                         merged_df["synonyms"][mer] = c.synonyms
#                         merged_df.loc[mer, "IUPAC"] = c.iupac_name
#                         merged_df.loc[mer, "Formula"] = c.molecular_formula
#             except Exception:
#                 pass
#             merged_df["superclass"][mer] = np.nan
#             merged_df["class"][mer] = np.nan
#             merged_df["subclass"][mer] = np.nan
#             merged_df["ClassificationSource"][mer] = np.nan
#         else:
#             pass
    
    
    
    return merged_df


# In[30]:


def sources_4(candidates_with_counts, merged_df, mer, sirius_df):

    """if only 3 sources have confirmed the presence of a certain SMILES.
    This holds true when each candidate SMILES has only 3 sources. The
    function selects the best candidate and adds the 3 sources as
    annotation sources

    Parameters:
    candidates_with_counts: this is the result from the function add_count_column
    and contains a ordered dataframe, with the most sourced SMILES at top.
    merged_df: dataframe that contains all features from the input mzML file

    Returns:
    merged_df: with added top SMILES, Annotation Sources, Annotation Count, and
    MSI-Level

    Usage:
    sources_2(candidates_with_counts, merged_df, mer)

    """

    df_count_4 = candidates_with_counts[candidates_with_counts["Count"] == 4]
    #df_count_4["rank_num"] = [counts.split("_")[1] for counts in df_count_4["ranks"]]
    #df_count_4["rank_num"] = [int(x) for x in df_count_4["rank_num"]]
    #df_count_4 = df_count_4.sort_values(by="rank_num")
    df_count_4 = df_count_4[df_count_4["rank_num"] == min(df_count_4["rank_num"])]
    df_count_4["count_min"] = [
        str(df_count_4["SIRIUS"][x])
        + str(df_count_4["GNPS"][x])
        + str(df_count_4["MassBank"][x])
        + str(df_count_4["HMDB"][x])
        for x, row in df_count_4.iterrows()
    ]
    df_count_4["count_max"] = [x.count("_") for x in df_count_4["count_min"]]
    df_count_4 = df_count_4.sort_values(by="count_max", ascending=False)

    df_count_4.reset_index(drop=True, inplace=True)

    merged_df.loc[mer, "AnnotationCount"] = df_count_4["Count"][0]
    
    heavy_atoms = ["C", "N", "P", "O", "S"]
    
    Mol = []

    for j in list(df_count_4["SMILES"]):
        if not isNaN(j):
            # print(type(j))
            mol2 = Chem.MolFromSmiles(j)
            Mol.append(mol2)

    if len(Mol) >= 2:
        res = rdFMCS.FindMCS(Mol, timeout=60)
        sm_res = Chem.MolToSmiles(Chem.MolFromSmarts(res.smartsString))
        # if there are atleast 3 heavy atoms in the MCSS, then add it to the result file
        elem = [ele for ele in heavy_atoms if (ele in sm_res)]
        if elem and len(sm_res) >= 3:
            merged_df.loc[mer, "MCSS"] = Chem.MolToSmiles(
                Chem.MolFromSmarts(res.smartsString)
            )

    gnps_indices = list(df_count_4[(df_count_4["GNPS"].notnull())].index)
    mbank_indices = list(df_count_4[(df_count_4["MassBank"].notnull())].index)
    hmdb_indices = list(df_count_4[(df_count_4["HMDB"].notnull())].index)
    sirius_indices = list(df_count_4[(df_count_4["SIRIUS"].notnull())].index)

    if 0 in sirius_indices:
        # print("sirius")
        merged_df.loc[mer, "AnnotationSources"] = (
            str(merged_df["AnnotationSources"][mer]) + "|SIRIUS"
        )
    if 0 in mbank_indices:
        merged_df.loc[mer, "AnnotationSources"] = (
            str(merged_df["AnnotationSources"][mer]) + "|Massbank"
        )
        # print("mbank")
    if 0 in hmdb_indices:
        merged_df.loc[mer, "AnnotationSources"] = (
            str(merged_df["AnnotationSources"][mer]) + "|HMDB"
        )
        # print("hmdb")
    if 0 in gnps_indices:
        merged_df.loc[mer, "AnnotationSources"] = (
            str(merged_df["AnnotationSources"][mer]) + "|GNPS"
        )
        # print("gnps")
    
    if "nan|SIRIUS" == merged_df["AnnotationSources"][mer]:
        merged_df.loc[mer, "MSILevel"] = 3

    merged_df["AnnotationSources"][mer] = merged_df["AnnotationSources"][mer].replace(
        "nan|", ""
    )
    
#     if (
#         "HMDB" in merged_df["AnnotationSources"][mer]
#         or "GNPS" in merged_df["AnnotationSources"][mer]
#         or "MassBank" in merged_df["AnnotationSources"][mer]
#     ):
    merged_df.loc[mer, "MSILevel"] = 2
    
    merged_df.loc[mer, "SMILES"] = df_count_4["SMILES"][0]
    comp = pcp.get_compounds(df_count_4["SMILES"][0], 'smiles')
    try:
        if comp:
            for c in comp:
                merged_df["synonyms"][mer] = c.synonyms
                merged_df.loc[mer, "IUPAC"] = c.iupac_name
                merged_df.loc[mer, "Formula"] = c.molecular_formula
                merged_df.loc[mer, "PubChemID"] = c.cid
    except Exception:
        pass
    if "SIRIUS" not in merged_df["AnnotationSources"][mer]:
        merged_df["superclass"][mer] = np.nan
        merged_df["class"][mer] = np.nan
        merged_df["subclass"][mer] = np.nan
        merged_df["ClassificationSource"][mer] = np.nan
    

    
#     if "GNPS" in merged_df["AnnotationSources"][mer]:
#         if len(df_count_4.loc[df_count_4["Source"] == "GNPS"]) == 1:
#             new = df_count_4.loc[df_count_4["Source"] == "GNPS"]
#             new.reset_index(drop=True, inplace=True)
#             merged_df.loc[mer, "SMILES"] = new["SMILES"][0]
#             comp = pcp.get_compounds(new["SMILES"][0], 'smiles')
#             try:
#                 if comp:
#                     for c in comp:
#                         merged_df["synonyms"][mer] = c.synonyms
#                         merged_df.loc[mer, "IUPAC"] = c.iupac_name
#                         merged_df.loc[mer, "Formula"] = c.molecular_formula
#             except Exception:
#                 pass
#             merged_df["superclass"][mer] = np.nan
#             merged_df["class"][mer] = np.nan
#             merged_df["subclass"][mer] = np.nan
#             merged_df["ClassificationSource"][mer] = np.nan
#         elif len(df_count_4.loc[df_count_4["Source"] == "GNPS"]) > 1:
#             new = df_count_4.loc[df_count_4["Source"] == "GNPS"]
#             new.reset_index(drop=True, inplace=True)
#             merged_df.loc[mer, "SMILES"] = new["SMILES"][0]
#             comp = pcp.get_compounds(new["SMILES"][0], 'smiles')
#             try:
#                 if comp:
#                     for c in comp:
#                         merged_df["synonyms"][mer] = c.synonyms
#                         merged_df.loc[mer, "IUPAC"] = c.iupac_name
#                         merged_df.loc[mer, "Formula"] = c.molecular_formula
#             except Exception:
#                 pass
#             merged_df["superclass"][mer] = np.nan
#             merged_df["class"][mer] = np.nan
#             merged_df["subclass"][mer] = np.nan
#             merged_df["ClassificationSource"][mer] = np.nan
#         else:
#             pass
    
#     elif "SIRIUS" in merged_df["AnnotationSources"][mer]:
#         if len(df_count_4.loc[df_count_4["Source"] == "SIRIUS"]) == 1:
#             new = df_count_4.loc[df_count_4["Source"] == "SIRIUS"]
#             new.reset_index(drop=True, inplace=True)
#             merged_df.loc[mer, "SMILES"] = new["SMILES"][0]
#             merged_df["Formula"] = sirius_df["molecularFormula"][0]
#             if "class" in sirius_df.columns:
#                 merged_df["superclass"] = sirius_df["superclass"][0]
#                 merged_df["class"] = sirius_df["class"][0]
#                 merged_df["subclass"] = sirius_df["subclass"][0]
#                 merged_df["ClassificationSource"] = "CANOPUS"
#         elif len(df_count_4.loc[df_count_4["Source"] == "SIRIUS"]) > 1:
#             new = df_count_4.loc[df_count_4["Source"] == "SIRIUS"]
#             new.reset_index(drop=True, inplace=True)
#             merged_df.loc[mer, "SMILES"] = new["SMILES"][0]
#             merged_df["Formula"] = sirius_df["molecularFormula"][0]
#             if "class" in sirius_df.columns:
#                 merged_df["superclass"] = sirius_df["superclass"][0]
#                 merged_df["class"] = sirius_df["class"][0]
#                 merged_df["subclass"] = sirius_df["subclass"][0]
#                 merged_df["ClassificationSource"] = "CANOPUS"
#         else:
#             pass
    
    
#     elif "MassBank" in merged_df["AnnotationSources"][mer]:
#         if len(df_count_4.loc[df_count_4["Source"] == "MassBank"]) == 1:
#             new = df_count_4.loc[df_count_4["Source"] == "MassBank"]
#             new.reset_index(drop=True, inplace=True)
#             merged_df.loc[mer, "SMILES"] = new["SMILES"][0]
#             comp = pcp.get_compounds(new["SMILES"][0], 'smiles')
#             try:
#                 if comp:
#                     for c in comp:
#                         merged_df["synonyms"][mer] = c.synonyms
#                         merged_df.loc[mer, "IUPAC"] = c.iupac_name
#                         merged_df.loc[mer, "Formula"] = c.molecular_formula
#             except Exception:
#                 pass
#             merged_df["superclass"][mer] = np.nan
#             merged_df["class"][mer] = np.nan
#             merged_df["subclass"][mer] = np.nan
#             merged_df["ClassificationSource"][mer] = np.nan
#         elif len(df_count_4.loc[df_count_4["Source"] == "MassBank"]) > 1:
#             new = df_count_4.loc[df_count_4["Source"] == "MassBank"]
#             new.reset_index(drop=True, inplace=True)
#             merged_df.loc[mer, "SMILES"] = new["SMILES"][0]
#             comp = pcp.get_compounds(new["SMILES"][0], 'smiles')
#             try:
#                 if comp:
#                     for c in comp:
#                         merged_df["synonyms"][mer] = c.synonyms
#                         merged_df.loc[mer, "IUPAC"] = c.iupac_name
#                         merged_df.loc[mer, "Formula"] = c.molecular_formula
#             except Exception:
#                 pass
#             merged_df["superclass"][mer] = np.nan
#             merged_df["class"][mer] = np.nan
#             merged_df["subclass"][mer] = np.nan
#             merged_df["ClassificationSource"][mer] = np.nan
#         else:
#             pass
#     elif "HMDB" in merged_df["AnnotationSources"][mer]:
#         if len(df_count_4.loc[df_count_4["Source"] == "HMDB"]) == 1:
#             new = df_count_4.loc[df_count_4["Source"] == "HMDB"]
#             new.reset_index(drop=True, inplace=True)
#             merged_df.loc[mer, "SMILES"] = new["SMILES"][0]
#             comp = pcp.get_compounds(new["SMILES"][0], 'smiles')
#             try:
#                 if comp:
#                     for c in comp:
#                         merged_df["synonyms"][mer] = c.synonyms
#                         merged_df.loc[mer, "IUPAC"] = c.iupac_name
#                         merged_df.loc[mer, "Formula"] = c.molecular_formula
#             except Exception:
#                 pass
#             merged_df["superclass"][mer] = np.nan
#             merged_df["class"][mer] = np.nan
#             merged_df["subclass"][mer] = np.nan
#             merged_df["ClassificationSource"][mer] = np.nan
#         elif len(df_count_4.loc[df_count_4["Source"] == "HMDB"]) > 1:
#             new = df_count_4.loc[df_count_4["Source"] == "HMDB"]
#             new.reset_index(drop=True, inplace=True)
#             merged_df.loc[mer, "SMILES"] = new["SMILES"][0]
#             comp = pcp.get_compounds(new["SMILES"][0], 'smiles')
#             try:
#                 if comp:
#                     for c in comp:
#                         merged_df["synonyms"][mer] = c.synonyms
#                         merged_df.loc[mer, "IUPAC"] = c.iupac_name
#                         merged_df.loc[mer, "Formula"] = c.molecular_formula
#             except Exception:
#                 pass
#             merged_df["superclass"][mer] = np.nan
#             merged_df["class"][mer] = np.nan
#             merged_df["subclass"][mer] = np.nan
#             merged_df["ClassificationSource"][mer] = np.nan
#         else:
#             pass
    
    
    

    
    return merged_df


# In[ ]:


def checkSMILES_validity(input_dir, resultcsv):
    

    """checkSMILES_validity does exactly as the name says, using
    RDKit, whether the SMILES are invalid or have invalid
    chemistry

    Parameters:
    input_dir (str): This is the input directory where all the .mzML
    files and their respective result directories are stored.

    results: df from combine_CuratedR

    Returns:
    dataframe: with valid SMILES
    csv: "MetabolomicsResults/final_curation_with_validSMILES.csv"

    Usage:
    checkSMILES_validity(input_dir = "usr/project/", results)

    """
    results = pd.read_csv(resultcsv)
    # check validity of SMILES
    for i, row in results.iterrows():
        if not isNaN(results["SMILES"][i]):
            m = Chem.MolFromSmiles(results["SMILES"][i], sanitize=False)
            if m is None:
                results["SMILES"][i] = "invalid_SMILES"
            else:
                try:
                    Chem.SanitizeMol(m)
                except Exception:
                    results["SMILES"][i] = "invalid_chemistry"
#     results.to_csv(
#         input_dir + "MetabolomicsResults/final_curation_with_validSMILES.csv"
#     )
    return results



# In[5]:


def CandidateSelection_SimilarityandIdentity(input_dir, standards = False):

    
    # entry is all files and folders in input_dir
    for entry in os.listdir(input_dir):
        # if the entry is also a directory
        if os.path.isdir(os.path.join(input_dir, entry)):
            print(entry)

            # reach spectra_dereplication folder
            sub_dir_spec = input_dir + "/" + entry + "/spectral_dereplication/"
            # reach SIRIUS results
            sub_dir_sir = input_dir + "/" + entry + "/insilico/SIRIUS/"
            # new line
            if os.path.exists(sub_dir_spec) and os.path.exists(sub_dir_sir):

                # list of all csv files in the spectral dereplication foler
                spec_msp_csv = glob.glob(
                    input_dir + "/" + entry + "/spectral_dereplication" + "/*.csv"
                )
                # Sirius csv result file
                sir_msp_csv = input_dir + "/" + entry + "/insilico/MS1DATA.csv"

                # if both exist; which should be the case, even in case of 0 results
                if os.path.exists(sir_msp_csv) and os.path.exists(spec_msp_csv[0]):

                    # read both csv files
                    spec_msv = pd.read_csv(spec_msp_csv[0])
                    sir_msv = pd.read_csv(sir_msp_csv)

                    spec_msv = spec_msv[
                        [
                            "id_X",
                            "premz",
                            "rtmin",
                            "rtmax",
                            "rtmed",
                            "rtmean",
                            "col_eng",
                            "pol",
                            "int",
                            "source_file",
                            "mbank_results_csv",
                        ]
                    ]
                    sir_msv = sir_msv[
                        [
                            "id_X",
                            "premz",
                            "rtmed",
                            "rtmean",
                            "int",
                            "col_eng",
                            "pol",
                            "ms2Peaks",
                            "ms1Peaks",
                            "sirius_result_dir",
                        ]
                    ]

                    merged_df = sir_msv.merge(
                        spec_msv,
                        how="inner",
                        left_on=["premz", "rtmed", "rtmean", "int", "col_eng", "pol"],
                        right_on=["premz", "rtmed", "rtmean", "int", "col_eng", "pol"],
                    )
                    merged_df["Formula"] = np.nan
                    merged_df["SMILES"] = np.nan
                    merged_df["PubChemID"] = np.nan
                    merged_df["IUPAC"] = np.nan
                    merged_df["synonyms"] = np.nan
                    merged_df["AnnotationSources"] = np.nan
                    merged_df["AnnotationCount"] = np.nan
                    merged_df["MSILevel"] = np.nan
                    merged_df["MCSS"] = np.nan
                    merged_df["superclass"] = np.nan
                    merged_df["class"] = np.nan
                    merged_df["subclass"] = np.nan
                    merged_df["ClassificationSource"] = np.nan
                    
                    for_only_formula = []
                    for_formula_canopus = []
                    
                    can_selec_dir = input_dir + "/" + entry + "/Candidate_Selection"
                    if not os.path.isdir(can_selec_dir):
                        os.mkdir(can_selec_dir)
#                     cmn_dir = can_selec_dir + "/ChemicalMN"
#                     os.mkdir(cmn_dir)
#                     cand_list_dir = can_selec_dir + "/Candidate_Lists"
#                     os.mkdir(cmn_dir)      

                    for mer, rows in merged_df.iterrows():
                        print(merged_df["premz"][mer])
                        if not isNaN(merged_df["sirius_result_dir"][mer]):
                            sirius_csv = merged_df["sirius_result_dir"][mer].replace(
                                "./", input_dir + "/"
                            )
                            #print(sirius_csv)
                        else:
                            df = pd. DataFrame(list())
                            df. to_csv(input_dir + '/empty_csv.csv')
                            sirius_csv = input_dir + '/empty_csv.csv'

                        mbank_csv = merged_df["mbank_results_csv"][mer].replace(
                            "./", input_dir + "/"
                        )
                        gnps_csv = (
                            merged_df["mbank_results_csv"][mer]
                            .replace("./", input_dir + "/")
                            .replace("mbank", "gnps")
                            .replace("MassBank", "GNPS")
                        )
                        hmdb_csv = (
                            merged_df["mbank_results_csv"][mer]
                            .replace("./", input_dir + "/")
                            .replace("mbank", "hmdb")
                            .replace("MassBank", "HMDB")
                        )

                        if (
                            os.path.exists(sirius_csv)
                            and os.path.exists(gnps_csv)
                            and os.path.exists(mbank_csv)
                            and os.path.exists(hmdb_csv)
                        ):

                            sirius_df = pd.read_csv(sirius_csv)
                            if "formula" in sirius_csv:
                            
                                if len(sirius_df) > 0:
                                    
                                    if "smiles" not in sirius_df.columns and "molecularFormula" in sirius_df.columns:
                
                                        #merged_df.loc[mer, "Formula"] = sirius_df["molecularFormula"][0]
                                        #merged_df["AnnotationSources"][mer] = "SIRIUS-Formula"    

                                        index = mer
                                        Formula = sirius_df["molecularFormula"][0]
                                        for_only_formula.append(
                                            {
                                                "index": index,
                                                "Formula": Formula,
                                                "AnnotationSources": "SIRIUS-Formula"
                                            }
                                        )

                                        if "class" in sirius_df.columns:
                                            #merged_df.loc[mer, "superclass"] = sirius_df["superclass"][0]
                                            #merged_df.loc[mer, "class"] = sirius_df["class"][0]
                                            #merged_df.loc[mer, "subclass"] = sirius_df["subclass"][0]
                                            #merged_df.loc[mer, "ClassificationSource"] = "CANOPUS"
                                            #merged_df.loc[mer, "AnnotationSources"] = "SIRIUS-Formula|CANOPUS"
                                            index = mer
                                            Formula = sirius_df["molecularFormula"][0]
                                            superclass = sirius_df["superclass"][0]
                                            classes = sirius_df["class"][0]
                                            subclass = sirius_df["subclass"][0]
                                            for_formula_canopus.append(
                                                {
                                                    "index": index,
                                                    "Formula": Formula,
                                                    "superclass": superclass,
                                                    "class": classes,
                                                    "ClassificationSource": "CANOPUS",
                                                    "AnnotationSources": "SIRIUS-Formula|CANOPUS"
                                                }
                                            )

                                        sirius_df = []

                            if "structure" in sirius_csv:
                            
                                if len(sirius_df) > 0:

                                    if len(sirius_df) > 50:
                                        sirius_df = sirius_df[0:50]
                                    if "smiles" in sirius_df.columns:
                                        
                                        sirius_df = sirius_df.drop_duplicates("smiles")
                                        sirius_df = sirius_df.dropna(subset=["smiles"])
                                        merged_df["Formula"][mer] = sirius_df["molecularFormula"][0]
                                        merged_df.loc[mer, "PubChemID"] = sirius_df["pubchemids"][0]
                                        if "class" in sirius_df.columns:
                                            merged_df["superclass"][mer] = sirius_df["superclass"][0]
                                            merged_df["class"][mer] = sirius_df["class"][0]
                                            merged_df["subclass"][mer] = sirius_df["subclass"][0]
                                            merged_df["ClassificationSource"][mer] = "CANOPUS"
                                            
                            elif len(sirius_df) == 0:
                                #print("NO Structures")
                                merged_df["Formula"][mer] = np.nan
                                merged_df["superclass"][mer] = np.nan
                                merged_df["class"][mer] = np.nan
                                merged_df["subclass"][mer] = np.nan
                                merged_df["ClassificationSource"][mer] = np.nan
                    


                            mbank_df = pd.read_csv(mbank_csv)
                            if len(mbank_df) > 0:
                                mbank_df = mbank_df.drop_duplicates("MBSMILES")
                                mbank_df = mbank_df.dropna(subset=["MBSMILES"])

                            gnps_df = pd.read_csv(gnps_csv)
                            if len(gnps_df) > 0:
                                gnps_df = gnps_df.drop_duplicates("GNPSSMILES")
                                gnps_df = gnps_df.dropna(subset=["GNPSSMILES"])

                            hmdb_df = pd.read_csv(hmdb_csv)
                            if len(hmdb_df) > 0:
                                hmdb_df = hmdb_df.drop_duplicates("HMDBSMILES")
                                hmdb_df = hmdb_df.dropna(subset=["HMDBSMILES"])

                            # 1 SGHM
                            if (
                                len(sirius_df) > 0
                                and len(gnps_df) > 0
                                and len(mbank_df) > 0
                                and len(hmdb_df) > 0
                            ):
                                mbank_df["rank_ids"] = [
                                    "M_" + str(s + 1) for s in range(len(mbank_df))
                                ]

                                gnps_df["rank_ids"] = [
                                    "G_" + str(s + 1) for s in range(len(gnps_df))
                                ]

                                hmdb_df["rank_ids"] = [
                                    "H_" + str(s + 1) for s in range(len(hmdb_df))
                                ]

                                sirius_df["rank_ids"] = [
                                    "S_" + str(s) for s in sirius_df["rank"]
                                ]
                                sirius_df["Source"] = "SIRIUS"

                                source_l1 = [
                                    *(list(sirius_df["Source"])),
                                    *(list(gnps_df["Source"])),
                                    *(list(mbank_df["Source"])),
                                    *(list(hmdb_df["Source"])),
                                ]

                                rank_l2 = [
                                    *(list(sirius_df["rank_ids"])),
                                    *(list(gnps_df["rank_ids"])),
                                    *(list(mbank_df["rank_ids"])),
                                    *(list(hmdb_df["rank_ids"])),
                                ]

                                smiles_l3 = [
                                    *(list(sirius_df["smiles"])),
                                    *(list(gnps_df["GNPSSMILES"])),
                                    *(list(mbank_df["MBSMILES"])),
                                    *(list(hmdb_df["HMDBSMILES"])),
                                ]

                                sm = pd.DataFrame(
                                    list(zip(source_l1, rank_l2, smiles_l3)),
                                    columns=["Source", "ranks", "SMILES"],
                                )

                                df_edge = chemMN_CandidateSelection(sm)
                                if df_edge is not None:


                                    df_edge.to_csv(
                                        input_dir
                                        + "/"
                                        + entry
                                        + "/"
                                        + 'Candidate_Selection'
                                        + "/"
                                        + str(merged_df["premz"][mer])
                                        + "_ChemMNedges.tsv",
                                        sep="\t",
                                    )
                                one_candidate = one_candidate_selection(
                                    sm,
                                    sirius_df=sirius_df,
                                    mbank_df=mbank_df,
                                    gnps_df=gnps_df,
                                    hmdb_df=hmdb_df,
                                    Source="SGHM",
                                )

                                candidates_with_counts = add_count_column(one_candidate)
                                candidates_with_counts.to_csv(
                                    input_dir
                                    + "/"
                                    + entry
                                    + "/"
                                    + 'Candidate_Selection'
                                    + "/"
                                    + str(merged_df["premz"][mer])
                                    + "sorted_candidate_list.tsv",
                                    sep="\t",
                                )
                                if max(candidates_with_counts["Count"]) == 4:
                                    sources_4(candidates_with_counts, merged_df, mer, sirius_df)
                                if max(candidates_with_counts["Count"]) == 3:
                                    sources_3(candidates_with_counts, merged_df, mer, sirius_df)
                                if max(candidates_with_counts["Count"]) == 2:
                                    sources_2(candidates_with_counts, merged_df, mer, sirius_df)
                                if max(candidates_with_counts["Count"]) == 1:
                                    sources_1(candidates_with_counts, merged_df, mer, sirius_df)

                            # 2 SGM
                            elif (
                                len(sirius_df) > 0
                                and len(gnps_df) > 0
                                and len(mbank_df) > 0
                                and len(hmdb_df) == 0
                            ):

                                mbank_df["rank_ids"] = [
                                    "M_" + str(s + 1) for s in range(len(mbank_df))
                                ]

                                gnps_df["rank_ids"] = [
                                    "G_" + str(s + 1) for s in range(len(gnps_df))
                                ]

                                sirius_df["rank_ids"] = [
                                    "S_" + str(s) for s in sirius_df["rank"]
                                ]
                                sirius_df["Source"] = "SIRIUS"

                                source_l1 = [
                                    *(list(sirius_df["Source"])),
                                    *(list(gnps_df["Source"])),
                                    *(list(mbank_df["Source"])),
                                ]

                                rank_l2 = [
                                    *(list(sirius_df["rank_ids"])),
                                    *(list(gnps_df["rank_ids"])),
                                    *(list(mbank_df["rank_ids"])),
                                ]

                                smiles_l3 = [
                                    *(list(sirius_df["smiles"])),
                                    *(list(gnps_df["GNPSSMILES"])),
                                    *(list(mbank_df["MBSMILES"])),
                                ]

                                sm = pd.DataFrame(
                                    list(zip(source_l1, rank_l2, smiles_l3)),
                                    columns=["Source", "ranks", "SMILES"],
                                )

                                df_edge = chemMN_CandidateSelection(sm)
                                if df_edge is not None:


                                    df_edge.to_csv(
                                        input_dir
                                        + "/"
                                        + entry
                                        + "/"
                                        + 'Candidate_Selection'
                                        + "/"
                                        + str(merged_df["premz"][mer])
                                        + "_ChemMNedges.tsv",
                                        sep="\t",
                                    )
                                one_candidate = one_candidate_selection(
                                    sm,
                                    sirius_df=sirius_df,
                                    mbank_df=mbank_df,
                                    gnps_df=gnps_df,
                                    # hmdb_df = hmdb_df,
                                    Source="SGM",
                                )
                                candidates_with_counts = add_count_column(one_candidate)
                                candidates_with_counts.to_csv(
                                    input_dir
                                    + "/"
                                    + entry
                                    + "/"
                                    + 'Candidate_Selection'
                                    + "/"
                                    + str(merged_df["premz"][mer])
                                    + "sorted_candidate_list.tsv",
                                    sep="\t",
                                )

                                if max(candidates_with_counts["Count"]) == 3:
                                    sources_3(candidates_with_counts, merged_df, mer, sirius_df)
                                if max(candidates_with_counts["Count"]) == 2:
                                    sources_2(candidates_with_counts, merged_df, mer, sirius_df)
                                if max(candidates_with_counts["Count"]) == 1:
                                    sources_1(candidates_with_counts, merged_df, mer, sirius_df)

                            # 3 SHM
                            elif (
                                len(sirius_df) > 0
                                and len(gnps_df) == 0
                                and len(mbank_df) > 0
                                and len(hmdb_df) > 0
                            ):

                                mbank_df["rank_ids"] = [
                                    "M_" + str(s + 1) for s in range(len(mbank_df))
                                ]

                                # gnps_df["rank_ids"] = ["G_" + str(s+1) for s in range(len(gnps_df))]

                                hmdb_df["rank_ids"] = [
                                    "H_" + str(s + 1) for s in range(len(hmdb_df))
                                ]

                                sirius_df["rank_ids"] = [
                                    "S_" + str(s) for s in sirius_df["rank"]
                                ]
                                sirius_df["Source"] = "SIRIUS"

                                source_l1 = [
                                    *(list(sirius_df["Source"])),
                                    *(list(mbank_df["Source"])),
                                    *(list(hmdb_df["Source"])),
                                ]

                                rank_l2 = [
                                    *(list(sirius_df["rank_ids"])),
                                    *(list(mbank_df["rank_ids"])),
                                    *(list(hmdb_df["rank_ids"])),
                                ]

                                smiles_l3 = [
                                    *(list(sirius_df["smiles"])),
                                    *(list(mbank_df["MBSMILES"])),
                                    *(list(hmdb_df["HMDBSMILES"])),
                                ]

                                sm = pd.DataFrame(
                                    list(zip(source_l1, rank_l2, smiles_l3)),
                                    columns=["Source", "ranks", "SMILES"],
                                )

                                df_edge = chemMN_CandidateSelection(sm)
                                if df_edge is not None:


                                    df_edge.to_csv(
                                        input_dir
                                        + "/"
                                        + entry
                                        + "/"
                                        + 'Candidate_Selection'
                                        + "/"
                                        + str(merged_df["premz"][mer])
                                        + "_ChemMNedges.tsv",
                                        sep="\t",
                                    )
                                one_candidate = one_candidate_selection(
                                    sm,
                                    sirius_df=sirius_df,
                                    mbank_df=mbank_df,
                                    # gnps_df = gnps_df ,
                                    hmdb_df=hmdb_df,
                                    Source="SHM",
                                )
                                candidates_with_counts = add_count_column(one_candidate)
                                candidates_with_counts.to_csv(
                                    input_dir
                                    + "/"
                                    + entry
                                    + "/"
                                    + 'Candidate_Selection'
                                    + "/"
                                    + str(merged_df["premz"][mer])
                                    + "sorted_candidate_list.tsv",
                                    sep="\t",
                                )

                                if max(candidates_with_counts["Count"]) == 3:
                                    sources_3(candidates_with_counts, merged_df, mer, sirius_df)
                                if max(candidates_with_counts["Count"]) == 2:
                                    sources_2(candidates_with_counts, merged_df, mer, sirius_df)
                                if max(candidates_with_counts["Count"]) == 1:
                                    sources_1(candidates_with_counts, merged_df, mer, sirius_df)

                            # 4 SGH
                            elif (
                                len(sirius_df) > 0
                                and len(gnps_df) > 0
                                and len(mbank_df) == 0
                                and len(hmdb_df) > 0
                            ):
                                # mbank_df["rank_ids"] = ["M_" + str(s+1) for s in range(len(mbank_df))]

                                gnps_df["rank_ids"] = [
                                    "G_" + str(s + 1) for s in range(len(gnps_df))
                                ]

                                hmdb_df["rank_ids"] = [
                                    "H_" + str(s + 1) for s in range(len(hmdb_df))
                                ]

                                sirius_df["rank_ids"] = [
                                    "S_" + str(s) for s in sirius_df["rank"]
                                ]
                                sirius_df["Source"] = "SIRIUS"

                                source_l1 = [
                                    *(list(sirius_df["Source"])),
                                    *(list(gnps_df["Source"])),
                                    *(list(hmdb_df["Source"])),
                                ]

                                rank_l2 = [
                                    *(list(sirius_df["rank_ids"])),
                                    *(list(gnps_df["rank_ids"])),
                                    *(list(hmdb_df["rank_ids"])),
                                ]

                                smiles_l3 = [
                                    *(list(sirius_df["smiles"])),
                                    *(list(gnps_df["GNPSSMILES"])),
                                    *(list(hmdb_df["HMDBSMILES"])),
                                ]

                                sm = pd.DataFrame(
                                    list(zip(source_l1, rank_l2, smiles_l3)),
                                    columns=["Source", "ranks", "SMILES"],
                                )

                                df_edge = chemMN_CandidateSelection(sm)
                                if df_edge is not None:




                                    df_edge.to_csv(
                                        input_dir
                                        + "/"
                                        + entry
                                        + "/"
                                        + 'Candidate_Selection'
                                        + "/"
                                        + str(merged_df["premz"][mer])
                                        + "_ChemMNedges.tsv",
                                        sep="\t",
                                    )
                                one_candidate = one_candidate_selection(
                                    sm,
                                    sirius_df=sirius_df,
                                    # mbank_df = mbank_df,
                                    gnps_df=gnps_df,
                                    hmdb_df=hmdb_df,
                                    Source="SGH",
                                )
                                candidates_with_counts = add_count_column(one_candidate)
                                candidates_with_counts.to_csv(
                                    input_dir
                                    + "/"
                                    + entry
                                    + "/"
                                    + 'Candidate_Selection'
                                    + "/"
                                    + str(merged_df["premz"][mer])
                                    + "sorted_candidate_list.tsv",
                                    sep="\t",
                                )
                                if max(candidates_with_counts["Count"]) == 3:
                                    sources_3(candidates_with_counts, merged_df, mer, sirius_df)
                                if max(candidates_with_counts["Count"]) == 2:
                                    sources_2(candidates_with_counts, merged_df, mer, sirius_df)
                                if max(candidates_with_counts["Count"]) == 1:
                                    sources_1(candidates_with_counts, merged_df, mer, sirius_df)


                            # 5 GHM
                            elif (
                                len(sirius_df) == 0
                                and len(gnps_df) > 0
                                and len(mbank_df) > 0
                                and len(hmdb_df) > 0
                            ):
                                mbank_df["rank_ids"] = [
                                    "M_" + str(s + 1) for s in range(len(mbank_df))
                                ]

                                gnps_df["rank_ids"] = [
                                    "G_" + str(s + 1) for s in range(len(gnps_df))
                                ]

                                hmdb_df["rank_ids"] = [
                                    "H_" + str(s + 1) for s in range(len(hmdb_df))
                                ]

                                # sirius_df["rank_ids"] = ["S_" + str(s) for s in sirius_df["rank"]]
                                # sirius_df["Source"] = "SIRIUS"

                                source_l1 = [
                                    *(list(gnps_df["Source"])),
                                    *(list(mbank_df["Source"])),
                                    *(list(hmdb_df["Source"])),
                                ]

                                rank_l2 = [
                                    *(list(gnps_df["rank_ids"])),
                                    *(list(mbank_df["rank_ids"])),
                                    *(list(hmdb_df["rank_ids"])),
                                ]

                                smiles_l3 = [
                                    *(list(gnps_df["GNPSSMILES"])),
                                    *(list(mbank_df["MBSMILES"])),
                                    *(list(hmdb_df["HMDBSMILES"])),
                                ]

                                sm = pd.DataFrame(
                                    list(zip(source_l1, rank_l2, smiles_l3)),
                                    columns=["Source", "ranks", "SMILES"],
                                )

                                df_edge = chemMN_CandidateSelection(sm)
                                if df_edge is not None:


                                    df_edge.to_csv(
                                        input_dir
                                        + "/"
                                        + entry
                                        + "/"
                                        + 'Candidate_Selection'
                                        + "/"
                                        + str(merged_df["premz"][mer])
                                        + "_ChemMNedges.tsv",
                                        sep="\t",
                                    )
                                one_candidate = one_candidate_selection(
                                    sm,
                                    # sirius_df = sirius_df,
                                    mbank_df=mbank_df,
                                    gnps_df=gnps_df,
                                    hmdb_df=hmdb_df,
                                    Source="GHM",
                                )
                                candidates_with_counts = add_count_column(one_candidate)
                                candidates_with_counts.to_csv(
                                    input_dir
                                    + "/"
                                    + entry
                                    + "/"
                                    + 'Candidate_Selection'
                                    + "/"
                                    + str(merged_df["premz"][mer])
                                    + "sorted_candidate_list.tsv",
                                    sep="\t",
                                )

                                if max(candidates_with_counts["Count"]) == 3:
                                    sources_3(candidates_with_counts, merged_df, mer, sirius_df=None)
                                if max(candidates_with_counts["Count"]) == 2:
                                    sources_2(candidates_with_counts, merged_df, mer, sirius_df=None)
                                if max(candidates_with_counts["Count"]) == 1:
                                    sources_1(candidates_with_counts, merged_df, mer, sirius_df=None)

                            # 6 SG
                            elif (
                                len(sirius_df) > 0
                                and len(gnps_df) > 0
                                and len(mbank_df) == 0
                                and len(hmdb_df) == 0
                            ):
                                # mbank_df["rank_ids"] = ["M_" + str(s+1) for s in range(len(mbank_df))]

                                gnps_df["rank_ids"] = [
                                    "G_" + str(s + 1) for s in range(len(gnps_df))
                                ]

                                sirius_df["rank_ids"] = [
                                    "S_" + str(s) for s in sirius_df["rank"]
                                ]
                                sirius_df["Source"] = "SIRIUS"

                                source_l1 = [
                                    *(list(sirius_df["Source"])),
                                    *(list(gnps_df["Source"])),
                                ]
                                # ,*(list(mbank_df["Source"]))]

                                rank_l2 = [
                                    *(list(sirius_df["rank_ids"])),
                                    *(list(gnps_df["rank_ids"])),
                                ]
                                # ,*(list(mbank_df["rank_ids"]))]

                                smiles_l3 = [
                                    *(list(sirius_df["smiles"])),
                                    *(list(gnps_df["GNPSSMILES"])),
                                ]
                                # ,*(list(mbank_df["MBSMILES"]))]

                                sm = pd.DataFrame(
                                    list(zip(source_l1, rank_l2, smiles_l3)),
                                    columns=["Source", "ranks", "SMILES"],
                                )

                                df_edge = chemMN_CandidateSelection(sm)
                                if df_edge is not None:


                                    df_edge.to_csv(
                                        input_dir
                                        + "/"
                                        + entry
                                        + "/"
                                        + 'Candidate_Selection'
                                        + "/"
                                        + str(merged_df["premz"][mer])
                                        + "_ChemMNedges.tsv",
                                        sep="\t",
                                    )

                                one_candidate = one_candidate_selection(
                                    sm,
                                    sirius_df=sirius_df,
                                    # mbank_df = mbank_df,
                                    gnps_df=gnps_df,
                                    # hmdb_df = hmdb_df,
                                    Source="SG",
                                )
                                candidates_with_counts = add_count_column(one_candidate)
                                candidates_with_counts.to_csv(
                                    input_dir
                                    + "/"
                                    + entry
                                    + "/"
                                    + 'Candidate_Selection'
                                    + "/"
                                    + str(merged_df["premz"][mer])
                                    + "sorted_candidate_list.tsv",
                                    sep="\t",
                                )


                                if max(candidates_with_counts["Count"]) == 2:
                                    sources_2(candidates_with_counts, merged_df, mer, sirius_df)
                                if max(candidates_with_counts["Count"]) == 1:
                                    sources_1(candidates_with_counts, merged_df, mer, sirius_df)

                            # 7 SH
                            elif (
                                len(sirius_df) > 0
                                and len(gnps_df) == 0
                                and len(mbank_df) == 0
                                and len(hmdb_df) > 0
                            ):
                                # mbank_df["rank_ids"] = ["M_" + str(s+1) for s in range(len(mbank_df))]

                                # gnps_df["rank_ids"] = ["G_" + str(s+1) for s in range(len(gnps_df))]

                                hmdb_df["rank_ids"] = [
                                    "H_" + str(s + 1) for s in range(len(hmdb_df))
                                ]

                                sirius_df["rank_ids"] = [
                                    "S_" + str(s) for s in sirius_df["rank"]
                                ]
                                sirius_df["Source"] = "SIRIUS"

                                source_l1 = [
                                    *(list(sirius_df["Source"])),
                                    *(list(hmdb_df["Source"])),
                                ]

                                rank_l2 = [
                                    *(list(sirius_df["rank_ids"])),
                                    *(list(hmdb_df["rank_ids"])),
                                ]

                                smiles_l3 = [
                                    *(list(sirius_df["smiles"])),
                                    *(list(hmdb_df["HMDBSMILES"])),
                                ]

                                sm = pd.DataFrame(
                                    list(zip(source_l1, rank_l2, smiles_l3)),
                                    columns=["Source", "ranks", "SMILES"],
                                )
                                df_edge = chemMN_CandidateSelection(sm)
                                if df_edge is not None:


                                    df_edge.to_csv(
                                        input_dir
                                        + "/"
                                        + entry
                                        + "/"
                                        + 'Candidate_Selection'
                                        + "/"
                                        + str(merged_df["premz"][mer])
                                        + "_ChemMNedges.tsv",
                                        sep="\t",
                                    )

                                one_candidate = one_candidate_selection(
                                    sm,
                                    sirius_df=sirius_df,
                                    # mbank_df = mbank_df,
                                    # gnps_df = gnps_df ,
                                    hmdb_df=hmdb_df,
                                    Source="SH",
                                )
                                candidates_with_counts = add_count_column(one_candidate)
                                candidates_with_counts.to_csv(
                                    input_dir
                                    + "/"
                                    + entry
                                    + "/"
                                    + 'Candidate_Selection'
                                    + "/"
                                    + str(merged_df["premz"][mer])
                                    + "sorted_candidate_list.tsv",
                                    sep="\t",
                                )

                                if max(candidates_with_counts["Count"]) == 2:
                                    sources_2(candidates_with_counts, merged_df, mer, sirius_df)
                                if max(candidates_with_counts["Count"]) == 1:
                                    sources_1(candidates_with_counts, merged_df, mer, sirius_df)

                            # 8 SM
                            elif (
                                len(sirius_df) > 0
                                and len(gnps_df) == 0
                                and len(mbank_df) > 0
                                and len(hmdb_df) == 0
                            ):
                                mbank_df["rank_ids"] = [
                                    "M_" + str(s + 1) for s in range(len(mbank_df))
                                ]

                                # gnps_df["rank_ids"] = ["G_" + str(s+1) for s in range(len(gnps_df))]

                                # hmdb_df["rank_ids"] = ["H_" + str(s+1) for s in range(len(hmdb_df))]

                                sirius_df["rank_ids"] = [
                                    "S_" + str(s) for s in sirius_df["rank"]
                                ]
                                sirius_df["Source"] = "SIRIUS"

                                source_l1 = [
                                    *(list(sirius_df["Source"])),
                                    *(list(mbank_df["Source"])),
                                ]

                                rank_l2 = [
                                    *(list(sirius_df["rank_ids"])),
                                    *(list(mbank_df["rank_ids"])),
                                ]

                                smiles_l3 = [
                                    *(list(sirius_df["smiles"])),
                                    *(list(mbank_df["MBSMILES"])),
                                ]

                                sm = pd.DataFrame(
                                    list(zip(source_l1, rank_l2, smiles_l3)),
                                    columns=["Source", "ranks", "SMILES"],
                                )

                                df_edge = chemMN_CandidateSelection(sm)
                                if df_edge is not None:


                                    df_edge.to_csv(
                                        input_dir
                                        + "/"
                                        + entry
                                        + "/"
                                        + 'Candidate_Selection'
                                        + "/"
                                        + str(merged_df["premz"][mer])
                                        + "_ChemMNedges.tsv",
                                        sep="\t",
                                    )

                                one_candidate = one_candidate_selection(
                                    sm,
                                    sirius_df=sirius_df,
                                    mbank_df=mbank_df,
                                    # gnps_df = gnps_df ,
                                    # hmdb_df = hmdb_df,
                                    Source="SM",
                                )
                                candidates_with_counts = add_count_column(one_candidate)
                                candidates_with_counts.to_csv(
                                    input_dir
                                    + "/"
                                    + entry
                                    + "/"
                                    + 'Candidate_Selection'
                                    + "/"
                                    + str(merged_df["premz"][mer])
                                    + "sorted_candidate_list.tsv",
                                    sep="\t",
                                )
                                if max(candidates_with_counts["Count"]) == 2:
                                    sources_2(candidates_with_counts, merged_df, mer, sirius_df)
                                if max(candidates_with_counts["Count"]) == 1:
                                    sources_1(candidates_with_counts, merged_df, mer, sirius_df)

                            # 9 GM
                            elif (
                                len(sirius_df) == 0
                                and len(gnps_df) > 0
                                and len(mbank_df) > 0
                                and len(hmdb_df) == 0
                            ):
                                mbank_df["rank_ids"] = [
                                    "M_" + str(s + 1) for s in range(len(mbank_df))
                                ]

                                gnps_df["rank_ids"] = [
                                    "G_" + str(s + 1) for s in range(len(gnps_df))
                                ]

                                # hmdb_df["rank_ids"] = ["H_" + str(s+1) for s in range(len(hmdb_df))]

                                # sirius_df["rank_ids"] = ["S_" + str(s) for s in sirius_df["rank"]]
                                # sirius_df["Source"] = "SIRIUS"

                                source_l1 = [
                                    *(list(mbank_df["Source"])),
                                    *(list(gnps_df["Source"])),
                                ]

                                rank_l2 = [
                                    *(list(mbank_df["rank_ids"])),
                                    *(list(gnps_df["rank_ids"])),
                                ]

                                smiles_l3 = [
                                    *(list(mbank_df["MBSMILES"])),
                                    *(list(gnps_df["GNPSSMILES"])),
                                ]

                                sm = pd.DataFrame(
                                    list(zip(source_l1, rank_l2, smiles_l3)),
                                    columns=["Source", "ranks", "SMILES"],
                                )

                                df_edge = chemMN_CandidateSelection(sm)
                                if df_edge is not None:


                                    df_edge.to_csv(
                                        input_dir
                                        + "/"
                                        + entry
                                        + "/"
                                        + 'Candidate_Selection'
                                        + "/"
                                        + str(merged_df["premz"][mer])
                                        + "_ChemMNedges.tsv",
                                        sep="\t",
                                    )

                                one_candidate = one_candidate_selection(
                                    sm,
                                    # sirius_df = sirius_df,
                                    mbank_df=mbank_df,
                                    gnps_df=gnps_df,
                                    # hmdb_df = hmdb_df,
                                    Source="GM",
                                )
                                candidates_with_counts = add_count_column(one_candidate)
                                candidates_with_counts.to_csv(
                                    input_dir
                                    + "/"
                                    + entry
                                    + "/"
                                    + 'Candidate_Selection'
                                    + "/"
                                    + str(merged_df["premz"][mer])
                                    + "sorted_candidate_list.tsv",
                                    sep="\t",
                                )
                                if max(candidates_with_counts["Count"]) == 2:
                                    sources_2(candidates_with_counts, merged_df, mer, sirius_df=None)
                                if max(candidates_with_counts["Count"]) == 1:
                                    sources_1(candidates_with_counts, merged_df, mer, sirius_df=None)

                            # 10 GH
                            elif (
                                len(sirius_df) == 0
                                and len(gnps_df) > 0
                                and len(mbank_df) == 0
                                and len(hmdb_df) > 0
                            ):
                                # mbank_df["rank_ids"] = ["M_" + str(s+1) for s in range(len(mbank_df))]

                                gnps_df["rank_ids"] = [
                                    "G_" + str(s + 1) for s in range(len(gnps_df))
                                ]

                                hmdb_df["rank_ids"] = [
                                    "H_" + str(s + 1) for s in range(len(hmdb_df))
                                ]

                                # sirius_df["rank_ids"] = ["S_" + str(s) for s in sirius_df["rank"]]
                                # sirius_df["Source"] = "SIRIUS"

                                source_l1 = [
                                    *(list(gnps_df["Source"])),
                                    *(list(hmdb_df["Source"])),
                                ]

                                rank_l2 = [
                                    *(list(gnps_df["rank_ids"])),
                                    *(list(hmdb_df["rank_ids"])),
                                ]

                                smiles_l3 = [
                                    *(list(gnps_df["GNPSSMILES"])),
                                    *(list(hmdb_df["HMDBSMILES"])),
                                ]

                                sm = pd.DataFrame(
                                    list(zip(source_l1, rank_l2, smiles_l3)),
                                    columns=["Source", "ranks", "SMILES"],
                                )
                                df_edge = chemMN_CandidateSelection(sm)
                                if df_edge is not None:


                                    df_edge.to_csv(
                                        input_dir
                                        + "/"
                                        + entry
                                        + "/"
                                        + 'Candidate_Selection'
                                        + "/"
                                        + str(merged_df["premz"][mer])
                                        + "_ChemMNedges.tsv",
                                        sep="\t",
                                    )
                                one_candidate = one_candidate_selection(
                                    sm,
                                    # sirius_df = sirius_df,
                                    # mbank_df = mbank_df,
                                    gnps_df=gnps_df,
                                    hmdb_df=hmdb_df,
                                    Source="GH",
                                )
                                candidates_with_counts = add_count_column(one_candidate)
                                candidates_with_counts.to_csv(
                                    input_dir
                                    + "/"
                                    + entry
                                    + "/"
                                    + 'Candidate_Selection'
                                    + "/"
                                    + str(merged_df["premz"][mer])
                                    + "sorted_candidate_list.tsv",
                                    sep="\t",
                                )

                                if max(candidates_with_counts["Count"]) == 2:
                                    sources_2(candidates_with_counts, merged_df, mer, sirius_df=None)
                                if max(candidates_with_counts["Count"]) == 1:
                                    sources_1(candidates_with_counts, merged_df, mer, sirius_df=None)
                            # 11 HM
                            elif (
                                len(sirius_df) == 0
                                and len(gnps_df) == 0
                                and len(mbank_df) > 0
                                and len(hmdb_df) > 0
                            ):
                                mbank_df["rank_ids"] = [
                                    "M_" + str(s + 1) for s in range(len(mbank_df))
                                ]

                                # gnps_df["rank_ids"] = ["G_" + str(s+1) for s in range(len(gnps_df))]

                                hmdb_df["rank_ids"] = [
                                    "H_" + str(s + 1) for s in range(len(hmdb_df))
                                ]

                                # sirius_df["rank_ids"] = ["S_" + str(s) for s in sirius_df["rank"]]
                                # sirius_df["Source"] = "SIRIUS"

                                source_l1 = [
                                    *(list(mbank_df["Source"])),
                                    *(list(hmdb_df["Source"])),
                                ]

                                rank_l2 = [
                                    *(list(mbank_df["rank_ids"])),
                                    *(list(hmdb_df["rank_ids"])),
                                ]

                                smiles_l3 = [
                                    *(list(mbank_df["MBSMILES"])),
                                    *(list(hmdb_df["HMDBSMILES"])),
                                ]

                                sm = pd.DataFrame(
                                    list(zip(source_l1, rank_l2, smiles_l3)),
                                    columns=["Source", "ranks", "SMILES"],
                                )

                                df_edge = chemMN_CandidateSelection(sm)
                                if df_edge is not None:

                                    df_edge.to_csv(
                                        input_dir
                                        + "/"
                                        + entry
                                        + "/"
                                        + 'Candidate_Selection'
                                        + "/"
                                        + str(merged_df["premz"][mer])
                                        + "_ChemMNedges.tsv",
                                        sep="\t",
                                    )
                                one_candidate = one_candidate_selection(
                                    sm,
                                    # sirius_df = sirius_df,
                                    mbank_df=mbank_df,
                                    # gnps_df = gnps_df ,
                                    hmdb_df=hmdb_df,
                                    Source="HM",
                                )
                                candidates_with_counts = add_count_column(one_candidate)
                                candidates_with_counts.to_csv(
                                    input_dir
                                    + "/"
                                    + entry
                                    + "/"
                                    + 'Candidate_Selection'
                                    + "/"
                                    + str(merged_df["premz"][mer])
                                    + "sorted_candidate_list.tsv",
                                    sep="\t",
                                )
                                if max(candidates_with_counts["Count"]) == 2:
                                    sources_2(candidates_with_counts, merged_df, mer, sirius_df=None)
                                if max(candidates_with_counts["Count"]) == 1:
                                    sources_1(candidates_with_counts, merged_df, mer, sirius_df=None)

                            # S
                            elif (
                                len(sirius_df) > 0
                                and len(gnps_df) == 0
                                and len(mbank_df) == 0
                                and len(hmdb_df) == 0
                            ):

                                sirius_df["rank_ids"] = [
                                    "S_" + str(s) for s in sirius_df["rank"]
                                ]
                                sirius_df["Source"] = "SIRIUS"

                                source_l1 = [*(list(sirius_df["Source"]))]

                                rank_l2 = [*(list(sirius_df["rank_ids"]))]

                                smiles_l3 = [*(list(sirius_df["smiles"]))]

                                sm = pd.DataFrame(
                                    list(zip(source_l1, rank_l2, smiles_l3)),
                                    columns=["Source", "ranks", "SMILES"],
                                )

                                df_edge = chemMN_CandidateSelection(sm)

                                if df_edge is not None:

                                    df_edge.to_csv(
                                        input_dir
                                        + "/"
                                        + entry
                                        + "/"
                                        + 'Candidate_Selection'
                                        + "/"
                                        + str(merged_df["premz"][mer])
                                        + "_ChemMNedges.tsv",
                                        sep="\t",
                                    )

                                one_candidate = one_candidate_selection(
                                    sm,
                                    sirius_df=sirius_df,
                                    # mbank_df = mbank_df,
                                    # gnps_df = gnps_df ,
                                    # hmdb_df = hmdb_df,
                                    Source="S",
                                )
                                candidates_with_counts = add_count_column(one_candidate)
                                candidates_with_counts.to_csv(
                                    input_dir
                                    + "/"
                                    + entry
                                    + "/"
                                    + 'Candidate_Selection'
                                    + "/"
                                    + str(merged_df["premz"][mer])
                                    + "sorted_candidate_list.tsv",
                                    sep="\t",
                                )

                                if max(candidates_with_counts["Count"]) == 1:
                                    sources_1(candidates_with_counts, merged_df, mer, sirius_df)
                            # G
                            elif (
                                len(sirius_df) == 0
                                and len(gnps_df) > 0
                                and len(mbank_df) == 0
                                and len(hmdb_df) == 0
                            ):
                                gnps_df["rank_ids"] = [
                                    "G_" + str(s + 1) for s in range(len(gnps_df))
                                ]

                                source_l1 = [*(list(gnps_df["Source"]))]

                                rank_l2 = [*(list(gnps_df["rank_ids"]))]

                                smiles_l3 = [*(list(gnps_df["GNPSSMILES"]))]

                                sm = pd.DataFrame(
                                    list(zip(source_l1, rank_l2, smiles_l3)),
                                    columns=["Source", "ranks", "SMILES"],
                                )

                                df_edge = chemMN_CandidateSelection(sm)
                                if df_edge is not None:

                                    df_edge.to_csv(
                                        input_dir
                                        + "/"
                                        + entry
                                        + "/"
                                        + 'Candidate_Selection'
                                        + "/"
                                        + str(merged_df["premz"][mer])
                                        + "_ChemMNedges.tsv",
                                        sep="\t",
                                    )

                                one_candidate = one_candidate_selection(
                                    sm,
                                    # sirius_df = sirius_df,
                                    # mbank_df = mbank_df,
                                    gnps_df=gnps_df,
                                    # hmdb_df = hmdb_df,
                                    Source="G",
                                )
                                candidates_with_counts = add_count_column(one_candidate)
                                candidates_with_counts.to_csv(
                                    input_dir
                                    + "/"
                                    + entry
                                    + "/"
                                    + 'Candidate_Selection'
                                    + "/"
                                    + str(merged_df["premz"][mer])
                                    + "sorted_candidate_list.tsv",
                                    sep="\t",
                                )
                                if max(candidates_with_counts["Count"]) == 1:
                                    sources_1(candidates_with_counts, merged_df, mer, sirius_df=None)
                            # M
                            elif (
                                len(sirius_df) == 0
                                and len(gnps_df) == 0
                                and len(mbank_df) > 0
                                and len(hmdb_df) == 0
                            ):
                                mbank_df["rank_ids"] = [
                                    "M_" + str(s + 1) for s in range(len(mbank_df))
                                ]

                                source_l1 = [*(list(mbank_df["Source"]))]

                                rank_l2 = [*(list(mbank_df["rank_ids"]))]

                                smiles_l3 = [*(list(mbank_df["MBSMILES"]))]

                                sm = pd.DataFrame(
                                    list(zip(source_l1, rank_l2, smiles_l3)),
                                    columns=["Source", "ranks", "SMILES"],
                                )

                                df_edge = chemMN_CandidateSelection(sm)
                                if df_edge is not None:

                                    df_edge.to_csv(
                                        input_dir
                                        + "/"
                                        + entry
                                        + "/"
                                        + 'Candidate_Selection'
                                        + "/"
                                        + str(merged_df["premz"][mer])
                                        + "_ChemMNedges.tsv",
                                        sep="\t",
                                    )

                                one_candidate = one_candidate_selection(
                                    sm,
                                    # sirius_df = sirius_df,
                                    mbank_df=mbank_df,
                                    # gnps_df = gnps_df ,
                                    # hmdb_df = hmdb_df,
                                    Source="M",
                                )
                                candidates_with_counts = add_count_column(one_candidate)
                                candidates_with_counts.to_csv(
                                    input_dir
                                    + "/"
                                    + entry
                                    + "/"
                                    + 'Candidate_Selection'
                                    + "/"
                                    + str(merged_df["premz"][mer])
                                    + "sorted_candidate_list.tsv",
                                    sep="\t",
                                )
                                if max(candidates_with_counts["Count"]) == 1:
                                    sources_1(candidates_with_counts, merged_df, mer, sirius_df=None)
                            # H
                            elif (
                                len(sirius_df) == 0
                                and len(gnps_df) == 0
                                and len(mbank_df) == 0
                                and len(hmdb_df) > 0
                            ):
                                hmdb_df["rank_ids"] = [
                                    "H_" + str(s + 1) for s in range(len(hmdb_df))
                                ]

                                source_l1 = [*(list(hmdb_df["Source"]))]

                                rank_l2 = [*(list(hmdb_df["rank_ids"]))]

                                smiles_l3 = [*(list(hmdb_df["HMDBSMILES"]))]

                                sm = pd.DataFrame(
                                    list(zip(source_l1, rank_l2, smiles_l3)),
                                    columns=["Source", "ranks", "SMILES"],
                                )

                                df_edge = chemMN_CandidateSelection(sm)
                                if df_edge is not None:


                                    df_edge.to_csv(
                                        input_dir
                                        + "/"
                                        + entry
                                        + "/"
                                        + 'Candidate_Selection'
                                        + "/"
                                        + str(merged_df["premz"][mer])
                                        + "_ChemMNedges.tsv",
                                        sep="\t",
                                    )

                                one_candidate = one_candidate_selection(
                                    sm,
                                    #                                                                         sirius_df = sirius_df,
                                    #                                                                         mbank_df = mbank_df,
                                    #                                                                         gnps_df = gnps_df ,
                                    hmdb_df=hmdb_df,
                                    Source="H",
                                )
                                candidates_with_counts = add_count_column(one_candidate)
                                candidates_with_counts.to_csv(
                                    input_dir
                                    + "/"
                                    + entry
                                    + "/"
                                    + 'Candidate_Selection'
                                    + "/"
                                    + str(merged_df["premz"][mer])
                                    + "sorted_candidate_list.tsv",
                                    sep="\t",
                                )

                                if max(candidates_with_counts["Count"]) == 1:
                                    sources_1(candidates_with_counts, merged_df, mer, sirius_df=None)

                        if standards:
                            if not isNaN(merged_df["SMILES"][mer]):

                                merged_df.loc[mer, "MSILevel"] = 1

                    
                    
                    #### code for formula and canopus classes
                    df_for_formula=pd.DataFrame(for_only_formula)
                    df_for_formula_n_canopus = pd.DataFrame(for_formula_canopus)
                    
                    
                    
                    for m, row in merged_df.iterrows():
                        for f, row in df_for_formula.iterrows():
                            if m == df_for_formula["index"][f]:
                                merged_df.loc[m, "Formula"] = df_for_formula["Formula"][f]
                                merged_df.loc[m, "AnnotationSources"] = df_for_formula["AnnotationSources"][f]

                        for c, row in df_for_formula_n_canopus.iterrows():
                            if m == df_for_formula_n_canopus["index"][f]:
                                merged_df.loc[m, "Formula"] = df_for_formula_n_canopus["Formula"][f]
                                merged_df.loc[m, "subclass"] = df_for_formula_n_canopus["subclass"][f]
                                merged_df.loc[m, "class"] = df_for_formula_n_canopus["class"][f]
                                merged_df.loc[m, "superclass"] = df_for_formula_n_canopus["superclass"][f]
                                merged_df.loc[m, "ClassificationSource"] = df_for_formula_n_canopus["ClassificationSource"][f]
                                merged_df.loc[m, "Formula"] = df_for_formula_n_canopus["Formula"][f]
                                merged_df.loc[m, "AnnotationSources"] = df_for_formula_n_canopus["AnnotationSources"][f]


                    
                    merged_df.to_csv(
                        input_dir
                        + "/"
                        + entry
                        + "/"
                        + "mergedResults-with-one-Candidates.csv"
                    )
                    #print(merged_df)
                    merged_df = checkSMILES_validity(input_dir, resultcsv = input_dir
                        + "/"
                        + entry
                        + "/"
                        + "mergedResults-with-one-Candidates.csv")
                    merged_df.to_csv(
                        input_dir
                        + "/"
                        + entry
                        + "/"
                        + "mergedResults-with-one-Candidates.csv"
                    )
                        


# In[32]:


def merge_all_results(input_dir):
    names = []
    
    # entry is all files and folders in input_dir
    for entry in os.listdir(input_dir):
        # if the entry is also a directory
        if os.path.isdir(os.path.join(input_dir, entry)):

            # reach spectra_dereplication folder
            merged_file_res = input_dir+ "/"+ entry+ "/mergedResults-with-one-Candidates.csv"

            if os.path.exists(merged_file_res):
                merged_csv = pd.read_csv(merged_file_res)
                names.append(merged_csv)
    merged = (pd.concat(names, ignore_index=True))
    # Select the ones you want
    merged = merged[['id_X_x', 'premz', 'rtmed', 'rtmean',
           'int', 'col_eng', 'pol', 'rtmin', 'rtmax', 'source_file', 'Formula',
                     'PubChemID',
                    'SMILES', 'IUPAC', 'synonyms', 'AnnotationSources', 'AnnotationCount',
           'MSILevel', 'MCSS', 'subclass', 'class', 'superclass', 'ClassificationSource']]
    merged.to_csv(input_dir + "/final_candidates.csv")
    return(merged)


# In[24]:


# Comparison with a list of SMILES from any Source
def SMILESscreening(input_dir, resultcsv, complist, listname):
   
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


# In[ ]:





# In[ ]:


def classification(input_dir, resultcsv):
   

    """classification function uses ClassyFire ChemONT

    Parameters:
    input_dir (str): This is the input directory where all the .mzML
    files and their respective result directories are stored.

    resultcsv: csv of df from combine_CuratedR or checkSMILES_validity

    Returns:
    dataframe: with classification
    csv: "MetabolomicsResults/final_curationList.csv"

    Usage:
    checkSMILES_validity(input_dir = "usr/project/", frame)

    """
    
    frame = pd.read_csv(resultcsv)
    inchis = []
    for i, row in frame.iterrows():
        if not isNaN(frame["SMILES"][i]):
            if "SIRIUS" not in frame["AnnotationSources"][i]:
                try:
                    InChI = Chem.MolToInchi(Chem.MolFromSmiles(frame["SMILES"][i]))
                    InChIKey = Chem.inchi.InchiToInchiKey(InChI)
                    inchis.append(
                        {
                            "index": i,
                            "smiles": frame["SMILES"][i],
                            "inchi": InChI,
                            "inchikey": InChIKey,
                        }
                    )
                except Exception:
                    pass
            elif "SIRIUS" in frame["AnnotationSources"][i]:
                if isNaN(frame["superclass"][i]):
                    try:
                        InChI = Chem.MolToInchi(Chem.MolFromSmiles(frame["SMILES"][i]))
                        InChIKey = Chem.inchi.InchiToInchiKey(InChI)
                        inchis.append(
                            {
                                "index": i,
                                "smiles": frame["SMILES"][i],
                                "inchi": InChI,
                                "inchikey": InChIKey,
                            }
                        )
                    except Exception:
                        pass
    inchis = pd.DataFrame(inchis)
    if len(inchis):
        inchis = inchis.loc[-isNaN(inchis["inchikey"])]
        # Retrieve ClassyFire classifications

        # This first step is done using inchikey and interrogation of the gnps classified structures
        """
        gnps_proxy = True
        url = "http://classyfire.wishartlab.com"
        proxy_url = "https://gnps-classyfire.ucsd.edu"
        chunk_size = 1000
        sleep_interval = 12
        """

        all_inchi_keys = list(inchis["inchikey"].drop_duplicates())

        resolved_ik_number_list = [0, 0]
        # total_inchikey_number = len(all_inchi_keys)

        while True:

            # start_time = time.time()

            # print('%s inchikey to resolve' % total_inchikey_number )
            get_classifications_cf_mod(all_inchi_keys, par_level=6)

            cleanse("all_json.json", "all_json.json")

            with open("all_json.json") as tweetfile:
                jsondic = json.loads(tweetfile.read())

            df = json_normalize(jsondic)
            df = df.drop_duplicates("inchikey")
            resolved_ik_number = len(df.drop_duplicates("inchikey").inchikey)
            resolved_ik_number_list.append(resolved_ik_number)
            # print('%s resolved inchikeys' % resolved_ik_number )
            # print("done in --- %s seconds ---" % (time.time() - start_time))

            if (
                resolved_ik_number_list[-1] < resolved_ik_number_list[-2]
                or resolved_ik_number_list[-1] == resolved_ik_number_list[-3]
            ):
                break
            cleanse("all_json.json", "all_json_cleaned.json")

            with open("all_json_cleaned.json") as tweetfile:
                jsondic = json.loads(tweetfile.read())

        flattened_classified_json = json_normalize(jsondic)
        flattened_df = flattened_classified_json.drop_duplicates("inchikey")
        flattened_df["inchikey"] = flattened_df["inchikey"].str.replace(
            r"InChIKey=", ""
        )
        df_merged = pd.merge(
            inchis, flattened_df, left_on="inchikey", right_on="inchikey", how="left"
        )

        for p, rowp in df_merged.iterrows():
            for q, rowq in frame.iterrows():
                if df_merged["smiles_x"][p] is frame["SMILES"][q]:
                    frame.loc[q, "subclass"] = df_merged["subclass.name"][p]
                    frame.loc[q, "class"] = df_merged["class.name"][p]
                    frame.loc[q, "superclass"] = df_merged["superclass.name"][p]
                    frame.loc[q, "ClassificationSource"] = "ClassyFire"

        frame.to_csv(input_dir + "/final_candidates_classes.csv")
        return frame


# In[19]:


# NP_Classifier classification
def Np_pathways(input_dir, resultcsv):
    df = pd.read_csv(resultcsv)
    npresults = []
    for i, row in df.iterrows():
        if not isNaN(df["SMILES"][i]):
            try:
                cvv = Chem.MolFromSmiles(df["SMILES"][i])
                cvv = Chem.MolToSmiles(cvv, isomericSmiles=False)
                c = urllib.parse.quote_plus(cvv, safe=" ")

                url = "https://npclassifier.ucsd.edu/classify?smiles=" + c
                names = str(df["id_X"][i])
                outx = str("myFile" + names + ".txt")
                file = wget.download(url, out=outx)
                a_dataframe = pd.read_csv(file, delimiter="]")
                xox = list(a_dataframe.columns.values)
                splitting0 = xox[0].split(":")
                xoc = re.sub('\ |\[|\]|"', " ", splitting0[1]).strip()
                splitting1 = xox[1].split(":")
                xos = re.sub('\ |\[|\]|"', " ", splitting1[1]).strip()
                # except:
                # splitting1 = xox[1].split(':')
                # xos = re.sub('\ |\[|\]|\"', '', splitting1[0])
                splitting2 = xox[2].split(":")
                xop = re.sub('\ |\[|\]|"', " ", splitting2[1]).strip()
                # df.loc[i, 'npclass'] = xoc
                # df.loc[i, 'npsuper_class'] = xos
                if not isNaN(df["class"][i]) and df["class"][i] in xoc:
                    df.loc[i, "np_pathway"] = xop
                os.remove(outx)
                time.sleep(0.5)

                npresults.append(
                    {
                        "index": i,
                        # 'id': df['file_id'][i],
                        "mz": df["premz"][i],
                        "rt": df["rtmed"][i],
                        "SMILES": df["SMILES"][i],
                        "class": xoc,
                        "subclass": xos,
                        "pathway": xop,
                    }
                )
                if df["class"][i] == xoc:
                    df.loc[i, "pathway"]
            except Exception:
                pass
    np_results = pd.DataFrame(npresults)
    np_results.to_csv(input_dir + "/NPClassifier_Results.csv")
    df.to_csv(input_dir + "/final_results_with_Pathways.csv")


# Chemical Similarity MN
def chemMN(input_dir, resultcsv):
    # read csv
    df = pd.read_csv(resultcsv)

    # define empty variable
    dbn = []

    # check the result csv
    for i, row in df.iterrows():
        # to compare each element with each opther element
        for j, row in df.iterrows():

            # if its not same id
            if df["SMILES"][i] != df["SMILES"][j]:

                if not isNaN(df["SMILES"][i]):
                    if not isNaN(df["SMILES"][j]):

                        try:
                            ms = [
                                Chem.MolFromSmiles(df["SMILES"][i]),
                                Chem.MolFromSmiles(df["SMILES"][j]),
                            ]
                            fps = [
                                AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=2048)
                                for x in ms
                            ]
                            tn = DataStructs.FingerprintSimilarity(fps[0], fps[1])
                            dbn.append(
                                {
                                    "Name_i": df["IUPAC"][i],
                                    "Name_j": df["IUPAC"][j],
                                    "i": df["SMILES"][i],
                                    "j": df["SMILES"][j],
                                    "Tanimoto": tn,
                                }
                            )
                        except Exception:
                            pass
    # save chemical similarities
    db_edgenode = pd.DataFrame(dbn)

    dfe = []
    heavy_atoms = ["C", "N", "P", "O", "S"]
    for i, row in db_edgenode.iterrows():
        if db_edgenode["Tanimoto"][i] >= 0.85:
            # list of mol used to calaculate the MCSS
            n = [
                Chem.MolFromSmiles(db_edgenode["i"][i]),
                Chem.MolFromSmiles(db_edgenode["j"][i]),
            ]
            res = rdFMCS.FindMCS(n, timeout=60)
            sm_res = Chem.MolToSmiles(Chem.MolFromSmarts(res.smartsString))
            # Check if the MCSS has one of the heavy atoms and whether they are
            # more than 3
            elem = [ele for ele in heavy_atoms if (ele in sm_res)]
            if elem and len(sm_res) >= 3:
                MCSS_SMILES = Chem.MolToSmiles(Chem.MolFromSmarts(res.smartsString))

            dfe.append(
                {
                    "Start": db_edgenode["Name_i"][i],
                    "End": db_edgenode["Name_j"][i],
                    "Tanimoto": db_edgenode["Tanimoto"][i],
                    "Start_SMILES": db_edgenode["i"][i],
                    "End_SMILES": db_edgenode["j"][i],
                    # 'Start_Source':db_edgenode['Source_i'][i],
                    # 'End_Source':db_edgenode['Source_j'][i],
                    "MCSS": MCSS_SMILES,
                }
            )
    if len(df_edge) > 0:
        
        df_edge = pd.DataFrame(dfe)
        df_edge["Start"] = df_edge["Start"].astype(str)
        df_edge["End"] = df_edge["End"].astype(str)
        df_edge["sorted_row"] = [sorted([a, b]) for a, b in zip(df_edge.Start, df_edge.End)]
        df_edge["sorted_row"] = df_edge["sorted_row"].astype(str)
        df_edge.drop_duplicates(subset=["sorted_row"], inplace=True)

        nodes = []
        for i, row in df.iterrows():
            n = df["IUPAC"][i]
            nodes.append({"nodes": n})

        node = pd.DataFrame(nodes)

        df_edge.to_csv(input_dir + "/ChemMNedges.tsv", sep="\t")
        node.to_csv(input_dir + "/ChemMNnodes.csv", index=False)

        newdf = df_edge
        newdf["StartAtt"] = np.nan
        newdf["EndAtt"] = np.nan
        for i, row in newdf.iterrows():
            for j, row in df.iterrows():
                if newdf["Start"][i] == df["IUPAC"][j]:
                    newdf.loc[i, "StartAtt"] = df["class"][j]
                if newdf["End"][i] == df["IUPAC"][j]:
                    newdf.loc[i, "EndAtt"] = df["class"][j]
        newdf.to_csv(input_dir + "/ChemMNcys.tsv", sep="\t")

        return newdf


def gnpsMNvsgnpsMAW(input_dir):
    
    """gnpsMNvsgnpsMAW checks with tanimoto similarity score, whether
    results from MAW GNPS and GNPS MN Masst results give same candidate

    Parameters:
    input_dir = input directory where you have stored the cytoscape file
    from GNPS MN results and have exported edge and node tables from cytoscape
    These two csv egde and node files must have "edge" and "node" in their name

    Returns:
    GNPS results with cluster index named
    GNPS MN results with a confirmation column if MAW detected same candidate,
    file named:

    Usage:
    gnpsMNvsgnpsMAW(input_dir)

    """
    # extract files with edges from MN results
    GMNfile_edge = [f for f in os.listdir(input_dir) if "edge" in f]
    # extract files with nodes from MN results
    GMNfile_node = [f for f in os.listdir(input_dir) if "node" in f]
    # read the files
    GMNdf_node = pd.read_csv(GMNfile_node[0])
    GMNdf_edge = pd.read_csv(GMNfile_edge[0])

    # extract only important columns from both csv files
    GMNdf_node = GMNdf_node[
        [
            "precursor mass",
            "RTMean",
            "UniqueFileSources",
            "charge",
            "cluster index",
            "componentindex",
            "Compound_Name",
            "Smiles",
            "SpectrumID",
        ]
    ]
    GMNdf_edge = GMNdf_edge[
        ["cosine_score", "EdgeAnnotation", "node1", "node2", "mass_difference"]
    ]

    # rename node1 to cluster index to merge nodes and edges results from MN
    GMNdf_edge = GMNdf_edge.rename(columns={"node1": "cluster index"})
    GMNdf = pd.merge(GMNdf_node, GMNdf_edge, on="cluster index")

    # Read results obtained from scoring_spec, named input_dir/MetabolomicsResults/scoredSpecDB.csv
    SDB = pd.read_csv(input_dir + "/MetabolomicsResults/scoredSpecDB.csv")
    # only keep GNPS resulst and remove other columns
    only_GNPS = SDB[SDB["annotation"].str.contains("GNPS")]
    only_GNPS = only_GNPS[
        [
            "id_X",
            "premz_x",
            "rtmean_x",
            "GNPSmax_similarity",
            "GNPSSMILES",
            #"GNPSspectrumID",
            "GNPScompound_name",
            "GNPSmirrorSpec",
        ]
    ]

    # from GNPS MAW results and GNPS MN results, calculate how many MAW results are same as MN:
    for i, row in only_GNPS.iterrows():
        for j, row in GMNdf.iterrows():
            if not isNaN(only_GNPS["GNPSSMILES"][i]) and not isNaN(GMNdf["Smiles"][j]):
                SKms = [
                    Chem.MolFromSmiles(only_GNPS["GNPSSMILES"][i]),
                    Chem.MolFromSmiles(GMNdf["Smiles"][j]),
                ]
                SKfps = [
                    AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=2048)
                    for x in SKms
                ]
                SKtn = DataStructs.FingerprintSimilarity(SKfps[0], SKfps[1])
                if SKtn == 1.0:
                    GMNdf.loc[j, "gnps_maw"] = "confirmed"
                    only_GNPS.loc[i, "index_MN_nodes"] = j
                elif SKtn < 1.0 and SKtn < 0.75:
                    GMNdf.loc[j, "gnps_maw"] = "similar"
                    only_GNPS.loc[i, "index_MN_nodes"] = j
    only_GNPS.to_csv(input_dir + "/only_GNPS.csv")
    GMNdf.to_csv(input_dir + "/GMNdf.csv")

