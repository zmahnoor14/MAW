#!/usr/bin/env python
# coding: utf-8

import yaml
import time
import provenance as p
from yaml.loader import SafeLoader

basic_config = {
    "blobstores": {
        "disk": {
            "type": "disk",
            "cachedir": "provenance-intro-artifacts",
            "read": True,
            "write": True,
            "delete": True,
        }
    },
    "artifact_repos": {
        "local": {
            "type": "memory"
        }
    },
    "default_repo": "local"
}
p.load_config(basic_config)

import glob
import json
import os
import re
import time
import wget
import urllib.parse
import argparse


import numpy as np
import pandas as pd
import pubchempy as pcp


from pybatchclassyfire import *
from pandas import json_normalize
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem import PandasTools

import plotly.express as px

# def argparser()

def isNaN(string):
    return string != string

# # SpecDB Post Processing
#@p.provenance()
def spec_postproc(entry, Source="all"):
    # currently only these subsets are removed from the names from GNPS
    matches = [
        "M+",
        "[M",
        "M-",
        "2M",
        "M*", 
        "20.0",
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
        "CollisionEnergy:"
    ]


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
    # # in case if we need HMDB later
    # if os.path.exists(entry.split("/DS")[0] + "/hmdb_dframe_str.csv"):
    #     extract_smiles = pd.read_csv(entry.split("/DS")[0] + "/hmdb_dframe_str.csv", low_memory=False)
   

    msp_file = glob.glob(
        entry + "/spectral_dereplication" + "/*.csv"
    )
    


    if len(msp_file) > 0:

        if os.path.exists(msp_file[0]):

            msp = pd.read_csv(msp_file[0])
            msp["mbank_results_csv"] = np.nan
            msp["gnps_results_csv"] = np.nan
            msp["hmdb_results_csv"] = np.nan
            # enter the directory with /spectral_dereplication/ results

            # enter the directory with /spectral_dereplication/ results
            # GNPS Results
            if Source == "gnps" or Source == "all":
                msp["gnps_results_csv"] = np.nan

                # print(entry)
                # enter the directory with /spectral_dereplication/ results
                sub_dir = (
                    entry + "/spectral_dereplication/GNPS/"
                )

                if os.path.exists(sub_dir):
                    files = glob.glob(sub_dir + "/*.csv")
                    # print(files)
                    files = [item for item in files if 'proc' not in item]

                    for mz, row in msp.iterrows():
                        for fls_g in files:

                            if msp["id_X"][mz] in fls_g:
                                
                                gnps_df = pd.read_csv(fls_g)
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
                                    gnps_df = gnps_df.drop_duplicates(
                                        subset=["GNPSSMILES"]
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

                                msp.loc[
                                    mz, "gnps_results_csv"
                                ] = csvname
                                
                                if not os.path.exists(csvname):
                                    #print("this is wrong?")
                                    #print(csvname)
                                    #print(os.path.splitext(fls_g)[0])
                                    gnps_df.to_csv(csvname)


            msp.to_csv(msp_file[0])
            # HMDB Results
            if Source == "hmdb" or Source == "all":
                sub_dir = (
                    entry + "/spectral_dereplication/HMDB/"
                )
                if os.path.exists(sub_dir):
                    files = glob.glob(sub_dir + "/*.csv")
                    files = [item for item in files if 'proc' not in item]
                    if os.path.exists(sub_dir):
                        # print(files)
                        for mz, row in msp.iterrows():
                            #print(mz)
                            # print(msp["id_X"][mz])
                            for fls_h in files:
                                 if msp["id_X"][mz] in fls_h:
                                        hmdb_df = pd.read_csv(fls_h)

                                        if len(hmdb_df) > 0:
                                            if "HMDBSMILES" in hmdb_df.columns:
                                                #print(hmdb_df)
                                                for i, row in hmdb_df.iterrows():
                                                # if compound name is present
                                                    if not HMDB_Scoring(hmdb_df, i):
                                                        hmdb_df.drop(i, inplace=True)
                                                hmdb_df = hmdb_df.drop_duplicates(
                                                    subset=["HMDBSMILES"]
                                                )


                                                csvname = (
                                                    (os.path.splitext(fls_h)[0])
                                                    + "proc"
                                                    + ".csv"
                                                )  
                                                msp.loc[
                                                    mz, "hmdb_results_csv"
                                                ] = csvname

                                                if not os.path.exists(csvname):
                                                    hmdb_df.to_csv(csvname) 
                                            else:    
                                                # merge on basis of id, frame and hmdb result files
                                                SmilesHM = pd.merge(
                                                    hmdb_df,
                                                    extract_smiles,
                                                    left_on=hmdb_df.HMDBcompoundID,
                                                    right_on=extract_smiles.DATABASE_ID,
                                                )
                                                hmdb_df["HMDBcompoundID"] = np.nan
                                                hmdb_df["HMDBSMILES"] = np.nan
                                                hmdb_df["HMDBformula"] = np.nan
                                                hmdb_df["HMDBcompound_name"] = np.nan
                                                for i, row in hmdb_df.iterrows():
                                                    #print(i)
                                                    # if compound name is present
                                                    if HMDB_Scoring(hmdb_df, i):

                                                        #hmdb_df.drop(i, inplace=True)
                                                        for j, row in SmilesHM.iterrows():
                                                            #print("SmilesHM")
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

                                                                ]
                                                                #print(hmdb_df["HMDBSMILES"][i])
                #                             
                                            #print(hmdb_df)
#                                         hmdb_df = hmdb_df.drop_duplicates(
#                                              subset=["HMDBSMILES"]
#                                          )
                                        csvname = (
                                            (os.path.splitext(fls_h)[0])
                                            + "proc"
                                            + ".csv"
                                        )  
                                        msp.loc[
                                            mz, "hmdb_results_csv"
                                        ] = csvname

                                        if not os.path.exists(csvname):
                                            hmdb_df.to_csv(csvname) 
            msp.to_csv(msp_file[0])
            # MASSBANK Results

            # enter the directory with /spectral_dereplication/ results
            if Source == "mbank" or Source == "all":

                sub_dir = (
                    
                    entry
                    + "/spectral_dereplication/MassBank/"
                )
                if os.path.exists(sub_dir):
                    files = glob.glob(sub_dir + "/*.csv")
                    files = [item for item in files if 'proc' not in item]
                    for mz, row in msp.iterrows():
                        # print(msp["id_X"][mz])
                        for fls_m in files:
                            if msp["id_X"][mz] in fls_m:
                                mbank_df = pd.read_csv(fls_m)
                                if len(mbank_df) > 0:

                                    for i, row in mbank_df.iterrows():
                                        # if compound name is present
                                         if not MB_Scoring(mbank_df, i):
                                            mbank_df.drop(i, inplace=True)
                                mbank_df = mbank_df.drop_duplicates(
                                    subset=["MBSMILES"]
                                )
                                csvname = (
                                    (os.path.splitext(fls_m)[0])
                                    + "proc"
                                    + ".csv"
                                )

                                msp.loc[
                                    mz, "mbank_results_csv"
                                ] = csvname

                                if not os.path.exists(csvname):

                                    mbank_df.to_csv(csvname)
            msp.to_csv(msp_file[0])
  
#@p.provenance()
def MCSS_for_SpecDB(entry, Source = "all"):
    # Describe the heavy atoms to be considered for MCSS
    heavy_atoms = ["C", "N", "P", "O", "S"]
    
    specdb_msp_file = glob.glob(entry + "/spectral_dereplication" + "/*.csv"
    )

    if len(specdb_msp_file) > 0:

        if os.path.exists(specdb_msp_file[0]):

            spec_msp = pd.read_csv(specdb_msp_file[0])

            for mz, row in spec_msp.iterrows():

                if Source == "gnps" or Source == "specdb" or Source == "all":

                    sub_dir = (
                        entry
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
                        entry
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
                        entry
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

def checkSMILES_validity(resultcsv):
    

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

    return results
 
def one_candidate_selection_metfrag(
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
    2. SGM (Metfrag, GNPS, MassBank)
    3. SHM (Metfrag, HMDB, MassBank)
    4. SGH (Metfrag, GNPS, HMDB)
    5. GHM (GNPS, HMDB, MassBank)
    6. SG (Metfrag, GNPS)
    7. SH (Metfrag, HMDB)
    8. SM (Metfrag, MassBank)
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

    df["MetFrag"] = np.nan
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
                        Chem.MolFromSmiles(sirius_df["SMILES"][sirius_i]),
                    ]
                    fps = [
                        AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=2048)
                        for x in ms
                    ]
                    tn = DataStructs.FingerprintSimilarity(fps[0], fps[1])
                    # since we are dealing with idenity here so tanimoto of 0.99 is appropriate
                    if tn >= tn_ident:

                        # if SIRIUS is blank, add the SIRIUS id
                        if isNaN(df["MetFrag"][smiles]):#changed

                            df.loc[smiles, "MetFrag"] = sirius_df["rank_ids"][sirius_i]#changed
                        # if not empty, add SIRIUS id, with a comma
                        else:
                            df.loc[smiles, "MetFrag"] = (
                                str(df["MetFrag"][smiles])
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

def add_count_column_metfrag(df_one_candidate):
    df_one_candidate = df_one_candidate.dropna(axis=0, how="all", subset = ['MetFrag', 'GNPS', 'MassBank', 'HMDB'])
    # create new df only with the Sources column
    df = pd.DataFrame(
        {
            "MetFrag": df_one_candidate["MetFrag"],
            "GNPS": df_one_candidate["GNPS"],
            "MassBank": df_one_candidate["MassBank"],
            "HMDB": df_one_candidate["HMDB"],
        }
    )

    # df_one_candidate = df_one_candidate.dropna(subset=["SIRIUS", "GNPS", "HMDB", "MassBank"], how='all', inplace=True)

    index_SIRIUS = [x for x, row in df.iterrows() if not isNaN(df["MetFrag"][x])]
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
            if df_one_candidate["Source"][r]=="MetFrag":
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
            if sorted_count_one_candidate["Source"][r]=="MetFrag":
                sorted_count_one_candidate.loc[r, "rank_db"] = 2
            if sorted_count_one_candidate["Source"][r]=="MassBank":
                sorted_count_one_candidate.loc[r, "rank_db"] = 3
            if sorted_count_one_candidate["Source"][r]=="HMDB":
                sorted_count_one_candidate.loc[r, "rank_db"] = 4

        sorted_count_one_candidate.sort_values(
            by=["Count", "rank_num", "rank_db"], ascending=[False, True, True],
            inplace=True)
        return sorted_count_one_candidate

def sources_1_metfrag(candidates_with_counts, merged_df, mer, sirius_df):
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


    df_count_1 = df_count_1[df_count_1["rank_num"] == min(df_count_1["rank_num"])]

    df_count_1["count_min"] = [
        str(df_count_1["MetFrag"][x])
        + str(df_count_1["GNPS"][x])
        + str(df_count_1["MassBank"][x])
        + str(df_count_1["HMDB"][x])
        for x, row in df_count_1.iterrows()
    ]

    df_count_1["count_max"] = [x.count("_") for x in df_count_1["count_min"]]

    df_count_1 = df_count_1.sort_values(by="count_max", ascending=False)

    df_count_1.reset_index(drop=True, inplace=True)

    merged_df.loc[mer, "AnnotationCount"] = df_count_1["Count"][0]

    gnps_indices = list(df_count_1[(df_count_1["GNPS"].notnull())].index)
    mbank_indices = list(df_count_1[(df_count_1["MassBank"].notnull())].index)
    hmdb_indices = list(df_count_1[(df_count_1["HMDB"].notnull())].index)
    sirius_indices = list(df_count_1[(df_count_1["MetFrag"].notnull())].index)

    if 0 in sirius_indices:
        merged_df.loc[mer, "AnnotationSources"] = (
            str(merged_df["AnnotationSources"][mer]) + "|MetFrag"
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
    if df_count_1["SMILES"][0] != "CC[S](=O)(O)O":
        try:
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
                if "MetFrag" == merged_df["AnnotationSources"][mer]:

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
        except:
            pass
    if "MetFrag" not in merged_df["AnnotationSources"][mer]:
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
    elif "MetFrag" == merged_df["AnnotationSources"][mer]:
        merged_df.loc[mer, "MSILevel"] = 3
            
    return merged_df

def sources_2_metfrag(candidates_with_counts, merged_df, mer, sirius_df):

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

    df_countnew = df_count_2[df_count_2["rank_num"] == min(df_count_2["rank_num"])]
    df_countnew["count_min"] = [
        str(df_countnew["MetFrag"][x])
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
    sirius_indices = list(df_countnew[(df_countnew["MetFrag"].notnull())].index)

    if 0 in sirius_indices:
        # print("sirius")
        merged_df.loc[mer, "AnnotationSources"] = (
            str(merged_df["AnnotationSources"][mer]) + "|MetFrag"
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

    merged_df.loc[mer, "MSILevel"] = 2


    merged_df["AnnotationSources"][mer] = merged_df["AnnotationSources"][mer].replace(
        "nan|", ""
    )
    merged_df.loc[mer, "AnnotationCount"] = df_countnew["Count"][0]
    
    merged_df.loc[mer, "SMILES"] = df_countnew["SMILES"][0]
    try:
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
    except:
        pass
    if "MetFrag" not in merged_df["AnnotationSources"][mer]:
        merged_df["superclass"][mer] = np.nan
        merged_df["class"][mer] = np.nan
        merged_df["subclass"][mer] = np.nan
        merged_df["ClassificationSource"][mer] = np.nan
    
    
    return merged_df

def sources_3_metfrag(candidates_with_counts, merged_df, mer, sirius_df):

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
        str(df_count_3["MetFrag"][x])
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
    sirius_indices = list(df_count_3[(df_count_3["MetFrag"].notnull())].index)

    if 0 in sirius_indices:
        # print("sirius")
        merged_df.loc[mer, "AnnotationSources"] = (
            str(merged_df["AnnotationSources"][mer]) + "|MetFrag"
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

    merged_df.loc[mer, "MSILevel"] = 2
    merged_df.loc[mer, "SMILES"] = df_count_3["SMILES"][0]
    
    try: 
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
    except:
        pass
    if "MetFrag" not in merged_df["AnnotationSources"][mer]:
        merged_df["superclass"][mer] = np.nan
        merged_df["class"][mer] = np.nan
        merged_df["subclass"][mer] = np.nan
        merged_df["ClassificationSource"][mer] = np.nan
       

    
    return merged_df

def sources_4_metfrag(candidates_with_counts, merged_df, mer, sirius_df):

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
        str(df_count_4["MetFrag"][x])
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
    sirius_indices = list(df_count_4[(df_count_4["MetFrag"].notnull())].index)

    if 0 in sirius_indices:
        # print("sirius")
        merged_df.loc[mer, "AnnotationSources"] = (
            str(merged_df["AnnotationSources"][mer]) + "|MetFrag"
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
    
    merged_df.loc[mer, "MSILevel"] = 2
    
    merged_df.loc[mer, "SMILES"] = df_count_4["SMILES"][0]
    try:
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
    except:
        pass
    if "MetFrag" not in merged_df["AnnotationSources"][mer]:
        merged_df["superclass"][mer] = np.nan
        merged_df["class"][mer] = np.nan
        merged_df["subclass"][mer] = np.nan
        merged_df["ClassificationSource"][mer] = np.nan  
     
    return merged_df

#@p.provenance()
def CandidateSelection_SimilarityandIdentity_Metfrag(entry, standards = False):

    sub_dir_spec = entry + "/spectral_dereplication/"
    sub_dir_sir = entry + "/insilico/MetFrag/" # change
    if os.path.exists(sub_dir_spec) and os.path.exists(sub_dir_sir):
        spec_msp_csv = glob.glob(entry + "/spectral_dereplication" + "/*.csv"
        )
        sir_msp_csv =  entry + "/insilico/MS1DATA.csv"
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
                    "gnps_results_csv",
                    "hmdb_results_csv"
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
                    "MetFragCSV", # change
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
            
            #for_only_formula = [] # change
            #for_formula_canopus = [] # change

            can_selec_dir = entry + "/Candidate_Selection"
            if not os.path.isdir(can_selec_dir):
                os.mkdir(can_selec_dir)
            for mer, rows in merged_df.iterrows():
                print(mer)
                # change starts here
    #             if not isNaN(merged_df["sirius_result_dir"][mer]):
    #                 sirius_csv = merged_df["sirius_result_dir"][mer]
    #                 # print(sirius_csv)
    #             else:
    #                 df = pd.DataFrame(list())
    #                 # print(df)
    #                 df.to_csv("./empty_csv.csv")
    #                 sirius_csv = "./empty_csv.csv"
                metfrag_csv = merged_df["MetFragCSV"][mer] # change 
                mbank_csv = merged_df["mbank_results_csv"][mer]
                gnps_csv = merged_df["gnps_results_csv"][mer]
                hmdb_csv = merged_df["hmdb_results_csv"][mer]
                if (
                    os.path.exists(metfrag_csv)
                    and os.path.exists(gnps_csv)
                    and os.path.exists(mbank_csv)
                    and os.path.exists(hmdb_csv)
                ):
                    metfrag_df = pd.read_csv(metfrag_csv) # change 
                    # here I change Metfrag to SIRIUS just for convinience as a lot of code is develped based on SIRIUS
                    # but same can be used for MetFrag
                    sirius_df= metfrag_df
                    if len(metfrag_df) > 0:
                        sirius_df = sirius_df.drop_duplicates("SMILES")
                        sirius_df = sirius_df.dropna(subset=["SMILES"])
                    if len(metfrag_df) == 0:
                        print("not gonna work")
                        
                    mbank_df = pd.read_csv(mbank_csv)
                    if len(mbank_df) > 0:
                        mbank_df = mbank_df.drop_duplicates("MBSMILES")
                        mbank_df = mbank_df.dropna(subset=["MBSMILES"])

                    gnps_df = pd.read_csv(gnps_csv)
                    if len(gnps_df) > 0:
                        gnps_df = gnps_df.drop_duplicates("GNPSSMILES")
                        gnps_df = gnps_df.dropna(subset=["GNPSSMILES"])

                    hmdb_df = pd.read_csv(hmdb_csv)
                    # print(hmdb_df)
                    if len(hmdb_df) > 0:
                        hmdb_df = hmdb_df.drop_duplicates("HMDBSMILES")
                        hmdb_df = hmdb_df.dropna(subset=["HMDBSMILES"])
                
                if (
                    os.path.exists(metfrag_csv)
                    and os.path.exists(gnps_csv)
                    and os.path.exists(mbank_csv)
                    and os.path.exists(hmdb_csv)
                ):
                    metfrag_df = pd.read_csv(metfrag_csv) # change 
                    # here I change Metfrag to SIRIUS just for convinience as a lot of code is develped based on SIRIUS
                    # but same can be used for MetFrag
                    sirius_df= metfrag_df
                    if len(metfrag_df) > 0:
                        sirius_df = sirius_df.drop_duplicates("SMILES")
                        sirius_df = sirius_df.dropna(subset=["SMILES"])
                        
                    mbank_df = pd.read_csv(mbank_csv)
                    if len(mbank_df) > 0:
                        mbank_df = mbank_df.drop_duplicates("MBSMILES")
                        mbank_df = mbank_df.dropna(subset=["MBSMILES"])

                    gnps_df = pd.read_csv(gnps_csv)
                    if len(gnps_df) > 0:
                        gnps_df = gnps_df.drop_duplicates("GNPSSMILES")
                        gnps_df = gnps_df.dropna(subset=["GNPSSMILES"])

                    hmdb_df = pd.read_csv(hmdb_csv)
                    # print(hmdb_df)
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
                        mbank_df["rank_ids"] = ["M_" + str(s + 1) for s in range(len(mbank_df))]

                        gnps_df["rank_ids"] = ["G_" + str(s + 1) for s in range(len(gnps_df))]

                        hmdb_df["rank_ids"] = ["H_" + str(s + 1) for s in range(len(hmdb_df))]
                        
                        sirius_df["rank_ids"] = ["E_" + str(s) for s in range(len(sirius_df))]
                        sirius_df["Source"] = "MetFrag"

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
                            *(list(sirius_df["SMILES"])),
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
                                entry
                                + "/"
                                + "Candidate_Selection"
                                + "/"
                                + str(merged_df["premz"][mer])
                                + "_ChemMNedges.tsv",
                                sep="\t",
                            )
                        one_candidate = one_candidate_selection_metfrag(
                            sm,
                            sirius_df=sirius_df,
                            mbank_df=mbank_df,
                            gnps_df=gnps_df,
                            hmdb_df=hmdb_df,
                            Source="SGHM",
                        )

                        candidates_with_counts = add_count_column_metfrag(one_candidate)
                        candidates_with_counts.to_csv(
                            entry
                            + "/"
                            + "Candidate_Selection"
                            + "/"
                            + str(merged_df["premz"][mer])
                            + "sorted_candidate_list.tsv",
                            sep="\t",
                        )
                        if max(candidates_with_counts["Count"]) == 4:
                            sources_4_metfrag(candidates_with_counts, merged_df, mer, sirius_df)
                        if max(candidates_with_counts["Count"]) == 3:
                            sources_3_metfrag(candidates_with_counts, merged_df, mer, sirius_df)
                        if max(candidates_with_counts["Count"]) == 2:
                            sources_2_metfrag(candidates_with_counts, merged_df, mer, sirius_df)
                        if max(candidates_with_counts["Count"]) == 1:
                            sources_1_metfrag(candidates_with_counts, merged_df, mer, sirius_df)
                    # 2 SGM
                    elif (
                        len(sirius_df) > 0
                        and len(gnps_df) > 0
                        and len(mbank_df) > 0
                        and len(hmdb_df) == 0
                    ):
                        mbank_df["rank_ids"] = ["M_" + str(s + 1) for s in range(len(mbank_df))]

                        gnps_df["rank_ids"] = ["G_" + str(s + 1) for s in range(len(gnps_df))]

                        sirius_df["rank_ids"] = ["E_" + str(s) for s in range(len(sirius_df))]
                        sirius_df["Source"] = "MetFrag"

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
                            *(list(sirius_df["SMILES"])),
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
                                entry
                                + "/"
                                + "Candidate_Selection"
                                + "/"
                                + str(merged_df["premz"][mer])
                                + "_ChemMNedges.tsv",
                                sep="\t",
                            )
                        one_candidate = one_candidate_selection_metfrag(
                            sm,
                            sirius_df=sirius_df,
                            mbank_df=mbank_df,
                            gnps_df=gnps_df,
                            # hmdb_df = hmdb_df,
                            Source="SGM",
                        )
                        candidates_with_counts = add_count_column_metfrag(one_candidate)
                        candidates_with_counts.to_csv(
                            entry
                            + "/"
                            + "Candidate_Selection"
                            + "/"
                            + str(merged_df["premz"][mer])
                            + "sorted_candidate_list.tsv",
                            sep="\t",
                        )

                        if max(candidates_with_counts["Count"]) == 3:
                            sources_3_metfrag(candidates_with_counts, merged_df, mer, sirius_df)
                        if max(candidates_with_counts["Count"]) == 2:
                            sources_2_metfrag(candidates_with_counts, merged_df, mer, sirius_df)
                        if max(candidates_with_counts["Count"]) == 1:
                            sources_1_metfrag(candidates_with_counts, merged_df, mer, sirius_df)
                    # 3 SHM
                    elif (
                        len(sirius_df) > 0
                        and len(gnps_df) == 0
                        and len(mbank_df) > 0
                        and len(hmdb_df) > 0
                    ):

                        mbank_df["rank_ids"] = ["M_" + str(s + 1) for s in range(len(mbank_df))]

                        # gnps_df["rank_ids"] = ["G_" + str(s+1) for s in range(len(gnps_df))]

                        hmdb_df["rank_ids"] = ["H_" + str(s + 1) for s in range(len(hmdb_df))]

                        sirius_df["rank_ids"] = ["E_" + str(s) for s in range(len(sirius_df))]
                        sirius_df["Source"] = "MetFrag"

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
                            *(list(sirius_df["SMILES"])),
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
                                entry
                                + "/"
                                + "Candidate_Selection"
                                + "/"
                                + str(merged_df["premz"][mer])
                                + "_ChemMNedges.tsv",
                                sep="\t",
                            )
                        one_candidate = one_candidate_selection_metfrag(
                            sm,
                            sirius_df=sirius_df,
                            mbank_df=mbank_df,
                            # gnps_df = gnps_df ,
                            hmdb_df=hmdb_df,
                            Source="SHM",
                        )
                        candidates_with_counts = add_count_column_metfrag(one_candidate)
                        candidates_with_counts.to_csv(
                            entry
                            + "/"
                            + "Candidate_Selection"
                            + "/"
                            + str(merged_df["premz"][mer])
                            + "sorted_candidate_list.tsv",
                            sep="\t",
                        )

                        if max(candidates_with_counts["Count"]) == 3:
                            sources_3_metfrag(candidates_with_counts, merged_df, mer, sirius_df)
                        if max(candidates_with_counts["Count"]) == 2:
                            sources_2_metfrag(candidates_with_counts, merged_df, mer, sirius_df)
                        if max(candidates_with_counts["Count"]) == 1:
                            sources_1_metfrag(candidates_with_counts, merged_df, mer, sirius_df)
                    # 4 SGH
                    elif (
                        len(sirius_df) > 0
                        and len(gnps_df) > 0
                        and len(mbank_df) == 0
                        and len(hmdb_df) > 0
                    ):
                        # mbank_df["rank_ids"] = ["M_" + str(s+1) for s in range(len(mbank_df))]

                        gnps_df["rank_ids"] = ["G_" + str(s + 1) for s in range(len(gnps_df))]

                        hmdb_df["rank_ids"] = ["H_" + str(s + 1) for s in range(len(hmdb_df))]
                        
                        sirius_df["rank_ids"] = ["E_" + str(s) for s in range(len(sirius_df))]
                        sirius_df["Source"] = "MetFrag"


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
                            *(list(sirius_df["SMILES"])),
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
                                entry
                                + "/"
                                + "Candidate_Selection"
                                + "/"
                                + str(merged_df["premz"][mer])
                                + "_ChemMNedges.tsv",
                                sep="\t",
                            )
                        one_candidate = one_candidate_selection_metfrag(
                            sm,
                            sirius_df=sirius_df,
                            # mbank_df = mbank_df,
                            gnps_df=gnps_df,
                            hmdb_df=hmdb_df,
                            Source="SGH",
                        )
                        candidates_with_counts = add_count_column_metfrag(one_candidate)
                        candidates_with_counts.to_csv(
                            entry
                            + "/"
                            + "Candidate_Selection"
                            + "/"
                            + str(merged_df["premz"][mer])
                            + "sorted_candidate_list.tsv",
                            sep="\t",
                        )
                        if max(candidates_with_counts["Count"]) == 3:
                            sources_3_metfrag(candidates_with_counts, merged_df, mer, sirius_df)
                        if max(candidates_with_counts["Count"]) == 2:
                            sources_2_metfrag(candidates_with_counts, merged_df, mer, sirius_df)
                        if max(candidates_with_counts["Count"]) == 1:
                            sources_1_metfrag(candidates_with_counts, merged_df, mer, sirius_df)
                    # 5 GHM
                    elif (
                        len(sirius_df) == 0
                        and len(gnps_df) > 0
                        and len(mbank_df) > 0
                        and len(hmdb_df) > 0
                    ):
                        mbank_df["rank_ids"] = ["M_" + str(s + 1) for s in range(len(mbank_df))]

                        gnps_df["rank_ids"] = ["G_" + str(s + 1) for s in range(len(gnps_df))]

                        hmdb_df["rank_ids"] = ["H_" + str(s + 1) for s in range(len(hmdb_df))]

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
                                entry
                                + "/"
                                + "Candidate_Selection"
                                + "/"
                                + str(merged_df["premz"][mer])
                                + "_ChemMNedges.tsv",
                                sep="\t",
                            )
                        one_candidate = one_candidate_selection_metfrag(
                            sm,
                            # sirius_df = sirius_df,
                            mbank_df=mbank_df,
                            gnps_df=gnps_df,
                            hmdb_df=hmdb_df,
                            Source="GHM",
                        )
                        candidates_with_counts = add_count_column_metfrag(one_candidate)
                        candidates_with_counts.to_csv(
                            entry
                            + "/"
                            + "Candidate_Selection"
                            + "/"
                            + str(merged_df["premz"][mer])
                            + "sorted_candidate_list.tsv",
                            sep="\t",
                        )

                        if max(candidates_with_counts["Count"]) == 3:
                            sources_3_metfrag(candidates_with_counts, merged_df, mer, sirius_df=None)
                        if max(candidates_with_counts["Count"]) == 2:
                            sources_2_metfrag(candidates_with_counts, merged_df, mer, sirius_df=None)
                        if max(candidates_with_counts["Count"]) == 1:
                            sources_1_metfrag(candidates_with_counts, merged_df, mer, sirius_df=None)
                    # 6 SG
                    elif (
                        len(sirius_df) > 0
                        and len(gnps_df) > 0
                        and len(mbank_df) == 0
                        and len(hmdb_df) == 0
                    ):
                        # mbank_df["rank_ids"] = ["M_" + str(s+1) for s in range(len(mbank_df))]

                        gnps_df["rank_ids"] = ["G_" + str(s + 1) for s in range(len(gnps_df))]

                        sirius_df["rank_ids"] = ["E_" + str(s) for s in range(len(sirius_df))]
                        sirius_df["Source"] = "MetFrag"

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
                            *(list(sirius_df["SMILES"])),
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
                                entry
                                + "/"
                                + "Candidate_Selection"
                                + "/"
                                + str(merged_df["premz"][mer])
                                + "_ChemMNedges.tsv",
                                sep="\t",
                            )

                        one_candidate = one_candidate_selection_metfrag(
                            sm,
                            sirius_df=sirius_df,
                            # mbank_df = mbank_df,
                            gnps_df=gnps_df,
                            # hmdb_df = hmdb_df,
                            Source="SG",
                        )
                        candidates_with_counts = add_count_column_metfrag(one_candidate)
                        candidates_with_counts.to_csv(
                            entry
                            + "/"
                            + "Candidate_Selection"
                            + "/"
                            + str(merged_df["premz"][mer])
                            + "sorted_candidate_list.tsv",
                            sep="\t",
                        )

                        if max(candidates_with_counts["Count"]) == 2:
                            sources_2_metfrag(candidates_with_counts, merged_df, mer, sirius_df)
                        if max(candidates_with_counts["Count"]) == 1:
                            sources_1_metfrag(candidates_with_counts, merged_df, mer, sirius_df)
                    # 7 SH
                    elif (
                        len(sirius_df) > 0
                        and len(gnps_df) == 0
                        and len(mbank_df) == 0
                        and len(hmdb_df) > 0
                    ):


                        hmdb_df["rank_ids"] = ["H_" + str(s + 1) for s in range(len(hmdb_df))]

                        sirius_df["rank_ids"] = ["E_" + str(s) for s in range(len(sirius_df))]
                        sirius_df["Source"] = "MetFrag"

                        source_l1 = [
                            *(list(sirius_df["Source"])),
                            *(list(hmdb_df["Source"])),
                        ]

                        rank_l2 = [
                            *(list(sirius_df["rank_ids"])),
                            *(list(hmdb_df["rank_ids"])),
                        ]

                        smiles_l3 = [
                            *(list(sirius_df["SMILES"])),
                            *(list(hmdb_df["HMDBSMILES"])),
                        ]

                        sm = pd.DataFrame(
                            list(zip(source_l1, rank_l2, smiles_l3)),
                            columns=["Source", "ranks", "SMILES"],
                        )
                        df_edge = chemMN_CandidateSelection(sm)
                        if df_edge is not None:

                            df_edge.to_csv(
                                entry
                                + "/"
                                + "Candidate_Selection"
                                + "/"
                                + str(merged_df["premz"][mer])
                                + "_ChemMNedges.tsv",
                                sep="\t",
                            )

                        one_candidate = one_candidate_selection_metfrag(
                            sm,
                            sirius_df=sirius_df,
                            # mbank_df = mbank_df,
                            # gnps_df = gnps_df ,
                            hmdb_df=hmdb_df,
                            Source="SH",
                        )
                        candidates_with_counts = add_count_column_metfrag(one_candidate)
                        candidates_with_counts.to_csv(
                            entry
                            + "/"
                            + "Candidate_Selection"
                            + "/"
                            + str(merged_df["premz"][mer])
                            + "sorted_candidate_list.tsv",
                            sep="\t",
                        )

                        if max(candidates_with_counts["Count"]) == 2:
                            sources_2_metfrag(candidates_with_counts, merged_df, mer, sirius_df)
                        if max(candidates_with_counts["Count"]) == 1:
                            sources_1_metfrag(candidates_with_counts, merged_df, mer, sirius_df)
                    # 8 SM
                    elif (
                        len(sirius_df) > 0
                        and len(gnps_df) == 0
                        and len(mbank_df) > 0
                        and len(hmdb_df) == 0
                    ):
                        mbank_df["rank_ids"] = ["M_" + str(s + 1) for s in range(len(mbank_df))]

                        sirius_df["rank_ids"] = ["E_" + str(s) for s in range(len(sirius_df))]
                        sirius_df["Source"] = "MetFrag"

                        source_l1 = [
                            *(list(sirius_df["Source"])),
                            *(list(mbank_df["Source"])),
                        ]

                        rank_l2 = [
                            *(list(sirius_df["rank_ids"])),
                            *(list(mbank_df["rank_ids"])),
                        ]

                        smiles_l3 = [
                            *(list(sirius_df["SMILES"])),
                            *(list(mbank_df["MBSMILES"])),
                        ]

                        sm = pd.DataFrame(
                            list(zip(source_l1, rank_l2, smiles_l3)),
                            columns=["Source", "ranks", "SMILES"],
                        )

                        df_edge = chemMN_CandidateSelection(sm)
                        if df_edge is not None:

                            df_edge.to_csv(
                                entry
                                + "/"
                                + "Candidate_Selection"
                                + "/"
                                + str(merged_df["premz"][mer])
                                + "_ChemMNedges.tsv",
                                sep="\t",
                            )

                        one_candidate = one_candidate_selection_metfrag(
                            sm,
                            sirius_df=sirius_df,
                            mbank_df=mbank_df,
                            # gnps_df = gnps_df ,
                            # hmdb_df = hmdb_df,
                            Source="SM",
                        )
                        candidates_with_counts = add_count_column_metfrag(one_candidate)
                        candidates_with_counts.to_csv(
                            entry
                            + "/"
                            + "Candidate_Selection"
                            + "/"
                            + str(merged_df["premz"][mer])
                            + "sorted_candidate_list.tsv",
                            sep="\t",
                        )
                        if max(candidates_with_counts["Count"]) == 2:
                            sources_2_metfrag(candidates_with_counts, merged_df, mer, sirius_df)
                        if max(candidates_with_counts["Count"]) == 1:
                            sources_1_metfrag(candidates_with_counts, merged_df, mer, sirius_df)
                    # 9 GM
                    elif (
                        len(sirius_df) == 0
                        and len(gnps_df) > 0
                        and len(mbank_df) > 0
                        and len(hmdb_df) == 0
                    ):
                        mbank_df["rank_ids"] = ["M_" + str(s + 1) for s in range(len(mbank_df))]

                        gnps_df["rank_ids"] = ["G_" + str(s + 1) for s in range(len(gnps_df))]

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
                                entry
                                + "/"
                                + "Candidate_Selection"
                                + "/"
                                + str(merged_df["premz"][mer])
                                + "_ChemMNedges.tsv",
                                sep="\t",
                            )

                        one_candidate = one_candidate_selection_metfrag(
                            sm,
                            # sirius_df = sirius_df,
                            mbank_df=mbank_df,
                            gnps_df=gnps_df,
                            # hmdb_df = hmdb_df,
                            Source="GM",
                        )
                        candidates_with_counts = add_count_column_metfrag(one_candidate)
                        candidates_with_counts.to_csv(
                            entry
                            + "/"
                            + "Candidate_Selection"
                            + "/"
                            + str(merged_df["premz"][mer])
                            + "sorted_candidate_list.tsv",
                            sep="\t",
                        )
                        if max(candidates_with_counts["Count"]) == 2:
                            sources_2_metfrag(candidates_with_counts, merged_df, mer, sirius_df=None)
                        if max(candidates_with_counts["Count"]) == 1:
                            sources_1_metfrag(candidates_with_counts, merged_df, mer, sirius_df=None)
                    # 10 GH
                    elif (
                        len(sirius_df) == 0
                        and len(gnps_df) > 0
                        and len(mbank_df) == 0
                        and len(hmdb_df) > 0
                    ):
                        # mbank_df["rank_ids"] = ["M_" + str(s+1) for s in range(len(mbank_df))]

                        gnps_df["rank_ids"] = ["G_" + str(s + 1) for s in range(len(gnps_df))]

                        hmdb_df["rank_ids"] = ["H_" + str(s + 1) for s in range(len(hmdb_df))]

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
                                entry
                                + "/"
                                + "Candidate_Selection"
                                + "/"
                                + str(merged_df["premz"][mer])
                                + "_ChemMNedges.tsv",
                                sep="\t",
                            )
                        one_candidate = one_candidate_selection_metfrag(
                            sm,
                            # sirius_df = sirius_df,
                            # mbank_df = mbank_df,
                            gnps_df=gnps_df,
                            hmdb_df=hmdb_df,
                            Source="GH",
                        )
                        candidates_with_counts = add_count_column_metfrag(one_candidate)
                        candidates_with_counts.to_csv(
                            entry
                            + "/"
                            + "Candidate_Selection"
                            + "/"
                            + str(merged_df["premz"][mer])
                            + "sorted_candidate_list.tsv",
                            sep="\t",
                        )

                        if max(candidates_with_counts["Count"]) == 2:
                            sources_2_metfrag(candidates_with_counts, merged_df, mer, sirius_df=None)
                        if max(candidates_with_counts["Count"]) == 1:
                            sources_1_metfrag(candidates_with_counts, merged_df, mer, sirius_df=None)
                    # 11 HM
                    elif (
                        len(sirius_df) == 0
                        and len(gnps_df) == 0
                        and len(mbank_df) > 0
                        and len(hmdb_df) > 0
                    ):
                        mbank_df["rank_ids"] = ["M_" + str(s + 1) for s in range(len(mbank_df))]

                        # gnps_df["rank_ids"] = ["G_" + str(s+1) for s in range(len(gnps_df))]

                        hmdb_df["rank_ids"] = ["H_" + str(s + 1) for s in range(len(hmdb_df))]

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
                                entry
                                + "/"
                                + "Candidate_Selection"
                                + "/"
                                + str(merged_df["premz"][mer])
                                + "_ChemMNedges.tsv",
                                sep="\t",
                            )
                        one_candidate = one_candidate_selection_metfrag(
                            sm,
                            # sirius_df = sirius_df,
                            mbank_df=mbank_df,
                            # gnps_df = gnps_df ,
                            hmdb_df=hmdb_df,
                            Source="HM",
                        )
                        candidates_with_counts = add_count_column_metfrag(one_candidate)
                        candidates_with_counts.to_csv(
                            entry
                            + "/"
                            + "Candidate_Selection"
                            + "/"
                            + str(merged_df["premz"][mer])
                            + "sorted_candidate_list.tsv",
                            sep="\t",
                        )
                        if max(candidates_with_counts["Count"]) == 2:
                            sources_2_metfrag(candidates_with_counts, merged_df, mer, sirius_df=None)
                        if max(candidates_with_counts["Count"]) == 1:
                            sources_1_metfrag(candidates_with_counts, merged_df, mer, sirius_df=None)
                    # sS
                    elif (
                        len(sirius_df) > 0
                        and len(gnps_df) == 0
                        and len(mbank_df) == 0
                        and len(hmdb_df) == 0
                    ):

                        sirius_df["rank_ids"] = ["E_" + str(s) for s in range(len(sirius_df))]
                        sirius_df["Source"] = "MetFrag"

                        source_l1 = [*(list(sirius_df["Source"]))]

                        rank_l2 = [*(list(sirius_df["rank_ids"]))]

                        smiles_l3 = [*(list(sirius_df["SMILES"]))]

                        sm = pd.DataFrame(
                            list(zip(source_l1, rank_l2, smiles_l3)),
                            columns=["Source", "ranks", "SMILES"],
                        )

                        df_edge = chemMN_CandidateSelection(sm)

                        if df_edge is not None:

                            df_edge.to_csv(
                                entry
                                + "/"
                                + "Candidate_Selection"
                                + "/"
                                + str(merged_df["premz"][mer])
                                + "_ChemMNedges.tsv",
                                sep="\t",
                            )

                        one_candidate = one_candidate_selection_metfrag(
                            sm,
                            sirius_df=sirius_df,
                            # mbank_df = mbank_df,
                            # gnps_df = gnps_df ,
                            # hmdb_df = hmdb_df,
                            Source="S",
                        )
                        candidates_with_counts = add_count_column_metfrag(one_candidate)
                        candidates_with_counts.to_csv(
                            entry
                            + "/"
                            + "Candidate_Selection"
                            + "/"
                            + str(merged_df["premz"][mer])
                            + "sorted_candidate_list.tsv",
                            sep="\t",
                        )
                        print(candidates_with_counts["SMILES"])
                        if max(candidates_with_counts["Count"]) == 1:
                            sources_1_metfrag(candidates_with_counts, merged_df, mer, sirius_df)
                    # G
                    elif (
                        len(sirius_df) == 0
                        and len(gnps_df) > 0
                        and len(mbank_df) == 0
                        and len(hmdb_df) == 0
                    ):
                        gnps_df["rank_ids"] = ["G_" + str(s + 1) for s in range(len(gnps_df))]

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
                                entry
                                + "/"
                                + "Candidate_Selection"
                                + "/"
                                + str(merged_df["premz"][mer])
                                + "_ChemMNedges.tsv",
                                sep="\t",
                            )

                        one_candidate = one_candidate_selection_metfrag(
                            sm,
                            # sirius_df = sirius_df,
                            # mbank_df = mbank_df,
                            gnps_df=gnps_df,
                            # hmdb_df = hmdb_df,
                            Source="G",
                        )
                        candidates_with_counts = add_count_column_metfrag(one_candidate)
                        candidates_with_counts.to_csv(
                            entry
                            + "/"
                            + "Candidate_Selection"
                            + "/"
                            + str(merged_df["premz"][mer])
                            + "sorted_candidate_list.tsv",
                            sep="\t",
                        )
                        if max(candidates_with_counts["Count"]) == 1:
                            sources_1_metfrag(candidates_with_counts, merged_df, mer, sirius_df=None)
                    # M
                    elif (
                        len(sirius_df) == 0
                        and len(gnps_df) == 0
                        and len(mbank_df) > 0
                        and len(hmdb_df) == 0
                    ):
                        mbank_df["rank_ids"] = ["M_" + str(s + 1) for s in range(len(mbank_df))]

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
                                entry
                                + "/"
                                + "Candidate_Selection"
                                + "/"
                                + str(merged_df["premz"][mer])
                                + "_ChemMNedges.tsv",
                                sep="\t",
                            )

                        one_candidate = one_candidate_selection_metfrag(
                            sm,
                            # sirius_df = sirius_df,
                            mbank_df=mbank_df,
                            # gnps_df = gnps_df ,
                            # hmdb_df = hmdb_df,
                            Source="M",
                        )
                        candidates_with_counts = add_count_column_metfrag(one_candidate)
                        candidates_with_counts.to_csv(
                            entry
                            + "/"
                            + "Candidate_Selection"
                            + "/"
                            + str(merged_df["premz"][mer])
                            + "sorted_candidate_list.tsv",
                            sep="\t",
                        )
                        if max(candidates_with_counts["Count"]) == 1:
                            sources_1_metfrag(candidates_with_counts, merged_df, mer, sirius_df=None)
                    # H
                    elif (
                        len(sirius_df) == 0
                        and len(gnps_df) == 0
                        and len(mbank_df) == 0
                        and len(hmdb_df) > 0
                    ):
                        hmdb_df["rank_ids"] = ["H_" + str(s + 1) for s in range(len(hmdb_df))]

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
                                entry
                                + "/"
                                + "Candidate_Selection"
                                + "/"
                                + str(merged_df["premz"][mer])
                                + "_ChemMNedges.tsv",
                                sep="\t",
                            )

                        one_candidate = one_candidate_selection_metfrag(
                            sm,
                            #                                                                         sirius_df = sirius_df,
                            #                                                                         mbank_df = mbank_df,
                            #                                                                         gnps_df = gnps_df ,
                            hmdb_df=hmdb_df,
                            Source="H",
                        )
                        candidates_with_counts = add_count_column_metfrag(one_candidate)
                        candidates_with_counts.to_csv(
                            entry
                            + "/"
                            + "Candidate_Selection"
                            + "/"
                            + str(merged_df["premz"][mer])
                            + "sorted_candidate_list.tsv",
                            sep="\t",
                        )

                        if max(candidates_with_counts["Count"]) == 1:
                            sources_1_metfrag(candidates_with_counts, merged_df, mer, sirius_df=None)

                        if standards:
                            if not isNaN(merged_df["SMILES"][mer]):

                                merged_df.loc[mer, "MSILevel"] = 1

            merged_df.to_csv(
                entry
                + "/"
                + "mergedResults-with-one-Candidates.csv"
            )
            #print(merged_df)
            merged_df = checkSMILES_validity(resultcsv = entry
                + "/"
                + "mergedResults-with-one-Candidates.csv")
            merged_df.to_csv(
                entry
                + "/"
                + "mergedResults-with-one-Candidates.csv"
            )

# doesnt apply anymore here. because the argument is a single directory with results
#@p.provenance()
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

    results.to_csv(
        input_dir + "MetabolomicsResults/final_curationListVS" + listname + ".csv"
    )
    return results


#@p.provenance()
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

        frame.to_csv(result_csv)
        return frame

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

# naming is any string used to name or put an id to the sunburst as many sunbursts will be created
def sunburst(input_dir, input_csv, naming):
    
    cl = pd.read_csv(input_csv)
    class_data = cl[['superclass', 'class', 'subclass']]
    spclass = list(class_data['superclass']) # all superclasses
    uniq_spclass = list(np.unique(list(class_data['superclass']))) # only unique super classes
    uniq_spc = [s for s in uniq_spclass if 'nan' not in s ] # only unique super classes with no NA values
    print(len(uniq_spclass))
    clss = list(class_data['class'])
    uniq_class = list(np.unique(list(class_data['class'])))
    uniq_c = [s for s in uniq_class if 'nan' not in s ]
    len(uniq_class)
    sbclass = list(class_data['subclass'])
    uniq_sbclass = list(np.unique(list(class_data['subclass'])))
    uniq_sbc = [s for s in uniq_sbclass if 'nan' not in s ]
    len(uniq_sbclass)

    #all characters
    Names = ['Organic Compounds'] + uniq_spclass+uniq_class+uniq_sbclass

    df = pd.DataFrame(Names)
    df['values'] = ''
    df['parents'] = ''

    df = df.rename(columns={0: 'characters'})
    
    if "nan" in np.unique(df["characters"]):
    
        #for i, row in df.iterrows():
        for i, row in df[0:len(df)-2].iterrows():
            if 'Organic Compounds' in df['characters'][i]:
                df.loc[i, 'values'] = 0
                df.loc[i, 'parents'] = ''

            elif df['characters'][i] in uniq_spclass:

                df.loc[i, 'values'] = spclass.count(df['characters'][i])
                df.loc[i, 'parents'] = 'Organic Compounds'

            elif df['characters'][i] in uniq_class:

                df.loc[i, 'values'] = clss.count(df['characters'][i])
                df.loc[i, 'parents'] = 'Organic Compounds'

                df.loc[i, 'values'] = clss.count(df['characters'][i])
                clsp = class_data['superclass'][class_data[class_data['class'] == df['characters'][i]].index.tolist()[0]]
                df.loc[i, 'parents'] = clsp


            elif df['characters'][i] in uniq_sbclass:
                df.loc[i, 'values'] = sbclass.count(df['characters'][i])
                sbclsp = class_data['class'][class_data[class_data['subclass'] == df['characters'][i]].index.tolist()[0]]
                df.loc[i, 'parents'] = sbclsp
    else:
        for i, row in df.iterrows():
            if 'Organic Compounds' in df['characters'][i]:
                df.loc[i, 'values'] = 0
                df.loc[i, 'parents'] = ''

            elif df['characters'][i] in uniq_spclass:

                df.loc[i, 'values'] = spclass.count(df['characters'][i])
                df.loc[i, 'parents'] = 'Organic Compounds'

            elif df['characters'][i] in uniq_class:

                df.loc[i, 'values'] = clss.count(df['characters'][i])
                df.loc[i, 'parents'] = 'Organic Compounds'

                df.loc[i, 'values'] = clss.count(df['characters'][i])
                clsp = class_data['superclass'][class_data[class_data['class'] == df['characters'][i]].index.tolist()[0]]
                df.loc[i, 'parents'] = clsp


            elif df['characters'][i] in uniq_sbclass:
                df.loc[i, 'values'] = sbclass.count(df['characters'][i])
                sbclsp = class_data['class'][class_data[class_data['subclass'] == df['characters'][i]].index.tolist()[0]]
                df.loc[i, 'parents'] = sbclsp
    data = dict(character = df['characters'], parents = df['parents'], values = df['values'])
    fig = px.sunburst(
        data,
        names='character',
        parents='parents',
        values='values',
        
    )
    fig.update_layout(margin = dict(t=0, l=0, r=0, b=0))
    name_html = input_dir+"/"+naming+"_sunburst.html"
    print(name_html)
    fig.write_html(name_html)
    fig.show()
    return data

def list_all_folders(input_dir, list_of_words_to_remove, starting_words_for_all_folders = None):
    #define the input directory
    path = input_dir
    # list all files
    file = os.listdir(path)
    if starting_words_for_all_folders:
        # list all folders that start with "DS" and arent mzML files
        folders = [x for x in os.listdir(path) if x.startswith(starting_words_for_all_folders)]
    else:
        # list all folders that start with "DS" and arent mzML files
        folders = [x for x in os.listdir(path)]
    
    for i in list_of_words_to_remove:
        folders = [x for x in folders if not i in x]
        
    return folders

def chemMN(input_dir, input_csv, naming, name_col):
    df = pd.read_csv(input_csv)
    dbn= []
    for i, row in df.iterrows():
        for j, row in df.iterrows():
            if df['SMILES'][i] != df['SMILES'][j]:
                try:
                    ms = [Chem.MolFromSmiles(df['SMILES'][i]), Chem.MolFromSmiles(df['SMILES'][j])]
                    fps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=1024) for x in ms]
                    tn = DataStructs.FingerprintSimilarity(fps[0],fps[1])
                    dbn.append({
                        'Namei':df[name_col][i],
                        'Namej':df[name_col][j],
                        'i': df['SMILES'][i],
                        'j': df['SMILES'][j],
                        'Tanimoto': tn
                    })
                except:
                    pass
        #print(i)
    db_edge = pd.DataFrame(dbn)
    db_edge.to_csv(input_dir+ "/"+ naming+ "_allVSall.csv")

    dfe = []
    x=0
    for i, row in db_edge.iterrows():        
        if db_edge['Tanimoto'][i] >= 0.85:
            x=x+1
            dfe.append({
                'Start':db_edge['Namei'][i],
                'End':db_edge['Namej'][i],
                'Tanimoto':db_edge['Tanimoto'][i]
            })
    new_df = pd.DataFrame(dfe)
    new_df['Start'] = new_df['Start'].astype(str)
    new_df['End'] = new_df['End'].astype(str)
    new_df['StartAtt']=np.nan
    new_df['EndAtt']=np.nan
    for i, row in new_df.iterrows():
        for j, row in df.iterrows():
            if new_df['Start'][i]==df[name_col][j]:
                new_df.loc[i, 'StartAtt'] = df['superclass'][j]
    for i, row in new_df.iterrows():
        for j, row in df.iterrows():
            if new_df['End'][i]==df[name_col][j]:
                new_df.loc[i, 'EndAtt'] = df['superclass'][j]

    new_df['sorted_names'] = new_df.apply(lambda row: '-'.join(sorted([row['Start'], row['End']])), axis=1)
    new_df = new_df.drop_duplicates(subset=["sorted_names"], keep="last")
    new_df.to_csv(input_dir + "/" + naming + "_chemMN_Cytoscape.tsv", sep='\t')
    return new_df


#Define input directory, keep all files in same directory and scripts so getwd works
# entry = sys.argv[1] # mzml_result
# sirius = sys.argv[2] == 'True'
# sub_dir = sys.argv[2]
# score_thresh = float(sys.argv[3])
#sirius_dirs = sys.argv[5:]

# provenance_result = spec_postproc(entry, Source = "all")

# if sirius:
#     sirius_postproc(entry, sirius_dirs)   
#     CandidateSelection_SimilarityandIdentity(entry, standards = False)
# else:    
# metfrag_postproc(entry, sub_dir, score_thresh) # add args.candidates
#CandidateSelection_SimilarityandIdentity_Metfrag(entry, standards = False)



# from ruamel.yaml.main import YAML

# with open("provenance_python.yaml", "w") as filehandle:
#     yaml = YAML()
#     yaml.default_flow_style = False
#     yaml.indent = 4
#     yaml.block_seq_indent = 2
#     yaml.dump(
#         { "fn_module": provenance_result.artifact.fn_module,
#           "fn_name": provenance_result.artifact.fn_name,
#           "run_info": provenance_result.artifact.run_info,
#           "inputs": provenance_result.artifact.inputs
#           },
#         filehandle,
#     )

# import provenance.vis as vis
# plot = vis.visualize_lineage(provenance_result)


def metfrag_postproc(metfrag_candidates, score_thresh):
    if os.path.exists(metfrag_candidates):
        metfrag_res = pd.read_csv(metfrag_candidates)
            if len(metfrag_res)>0:
                metfrag_res = metfrag_res[metfrag_res['Score'] >= score_thresh]
                metfrag_res.to_csv(metfrag_candidates)
    return metfrag_candidates

metfrag_candidate_list = []

def metfrag_append(ms1data, metfrag_candidate_list):
    if os.path.exists(ms1data):
        msp = pd.read_csv(ms1data)
        msp["MetFragCSV"] = np.nan
        # for each mz
        for mz, row in msp.iterrows():
            # make a list of files with this mz
            files_for_mz = []
            for file in metfrag_candidate_list:
                if str(msp["premz"][mz]) in file:
                    files_for_mz.append(file)
                            metfrag_res.to_csv(msp["MetFragCSV"][mz])
        msp.to_csv(ms1data)  


# args.candidates will be a python list
parser = argparse.ArgumentParser()
parser.add_argument("--ms1data")
parser.add_argument("--candidates", action = "append")
parser.add_argument("--score_thresh")

# parser.add_argument("--gnps")
# parser.add_argument("--hmdb")
# parser.add_argument("--mbank")

args = parser.parse_args()


metfrag_postproc(ms1data, candidates, score_thresh)
