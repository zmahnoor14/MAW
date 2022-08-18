#!/usr/bin/env python
# coding: utf-8

# In[1]:


from platform import python_version
print(python_version())


# In[2]:


import pandas as pd
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
import wget
import string
import urllib.parse
import openpyxl
import statistics
import sys
from itertools import chain


# In[3]:


from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem import PandasTools


# In[5]:


input_dir = "/Users/mahnoorzulfiqar/OneDriveUNI/MAW-data/StandardSMarinoi_Data"
input_dir


# In[23]:


def chemMN_CandidateSelection(df, tn_sim = 0.85):
    
    """chemMN_CandidateSelection function is used to generate a Cytoscape readable tsv file.
    This file contains start(starting SMILES) and end(target SMILES) nodes and the tanimoto 
    similarity scores between the nodes. User can visualize the structural similarity 
    between the given SMILES. It provides an "ALL against ALL" network.
    
    Parameters: 
    df: dataframe that contains "SMILES", "ranks", "Source". This function is specifically for 
    candidate selection and so these columns are necessary.
    
    
    Returns:
    dataframe: it returns a df with follwoing columns to be loaded into Cytoscape.
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
    #one_df = []
    # define empty variable to save the edges
    dbn= []
    # for each entry in the df
    for i, row in df.iterrows():
        # to compare each element with each other element of the df
        for j, row in df.iterrows():
            try:
                # calcultae tanimoto
                ms = [Chem.MolFromSmiles(df['SMILES'][i]), Chem.MolFromSmiles(df['SMILES'][j])]
                fps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in ms]
                tn = DataStructs.FingerprintSimilarity(fps[0],fps[1])
                
                # save all entries to a matrix
                dbn.append({
                    'Name_i':df['ranks'][i],
                    'Name_j':df['ranks'][j],
                    'i': df['SMILES'][i],
                    'j': df['SMILES'][j],
                    'Source_i':df['Source'][i],
                    'Source_j':df['Source'][j],
                    'Tanimoto': tn
                })

            except Exception as e:
                print(e.string)
                pass

    # save chemical similarities                    
    db_edgenode = pd.DataFrame(dbn)

    # another empty variable to store the results for final tsv file
    dfe = []
    
    # heavy atoms for MCSS Calculation
    heavy_atoms = ['C', 'N', 'P', 'O', 'S']
    
    # for the previous dataframe
    for i, row in db_edgenode.iterrows():      
        # if the tanimoto > 0.85 for high similarity
        if db_edgenode['Tanimoto'][i] >= tn:
            
            # calculate MCSS
            n = [Chem.MolFromSmiles(db_edgenode['i'][i]),Chem.MolFromSmiles(db_edgenode['j'][i])]
            res = rdFMCS.FindMCS(n)
            sm_res = Chem.MolToSmiles(Chem.MolFromSmarts(res.smartsString))
            
            
            # Check if the MCSS has one of the heavy atoms and whether they are
            # more than 3
            elem = [ele for ele in heavy_atoms if(ele in sm_res)]
            if elem and len(sm_res)>=3:
                MCSS_SMILES = Chem.MolToSmiles(Chem.MolFromSmarts(res.smartsString))
            
            # save everything into a dataframe
            dfe.append({
                'Start':db_edgenode['Name_i'][i],
                'End':db_edgenode['Name_j'][i],
                'Tanimoto':db_edgenode['Tanimoto'][i],
                'Start_SMILES':db_edgenode['i'][i],
                'End_SMILES':db_edgenode['j'][i],
                'Start_Source':db_edgenode['Source_i'][i],
                'End_Source':db_edgenode['Source_j'][i],
                'MCSS': MCSS_SMILES
            })
    df_edge = pd.DataFrame(dfe)
    # generate a column called sorted_row which contains ids of the start and end nodes as a list
    df_edge['Start'] = df_edge['Start'].astype(str)
    df_edge['End'] = df_edge['End'].astype(str)
    df_edge['sorted_row'] = [sorted([a,b]) for a,b in zip(df_edge.Start,df_edge.End)]
    df_edge['sorted_row'] = df_edge['sorted_row'].astype(str)
    df_edge.drop_duplicates(subset=['sorted_row'], inplace=True)

    #nodes= []
    #for i, row in df.iterrows():
        #n = df['ranks'][i]
        #nodes.append({
            #'nodes':n
        #})

    #node= pd.DataFrame(nodes)

    return(df_edge)


# In[24]:


def one_candidate_selection(df, Source = "SGHM", tn_ident = 0.99):
    
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
        if Source == "SGHM" or Source == "SGM" or Source == "SGH" or Source == "SHM" or Source == "SG" or Source == "SM" or Source == "SH":
            # sirius_df comes from within the function CandidateSelection_SimilarityandIdentity
            for sirius_i, row in sirius_df.iterrows():
                # calculate tanimoto
                try:
                    ms = [Chem.MolFromSmiles(df["SMILES"][smiles]), Chem.MolFromSmiles(sirius_df["smiles"][sirius_i])]
                    fps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in ms]
                    tn = DataStructs.FingerprintSimilarity(fps[0],fps[1])
                    # since we are dealing with idenity here so tanimoto of 0.99 is appropriate
                    if tn >= tn_ident:
                        
                        # if SIRIUS is blank, add the SIRIUS id
                        if isNaN(df["SIRIUS"][smiles]):

                            df.loc[smiles, "SIRIUS"] = sirius_df["rank_ids"][sirius_i]
                        # if not empty, add SIRIUS id, with a comma
                        else:
                            df.loc[smiles, "SIRIUS"] = str(df["SIRIUS"][smiles]) + ", "+ sirius_df["rank_ids"][sirius_i]
                            
                       
                except Exception as e:
                    print(e.string)
                    pass
        
        
        # If the Source contains GNPS
        if Source == "SGHM" or Source == "SGM" or Source == "SGH" or Source == "GHM" or Source == "SG" or Source == "GM" or Source == "GH":

            # gnps_df comes from within the function CandidateSelection_SimilarityandIdentity
            for gnps_i, row in gnps_df.iterrows():
                try:
                    # calculate tanimoto
                    ms = [Chem.MolFromSmiles(df["SMILES"][smiles]), Chem.MolFromSmiles(gnps_df["GNPSSMILES"][gnps_i])]
                    fps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in ms]
                    tn = DataStructs.FingerprintSimilarity(fps[0],fps[1])
                    
                    # since we are dealing with idenity here so tanimoto of 0.99 is appropriate
                    if tn >= tn_ident:
                        
                        # if GNPS is blank, add the GNPS id
                        if isNaN(df["GNPS"][smiles]):

                            df.loc[smiles, "GNPS"] = gnps_df["rank_ids"][gnps_i]
                        # if not empty, add GNPS id, with a comma
                        else:
                            df.loc[smiles, "GNPS"] = str(df["GNPS"][smiles]) + ", "+ gnps_df["rank_ids"][gnps_i]          
                        
                except Exception as e:
                    print(e.string)
                    pass

        # If the source contains MassBank
        if Source == "SGHM" or Source == "SGM" or Source == "SHM" or Source == "GHM" or Source == "SM" or Source == "GM" or Source == "HM":
            # mbank_df comes from within the function CandidateSelection_SimilarityandIdentity
            for mbank_i, row in mbank_df.iterrows():
                try:
                    # calculate tanimoto
                    ms = [Chem.MolFromSmiles(df["SMILES"][smiles]), Chem.MolFromSmiles(mbank_df["MBSMILES"][mbank_i])]
                    fps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in ms]
                    tn = DataStructs.FingerprintSimilarity(fps[0],fps[1])
                    
                    # since we are dealing with idenity here so tanimoto of 0.99 is appropriate
                    if tn >= tn_ident:
                        
                        # if MassBank is blank, add the MassBank id
                        if isNaN(df["MassBank"][smiles]):

                            df.loc[smiles, "MassBank"] = mbank_df["rank_ids"][mbank_i]
                        # if not empty, add MassBank id, with a comma
                        else:
                            df.loc[smiles, "MassBank"] = str(df["MassBank"][smiles]) + ", "+ mbank_df["rank_ids"][mbank_i]
                            
                       
                except Exception as e:
                    print(e.string)
                    pass


        # If the source contains HMDB
        if Source == "SGHM" or Source == "SGH" or Source == "SHM" or Source == "GHM" or Source == "SH" or Source == "GH" or Source == "HM":


            for hmdb_i, row in hmdb_df.iterrows():
                try:
                    # calculate tanimoto
                    ms = [Chem.MolFromSmiles(df["SMILES"][smiles]), Chem.MolFromSmiles(hmdb_df["HMDBSMILES"][hmdb_i])]
                    fps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=2048) for x in ms]
                    tn = DataStructs.FingerprintSimilarity(fps[0],fps[1])
                    
                    # since we are dealing with idenity here so tanimoto of 0.99 is appropriate
                    if tn >= tn_ident:
                        
                        # if HMDB is blank, add the HMDB id
                        if isNaN(df["HMDB"][smiles]):

                            df.loc[smiles, "HMDB"] = hmdb_df["rank_ids"][hmdb_i]
                        # if not empty, add HMDB id, with a comma
                        else:
                            df.loc[smiles, "HMDB"] = str(df["HMDB"][smiles]) + ", "+ hmdb_df["rank_ids"][hmdb_i]
                            
                       
                except Exception as e:
                    print(e.string)
                    pass
 
    return(df)


# In[25]:


def CandidateSelection_SimilarityandIdentity(input_dir):

    # entry is all files and folders in input_dir
    for entry in os.listdir(input_dir):
         #if the entry is also a directory
        if os.path.isdir(os.path.join(input_dir, entry)):


            # reach spectra_dereplication folder
            sub_dir_spec = input_dir + "/" + entry + '/spectral_dereplication/'
            # reach SIRIUS results
            sub_dir_sir = input_dir + "/" + entry + '/insilico/SIRIUS/'

            #list of all csv files in the spectral dereplication foler
            spec_msp_csv = (glob.glob(input_dir + "/" + entry + '/spectral_dereplication' +'/*.csv'))
            # Sirius csv result file
            sir_msp_csv = input_dir + "/" + entry + "/insilico/MS1DATA.csv"

            # if both exist; which should be the case, even in case of 0 results
            if os.path.exists(sir_msp_csv) and os.path.exists(spec_msp_csv[0]):

                # read both csv files
                spec_msv = pd.read_csv(spec_msp_csv[0])
                sir_msv = pd.read_csv(sir_msp_csv)

                spec_msv = spec_msv[[
                    'id_X',
                    'premz',
                    'rtmin',
                    'rtmax',
                    'rtmed',
                    'rtmean',
                    'col_eng',
                    'pol',
                    'int',
                    'source_file',
                    'mbank_results_csv',
                    ]]
                sir_msv = sir_msv[['id_X',
                    'premz',
                    'rtmed',
                    'rtmean',
                    'int',
                    'col_eng',
                    'pol',
                    'ms2Peaks',
                    'ms1Peaks',
                    'sirius_result_dir']]

                merged_df = sir_msv.merge(spec_msv, 
                              how='inner', 
                              left_on=['premz', 'rtmed','rtmean','int','col_eng','pol'], 
                              right_on=['premz','rtmed','rtmean','int','col_eng','pol'])

                for mer, rows in merged_df.iterrows():
                    print(mer)

                    sirius_csv = merged_df["sirius_result_dir"][mer].replace("./", input_dir+"/")
                    mbank_csv = merged_df["mbank_results_csv"][mer].replace("./", input_dir+"/")
                    gnps_csv = merged_df["mbank_results_csv"][mer].replace("./", input_dir+"/").replace('mbank', 'gnps').replace('MassBank', 'GNPS')
                    hmdb_csv = merged_df["mbank_results_csv"][mer].replace("./", input_dir+"/").replace('mbank', 'hmdb').replace('MassBank', 'HMDB')


                    if os.path.exists(sirius_csv) and os.path.exists(gnps_csv) and os.path.exists(mbank_csv) and os.path.exists(hmdb_csv):


                        sirius_df = pd.read_csv(sirius_csv)
                        sirius_df = sirius_df.drop_duplicates('smiles')

                        mbank_df = pd.read_csv(mbank_csv)
                        mbank_df = mbank_df.drop_duplicates('MBSMILES')

                        gnps_df = pd.read_csv(gnps_csv)
                        gnps_df = gnps_df.drop_duplicates('GNPSSMILES')

                        hmdb_df = pd.read_csv(hmdb_csv)
                        hmdb_df = hmdb_df.drop_duplicates('HMDBSMILES')

                        #1 SGHM
                        if len(sirius_df) > 0 and len(gnps_df) > 0 and len(mbank_df) > 0 and len(hmdb_df) > 0:
                            mbank_df["rank_ids"] = ["M_" + str(s+1) for s in range(len(mbank_df))]

                            gnps_df["rank_ids"] = ["G_" + str(s+1) for s in range(len(gnps_df))]

                            hmdb_df["rank_ids"] = ["H_" + str(s+1) for s in range(len(hmdb_df))]

                            sirius_df["rank_ids"] = ["S_" + str(s) for s in sirius_df["rank"]]
                            sirius_df["Source"] = "SIRIUS"


                            source_l1 = [*(list(sirius_df["Source"])), 
                               *(list(gnps_df["Source"]))
                               ,*(list(mbank_df["Source"])),
                                        *(list(hmdb_df["Source"]))]

                            rank_l2 = [*(list(sirius_df["rank_ids"])), 
                               *(list(gnps_df["rank_ids"]))
                               ,*(list(mbank_df["rank_ids"])), 
                                      *(list(hmdb_df["rank_ids"]))]

                            smiles_l3 = [*(list(sirius_df["smiles"])), 
                               *(list(gnps_df["GNPSSMILES"]))
                               ,*(list(mbank_df["MBSMILES"])),
                                        *(list(hmdb_df["HMDBSMILES"]))]

                            sm = pd.DataFrame(list(zip(source_l1, rank_l2, smiles_l3)), columns = ["Source", "ranks", "SMILES"])   
                            
                            df_edge = chemMN_CandidateSelection(sm)

                            df_edge.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_ChemMNedges.tsv", sep='\t')
                            one_candidate = one_candidate_selection(sm, Source = "SGHM")
                            one_candidate.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_one_candidate_list.tsv", sep='\t')
                            
                            
                            df_sources = pd.DataFrame({"SIRIUS": one_candidate["SIRIUS"], 
                                                       "GNPS": one_candidate["GNPS"],
                                                       "HMDB": one_candidate["HMDB"],
                                                       "MassBank": one_candidate["MassBank"]})
                            
                            
                            # now check which rows have a value

                            index_SIRIUS = [x for x, row in df_sources.iterrows() if not isNaN(df_sources["SIRIUS"][x])]
                            index_GNPS = [x for x, row in df_sources.iterrows() if not isNaN(df_sources["GNPS"][x])]
                            index_MassBank = [x for x, row in df_sources.iterrows() if not isNaN(df_sources["MassBank"][x])]
                            index_HMDB = [x for x, row in df_sources.iterrows() if not isNaN(df_sources["HMDB"][x])]

                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            

                        #2 SGM
                        elif  len(sirius_df) > 0 and len(gnps_df) > 0 and len(mbank_df) > 0 and len(hmdb_df) == 0:

                            mbank_df["rank_ids"] = ["M_" + str(s+1) for s in range(len(mbank_df))]

                            gnps_df["rank_ids"] = ["G_" + str(s+1) for s in range(len(gnps_df))]

                            sirius_df["rank_ids"] = ["S_" + str(s) for s in sirius_df["rank"]]
                            sirius_df["Source"] = "SIRIUS"


                            source_l1 = [*(list(sirius_df["Source"])), 
                               *(list(gnps_df["Source"]))
                               ,*(list(mbank_df["Source"]))]

                            rank_l2 = [*(list(sirius_df["rank_ids"])), 
                               *(list(gnps_df["rank_ids"]))
                               ,*(list(mbank_df["rank_ids"]))]

                            smiles_l3 = [*(list(sirius_df["smiles"])), 
                               *(list(gnps_df["GNPSSMILES"]))
                               ,*(list(mbank_df["MBSMILES"]))]

                            sm = pd.DataFrame(list(zip(source_l1, rank_l2, smiles_l3)), columns = ["Source", "ranks", "SMILES"])


                            df_edge = chemMN_CandidateSelection(sm)

                            df_edge.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_ChemMNedges.tsv", sep='\t')
                            one_candidate = one_candidate_selection(sm, Source = "SGM")
                            one_candidate.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_one_candidate_list.tsv", sep='\t')

                        #3 SHM
                        elif  len(sirius_df) > 0 and len(gnps_df) == 0 and len(mbank_df) > 0 and len(hmdb_df) > 0:

                            mbank_df["rank_ids"] = ["M_" + str(s+1) for s in range(len(mbank_df))]

                            #gnps_df["rank_ids"] = ["G_" + str(s+1) for s in range(len(gnps_df))]

                            hmdb_df["rank_ids"] = ["H_" + str(s+1) for s in range(len(hmdb_df))]

                            sirius_df["rank_ids"] = ["S_" + str(s) for s in sirius_df["rank"]]
                            sirius_df["Source"] = "SIRIUS"


                            source_l1 = [*(list(sirius_df["Source"]))
                               ,*(list(mbank_df["Source"])),
                                        *(list(hmdb_df["Source"]))]

                            rank_l2 = [*(list(sirius_df["rank_ids"]))
                               ,*(list(mbank_df["rank_ids"])), 
                                      *(list(hmdb_df["rank_ids"]))]

                            smiles_l3 = [*(list(sirius_df["smiles"]))
                               ,*(list(mbank_df["MBSMILES"])),
                                        *(list(hmdb_df["HMDBSMILES"]))]

                            sm = pd.DataFrame(list(zip(source_l1, rank_l2, smiles_l3)), columns = ["Source", "ranks", "SMILES"])


                            df_edge = chemMN_CandidateSelection(sm)

                            df_edge.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_ChemMNedges.tsv", sep='\t')
                            one_candidate = one_candidate_selection(sm, Source = "SHM")
                            one_candidate.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_one_candidate_list.tsv", sep='\t')

                        #4 SGH
                        elif  len(sirius_df) > 0 and len(gnps_df) > 0 and len(mbank_df) == 0 and len(hmdb_df) > 0:
                            #mbank_df["rank_ids"] = ["M_" + str(s+1) for s in range(len(mbank_df))]

                            gnps_df["rank_ids"] = ["G_" + str(s+1) for s in range(len(gnps_df))]

                            hmdb_df["rank_ids"] = ["H_" + str(s+1) for s in range(len(hmdb_df))]

                            sirius_df["rank_ids"] = ["S_" + str(s) for s in sirius_df["rank"]]
                            sirius_df["Source"] = "SIRIUS"


                            source_l1 = [*(list(sirius_df["Source"]))
                               ,*(list(gnps_df["Source"])),
                                        *(list(hmdb_df["Source"]))]

                            rank_l2 = [*(list(sirius_df["rank_ids"]))
                               ,*(list(gnps_df["rank_ids"])), 
                                      *(list(hmdb_df["rank_ids"]))]

                            smiles_l3 = [*(list(sirius_df["smiles"]))
                               ,*(list(gnps_df["GNPSSMILES"])),
                                        *(list(hmdb_df["HMDBSMILES"]))]

                            sm = pd.DataFrame(list(zip(source_l1, rank_l2, smiles_l3)), columns = ["Source", "ranks", "SMILES"])

                            df_edge = chemMN_CandidateSelection(sm)

                            df_edge.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_ChemMNedges.tsv", sep='\t')
                            one_candidate = one_candidate_selection(sm, Source = "SGH")
                            one_candidate.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_one_candidate_list.tsv", sep='\t')
                            
                        #5 GHM
                        elif  len(sirius_df) == 0 and len(gnps_df) > 0 and len(mbank_df) > 0 and len(hmdb_df) > 0:
                            mbank_df["rank_ids"] = ["M_" + str(s+1) for s in range(len(mbank_df))]

                            gnps_df["rank_ids"] = ["G_" + str(s+1) for s in range(len(gnps_df))]

                            hmdb_df["rank_ids"] = ["H_" + str(s+1) for s in range(len(hmdb_df))]

                            #sirius_df["rank_ids"] = ["S_" + str(s) for s in sirius_df["rank"]]
                            #sirius_df["Source"] = "SIRIUS"


                            source_l1 = [*(list(gnps_df["Source"]))
                               ,*(list(mbank_df["Source"])),
                                        *(list(hmdb_df["Source"]))]

                            rank_l2 = [*(list(gnps_df["rank_ids"]))
                               ,*(list(mbank_df["rank_ids"])), 
                                      *(list(hmdb_df["rank_ids"]))]

                            smiles_l3 = [*(list(gnps_df["GNPSSMILES"]))
                               ,*(list(mbank_df["MBSMILES"])),
                                        *(list(hmdb_df["HMDBSMILES"]))]

                            sm = pd.DataFrame(list(zip(source_l1, rank_l2, smiles_l3)), columns = ["Source", "ranks", "SMILES"])

                            df_edge = chemMN_CandidateSelection(sm)

                            df_edge.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_ChemMNedges.tsv", sep='\t')
                            one_candidate = one_candidate_selection(sm, Source = "GHM")
                            one_candidate.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_one_candidate_list.tsv", sep='\t')
                            
                        #6 SG
                        elif  len(sirius_df) > 0 and len(gnps_df) > 0 and len(mbank_df) == 0 and len(hmdb_df) == 0:
                            #mbank_df["rank_ids"] = ["M_" + str(s+1) for s in range(len(mbank_df))]

                            gnps_df["rank_ids"] = ["G_" + str(s+1) for s in range(len(gnps_df))]

                            sirius_df["rank_ids"] = ["S_" + str(s) for s in sirius_df["rank"]]
                            sirius_df["Source"] = "SIRIUS"


                            source_l1 = [*(list(sirius_df["Source"])), 
                               *(list(gnps_df["Source"]))]
                               #,*(list(mbank_df["Source"]))]

                            rank_l2 = [*(list(sirius_df["rank_ids"])), 
                               *(list(gnps_df["rank_ids"]))]
                               #,*(list(mbank_df["rank_ids"]))]

                            smiles_l3 = [*(list(sirius_df["smiles"])), 
                               *(list(gnps_df["GNPSSMILES"]))]
                               #,*(list(mbank_df["MBSMILES"]))]

                            sm = pd.DataFrame(list(zip(source_l1, rank_l2, smiles_l3)), columns = ["Source", "ranks", "SMILES"])


                            df_edge = chemMN_CandidateSelection(sm)

                            df_edge.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_ChemMNedges.tsv", sep='\t')
                            
                            one_candidate = one_candidate_selection(sm, Source = "SG")
                            one_candidate.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_one_candidate_list.tsv", sep='\t')

                        #7 SH
                        elif  len(sirius_df) > 0 and len(gnps_df) == 0 and len(mbank_df) == 0 and len(hmdb_df) > 0:
                            #mbank_df["rank_ids"] = ["M_" + str(s+1) for s in range(len(mbank_df))]

                            #gnps_df["rank_ids"] = ["G_" + str(s+1) for s in range(len(gnps_df))]

                            hmdb_df["rank_ids"] = ["H_" + str(s+1) for s in range(len(hmdb_df))]

                            sirius_df["rank_ids"] = ["S_" + str(s) for s in sirius_df["rank"]]
                            sirius_df["Source"] = "SIRIUS"


                            source_l1 = [*(list(sirius_df["Source"])),
                                        *(list(hmdb_df["Source"]))]

                            rank_l2 = [*(list(sirius_df["rank_ids"])), 
                                      *(list(hmdb_df["rank_ids"]))]

                            smiles_l3 = [*(list(sirius_df["smiles"])),
                                        *(list(hmdb_df["HMDBSMILES"]))]

                            sm = pd.DataFrame(list(zip(source_l1, rank_l2, smiles_l3)), columns = ["Source", "ranks", "SMILES"])
                            df_edge = chemMN_CandidateSelection(sm)

                            df_edge.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_ChemMNedges.tsv", sep='\t')
                            
                            one_candidate = one_candidate_selection(sm, Source = "SH")
                            one_candidate.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_one_candidate_list.tsv", sep='\t')
                        #8 SM
                        elif  len(sirius_df) > 0 and len(gnps_df) == 0 and len(mbank_df) > 0 and len(hmdb_df) == 0:
                            mbank_df["rank_ids"] = ["M_" + str(s+1) for s in range(len(mbank_df))]

                            #gnps_df["rank_ids"] = ["G_" + str(s+1) for s in range(len(gnps_df))]

                            #hmdb_df["rank_ids"] = ["H_" + str(s+1) for s in range(len(hmdb_df))]

                            sirius_df["rank_ids"] = ["S_" + str(s) for s in sirius_df["rank"]]
                            sirius_df["Source"] = "SIRIUS"


                            source_l1 = [*(list(sirius_df["Source"]))
                               ,*(list(mbank_df["Source"]))]

                            rank_l2 = [*(list(sirius_df["rank_ids"]))
                               ,*(list(mbank_df["rank_ids"]))]

                            smiles_l3 = [*(list(sirius_df["smiles"]))
                               ,*(list(mbank_df["MBSMILES"]))]

                            sm = pd.DataFrame(list(zip(source_l1, rank_l2, smiles_l3)), columns = ["Source", "ranks", "SMILES"])


                            df_edge = chemMN_CandidateSelection(sm)

                            df_edge.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_ChemMNedges.tsv", sep='\t')
                            
                            one_candidate = one_candidate_selection(sm, Source = "SM")
                            one_candidate.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_one_candidate_list.tsv", sep='\t')
                        #9 GM
                        elif  len(sirius_df) == 0 and len(gnps_df) > 0 and len(mbank_df) > 0 and len(hmdb_df) == 0:
                            mbank_df["rank_ids"] = ["M_" + str(s+1) for s in range(len(mbank_df))]

                            gnps_df["rank_ids"] = ["G_" + str(s+1) for s in range(len(gnps_df))]

                            #hmdb_df["rank_ids"] = ["H_" + str(s+1) for s in range(len(hmdb_df))]

                            #sirius_df["rank_ids"] = ["S_" + str(s) for s in sirius_df["rank"]]
                            #sirius_df["Source"] = "SIRIUS"


                            source_l1 = [*(list(mbank_df["Source"])),
                                        *(list(gnps_df["Source"]))]

                            rank_l2 = [*(list(mbank_df["rank_ids"])), 
                                      *(list(gnps_df["rank_ids"]))]

                            smiles_l3 = [*(list(mbank_df["MBSMILES"])),
                                        *(list(gnps_df["GNPSSMILES"]))]

                            sm = pd.DataFrame(list(zip(source_l1, rank_l2, smiles_l3)), columns = ["Source", "ranks", "SMILES"])


                            df_edge = chemMN_CandidateSelection(sm)

                            df_edge.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_ChemMNedges.tsv", sep='\t')
                            
                            one_candidate = one_candidate_selection(sm, Source = "GM")
                            one_candidate.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_one_candidate_list.tsv", sep='\t')
                        #10 GH
                        elif  len(sirius_df) == 0 and len(gnps_df) > 0 and len(mbank_df) == 0 and len(hmdb_df) > 0:
                            #mbank_df["rank_ids"] = ["M_" + str(s+1) for s in range(len(mbank_df))]

                            gnps_df["rank_ids"] = ["G_" + str(s+1) for s in range(len(gnps_df))]

                            hmdb_df["rank_ids"] = ["H_" + str(s+1) for s in range(len(hmdb_df))]

                            #sirius_df["rank_ids"] = ["S_" + str(s) for s in sirius_df["rank"]]
                            #sirius_df["Source"] = "SIRIUS"


                            source_l1 = [*(list(gnps_df["Source"])),
                                        *(list(hmdb_df["Source"]))]

                            rank_l2 = [*(list(gnps_df["rank_ids"])),
                                      *(list(hmdb_df["rank_ids"]))]

                            smiles_l3 = [*(list(gnps_df["GNPSSMILES"])),
                                        *(list(hmdb_df["HMDBSMILES"]))]

                            sm = pd.DataFrame(list(zip(source_l1, rank_l2, smiles_l3)), columns = ["Source", "ranks", "SMILES"])
                            df_edge = chemMN_CandidateSelection(sm)



                            df_edge.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_ChemMNedges.tsv", sep='\t')
                            one_candidate = one_candidate_selection(sm, Source = "GH")
                            one_candidate.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_one_candidate_list.tsv", sep='\t')
                        #11 HM
                        elif  len(sirius_df) == 0 and len(gnps_df) == 0 and len(mbank_df) > 0 and len(hmdb_df) > 0:
                            mbank_df["rank_ids"] = ["M_" + str(s+1) for s in range(len(mbank_df))]

                            #gnps_df["rank_ids"] = ["G_" + str(s+1) for s in range(len(gnps_df))]

                            hmdb_df["rank_ids"] = ["H_" + str(s+1) for s in range(len(hmdb_df))]

                            #sirius_df["rank_ids"] = ["S_" + str(s) for s in sirius_df["rank"]]
                            #sirius_df["Source"] = "SIRIUS"


                            source_l1 = [*(list(mbank_df["Source"])),
                                        *(list(hmdb_df["Source"]))]

                            rank_l2 = [*(list(mbank_df["rank_ids"])), 
                                      *(list(hmdb_df["rank_ids"]))]

                            smiles_l3 = [*(list(mbank_df["MBSMILES"])),
                                        *(list(hmdb_df["HMDBSMILES"]))]

                            sm = pd.DataFrame(list(zip(source_l1, rank_l2, smiles_l3)), columns = ["Source", "ranks", "SMILES"])


                            df_edge = chemMN_CandidateSelection(sm)

                            df_edge.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_ChemMNedges.tsv", sep='\t')
                            one_candidate = one_candidate_selection(sm, Source = "HM")
                            one_candidate.to_csv(input_dir + "/" + entry + "/" + str(merged_df["premz"][mer]) + "_one_candidate_list.tsv", sep='\t')

                        # S
                        #elif  len(sirius_df) > 0 and len(gnps_df) == 0 and len(mbank_df) == 0 and len(hmdb_df) == 0:

                        # G
                        #elif  len(sirius_df) == 0 and len(gnps_df) > 0 and len(mbank_df) == 0 and len(hmdb_df) == 0:

                        # M
                        #elif  len(sirius_df) == 0 and len(gnps_df) == 0 and len(mbank_df) > 0 and len(hmdb_df) == 0:

                        # H
                        #elif  len(sirius_df) == 0 and len(gnps_df) == 0 and len(mbank_df) == 0 and len(hmdb_df) > 0:
                merged_df.to_csv("mergedResults_withCandidates.csv")


# In[26]:


#read csv
one_can = pd.read_csv("/Users/mahnoorzulfiqar/OneDriveUNI/MAW-data/StandardSMarinoi_Data/VN_211016_acetyl_carnitine/204.122756958008_one_candidate_list.tsv", sep = "\t")


# In[27]:


# delete count column which is wrong
del one_can["Count"]


# In[28]:


def add_count_column(df_one_candidate):
    # create new df only with the Sources column
    df = pd.DataFrame({"SIRIUS": one_can["SIRIUS"], 
                   "GNPS": one_can["GNPS"],
                   "MassBank": one_can["MassBank"], 
                   "HMDB": one_can["HMDB"]})
    
    # now check which rows have a value

    index_SIRIUS = [x for x, row in df.iterrows() if not isNaN(df["SIRIUS"][x])]
    index_GNPS = [x for x, row in df.iterrows() if not isNaN(df["GNPS"][x])]
    index_MassBank = [x for x, row in df.iterrows() if not isNaN(df["MassBank"][x])]
    index_HMDB = [x for x, row in df.iterrows() if not isNaN(df["HMDB"][x])]
    
    # make a list of the rows
    list_of_indices = index_SIRIUS + index_GNPS + index_MassBank + index_HMDB

    # count how mnay times one of the rows is appearing and add count
    count_list = [[x, list_of_indices.count(x)] for x in set(list_of_indices)]
    # add this info to one_can
    one_can["Count"] = [count_list[x][1] for x in range(len(count_list))]
    # sort the list by count in descending order
    sorted_count_one_candidate = df_one_candidate.sort_values(by = "Count", ascending = False)
    return(sorted_count_one_candidate)


# In[29]:


count_df = add_count_column(df_one_candidate = one_can)
count_df


# In[30]:


def sources_4(df_count):
    frame = [df_count[df_count["Source"] == "GNPS"], df_count[df_count["Source"] == "MassBank"], 
                 df_count[df_count["Source"] == "SIRIUS"], df_count[df_count["Source"] == "HMDB"]]
    f = pd.concat(frame)
    # extract rank number
    f["rank_num"] = [int(item) for sublist in list(f["ranks"]) for item in sublist if item.isdigit()]
    f["MSI-Level"] = "Level-2: Probable Structure"
    # sort according to rank number
    f_ranked = f.sort_values(by = "rank_num")
    return(f_ranked)


# In[31]:


def sources_3(df_count, sources):
    if "HMDB" not in sources:
        frame = [df_count[df_count["Source"] == "GNPS"], df_count[df_count["Source"] == "MassBank"], 
                 df_count[df_count["Source"] == "SIRIUS"]]
        f = pd.concat(frame)
        # extract rank number
        f["rank_num"] = [int(item) for sublist in list(f["ranks"]) for item in sublist if item.isdigit()]
        f["MSI-Level"] = "Level-2: Probable Structure"
        # sort according to rank number
        f_ranked = f.sort_values(by = "rank_num")

    elif "SIRIUS" not in sources:
        frame = [df_count[df_count["Source"] == "GNPS"], df_count[df_count["Source"] == "MassBank"], 
                 df_count[df_count["Source"] == "HMDB"]]
        f = pd.concat(frame)
        # extract rank number
        f["rank_num"] = [int(item) for sublist in list(f["ranks"]) for item in sublist if item.isdigit()]
        f["MSI-Level"] = "Level-2: Probable Structure"
        # sort according to rank number
        f_ranked = f.sort_values(by = "rank_num")

    elif "MassBank" not in sources:
        frame = [df_count[df_count["Source"] == "GNPS"], df_count[df_count["Source"] == "SIRIUS"], 
                 df_count[df_count["Source"] == "HMDB"]]
        f = pd.concat(frame)
        # extract rank number
        f["rank_num"] = [int(item) for sublist in list(f["ranks"]) for item in sublist if item.isdigit()]
        f["MSI-Level"] = "Level-2: Probable Structure"
        # sort according to rank number
        f_ranked = f.sort_values(by = "rank_num")


    elif "GNPS" not in sources:
        frame = [df_count[df_count["Source"] == "MassBank"], df_count[df_count["Source"] == "SIRIUS"], 
                 df_count[df_count["Source"] == "HMDB"]]
        f = pd.concat(frame)
        # extract rank number
        f["rank_num"] = [int(item) for sublist in list(f["ranks"]) for item in sublist if item.isdigit()]
        f["MSI-Level"] = "Level-2: Probable Structure"
        # sort according to rank number
        f_ranked = f.sort_values(by = "rank_num")
        
    return(f_ranked)


# In[32]:


def sources_2(df_count, sources):
    if "SIRIUS" not in sources and "HMDB" not in sources:
        frame = [df_count[df_count["Source"] == "GNPS"], df_count[df_count["Source"] == "MassBank"]]
        f = pd.concat(frame)
        # extract rank number
        f["rank_num"] = [int(item) for sublist in list(f["ranks"]) for item in sublist if item.isdigit()]
        f["MSI-Level"] = "Level-2: Probable Structure"
        # sort according to rank number
        f_ranked = f.sort_values(by = "rank_num")

    elif "MassBank" not in sources_2 and "HMDB" not in sources_2:
        frame2 = [df_count_2[df_count_2["Source"] == "GNPS"], df_count_2[df_count_2["Source"] == "SIRIUS"]]
        f2 = pd.concat(frame2)
        # extract rank number
        f2["rank_num"] = [int(item) for sublist in list(f2["ranks"]) for item in sublist if item.isdigit()]
        f2["MSI-Level"] = "Level-2: Probable Structure"
        # sort according to rank number
        f2_ranked = f2.sort_values(by = "rank_num")

    elif "GNPS" not in sources_2 and "HMDB" not in sources_2:
        frame2 = [df_count_2[df_count_2["Source"] == "MassBank"], df_count_2[df_count_2["Source"] == "SIRIUS"]]
        f2 = pd.concat(frame2)
        # extract rank number
        f2["rank_num"] = [int(item) for sublist in list(f2["ranks"]) for item in sublist if item.isdigit()]
        f2["MSI-Level"] = "Level-2: Probable Structure"
        # sort according to rank number
        f2_ranked = f2.sort_values(by = "rank_num")

    elif "MassBank" not in sources_2 and "SIRIUS" not in sources_2:
        frame2 = [df_count_2[df_count_2["Source"] == "GNPS"], df_count_2[df_count_2["Source"] == "HMDB"]]
        f2 = pd.concat(frame2)
        # extract rank number
        f2["rank_num"] = [int(item) for sublist in list(f2["ranks"]) for item in sublist if item.isdigit()]
        f2["MSI-Level"] = "Level-2: Probable Structure"
        # sort according to rank number
        f2_ranked = f2.sort_values(by = "rank_num")

    elif "MassBank" not in sources_2 and "GNPS" not in sources_2:
        frame2 = [df_count_2[df_count_2["Source"] == "SIRIUS"], df_count_2[df_count_2["Source"] == "HMDB"]]
        f2 = pd.concat(frame2)
        # extract rank number
        f2["rank_num"] = [int(item) for sublist in list(f2["ranks"]) for item in sublist if item.isdigit()]
        f2["MSI-Level"] = "Level-2: Probable Structure"
        # sort according to rank number
        f2_ranked = f2.sort_values(by = "rank_num")

    elif "GNPS" not in sources_2 and "SIRIUS" not in sources_2:
        frame2 = [df_count_2[df_count_2["Source"] == "MassBank"], df_count_2[df_count_2["Source"] == "HMDB"]]
        f2 = pd.concat(frame2)
        # extract rank number
        f2["rank_num"] = [int(item) for sublist in list(f2["ranks"]) for item in sublist if item.isdigit()]
        f2["MSI-Level"] = "Level-2: Probable Structure"
        # sort according to rank number
        f2_ranked = f2.sort_values(by = "rank_num")
        
    return(f_ranked)


# In[33]:


def sources_1(df_count, sources):
    if "GNPS" in sources:
        frame = [df_count[df_count["Source"] == pr_dbs[0]]]
        
        f = pd.concat(frame)
        f["rank_num"] = [int(item) for sublist in list(f["ranks"]) for item in sublist if item.isdigit()]
        f["MSI-Level"] = "Level-2: Probable Structure"
        f_ranked = f.sort_values(by = "rank_num") 
        
        
    elif "MassBank" in sources:
        frame = [df_count[df_count["Source"] == pr_dbs[1]]]
        
        f = pd.concat(frame)
        f["rank_num"] = [int(item) for sublist in list(f["ranks"]) for item in sublist if item.isdigit()]
        f["MSI-Level"] = "Level-3: Tentative Candidate"
        f_ranked = f.sort_values(by = "rank_num") 
        
    elif "SIRIUS" in sources:
        frame = [df_count[df_count["Source"] == pr_dbs[2]]]
        
        f = pd.concat(frame)
        f["rank_num"] = [int(item) for sublist in list(f["ranks"]) for item in sublist if item.isdigit()]
        f["MSI-Level"] = "Level-3: Tentative Candidate"
        f_ranked = f.sort_values(by = "rank_num") 
        
    elif "HMDB" in sources:
        frame = [df_count[df_count["Source"] == pr_dbs[3]]]
        
        f = pd.concat(frame)
        f["rank_num"] = [int(item) for sublist in list(f["ranks"]) for item in sublist if item.isdigit()]
        f["MSI-Level"] = "Level-3: Tentative Candidate"
        f_ranked = f.sort_values(by = "rank_num") 
        
    return(f_ranked)


# In[34]:


list_of_no = list(np.unique(count_df["Count"]))


# In[35]:


# priority of DBs
pr_dbs = ["GNPS", "MassBank", "SIRIUS", "HMDB"]

# if there are only two counts/sources each time for a candidate
if len(list_of_no) == 2:
    
    #take the first part of the count_df with 2 as count
    df_count_2 = [count_df[count_df["Count"] == x ] for x in list_of_no[::-1]][0]
    
    sources_2 = (list(np.unique(df_count_2["Source"])))
    
    if len(sources_2) == 4:
        
        
    if len(sources_2) == 3:
        
            
    if len(sources_2) == 2:

            
            
    df_count_1 = [count_df[count_df["Count"] == x ] for x in list_of_no[::-1]][1]
    
    sources_1 = (list(np.unique(df_count_1["Source"])))
    
    
    
    
    


# In[36]:


pd.concat([f2_ranked, f1_ranked])


# In[37]:


pr_dbs = ["GNPS", "MassBank", "SIRIUS", "HMDB"]
if len(list_of_no) == 4:
    df_count_4 = [count_df[count_df["Count"] == x ] for x in list_of_no[::-1]][0]
    frame = [high_2[high_2["Source"] == pr_dbs[0]], high_2[high_2["Source"] == pr_dbs[1]], 
             high_2[high_2["Source"] == pr_dbs[2]], high_2[high_2["Source"] == pr_dbs[3]]]

    f = pd.concat(frame)
    f["rank_num"] = [int(item) for sublist in list(f["ranks"]) for item in sublist if item.isdigit()]
    f["MSI-Level"] = "Level-2: Probable Structure"
    f.sort_values(by = "rank_num")   
    


# In[38]:


pr_dbs = ["GNPS", "MassBank", "SIRIUS", "HMDB"]
if len(list_of_no) == 3:
    df_count_3 = [count_df[count_df["Count"] == x ] for x in list_of_no[::-1]][0]
    
    if "HMDB" not in (list(np.unique(df_count_3["Source"]))):
        
        frame2 = [high_2[high_2["Source"] == pr_dbs[0]], 
                 high_2[high_2["Source"] == pr_dbs[1]], high_2[high_2["Source"] == pr_dbs[2]]]
        f2 = pd.concat(frame)
        f2["rank_num"] = [int(item) for sublist in list(f2["ranks"]) for item in sublist if item.isdigit()]
        f2["MSI-Level"] = "Level-2: Probable Structure"
        f2.sort_values(by = "rank_num")  
        
    if "SIRIUS" not in pr_dbs:
        
        frame = [high_2[high_2["Source"] == pr_dbs[0]], 
                 high_2[high_2["Source"] == pr_dbs[1]], high_2[high_2["Source"] == pr_dbs[3]]]

        f = pd.concat(frame)
        f["rank_num"] = [int(item) for sublist in list(f["ranks"]) for item in sublist if item.isdigit()]
        f["MSI-Level"] = "Level-2: Probable Structure"
        f.sort_values(by = "rank_num") 
        
        
    if "MassBank" not in pr_dbs:
        
        frame = [high_2[high_2["Source"] == pr_dbs[0]], 
                 high_2[high_2["Source"] == pr_dbs[2]], high_2[high_2["Source"] == pr_dbs[3]]]

        f = pd.concat(frame)
        f["rank_num"] = [int(item) for sublist in list(f["ranks"]) for item in sublist if item.isdigit()]
        f["MSI-Level"] = "Level-2: Probable Structure"
        f.sort_values(by = "rank_num")  
    
    if "GNPS" not in pr_dbs:
        
        frame = [high_2[high_2["Source"] == pr_dbs[1]], 
                 high_2[high_2["Source"] == pr_dbs[2]], high_2[high_2["Source"] == pr_dbs[3]]]

        f = pd.concat(frame)
        f["rank_num"] = [int(item) for sublist in list(f["ranks"]) for item in sublist if item.isdigit()]
        f["MSI-Level"] = "Level-2: Probable Structure"
        f.sort_values(by = "rank_num")   
        


# In[39]:


pr_dbs = ["GNPS", "MassBank", "SIRIUS", "HMDB"]
if len(list_of_no) == 2:
    df_count_2 = [count_df[count_df["Count"] == x ] for x in list_of_no[::-1]][0]
    
    if "GNPS" in :
        


# In[40]:



    
    
        df_count_3 = [count_df[count_df["Count"] == x ] for x in list_of_no[::-1]][1]
    df_count_2 = [count_df[count_df["Count"] == x ] for x in list_of_no[::-1]][2]
    df_count_1 = [count_df[count_df["Count"] == x ] for x in list_of_no[::-1]][3]
    
    
elif len(list_of_no) == 3:
    df_count_3 = [count_df[count_df["Count"] == x ] for x in list_of_no[::-1]][0]
    df_count_2 = [count_df[count_df["Count"] == x ] for x in list_of_no[::-1]][1]
    df_count_1 = [count_df[count_df["Count"] == x ] for x in list_of_no[::-1]][2]
    
elif len(list_of_no) == 2:
    df_count_2 = [count_df[count_df["Count"] == x ] for x in list_of_no[::-1]][0]
    df_count_1 = [count_df[count_df["Count"] == x ] for x in list_of_no[::-1]][1]
    
elif len(list_of_no) == 4:
    df_count_1 = [count_df[count_df["Count"] == x ] for x in list_of_no[::-1]][0]


# In[41]:


len((list(np.unique(df_count_2["Source"]))))


# In[42]:


frame = [high_2[high_2["Source"] == a], high_2[high_2["Source"] == b], high_2[high_2["Source"] == c]]


# In[43]:


f = pd.concat(frame)


# In[44]:


f["rank_num"] = [int(item) for sublist in list(f["ranks"]) for item in sublist if item.isdigit()]


# In[45]:


f.sort_values(by = "rank_num")


# In[ ]:




