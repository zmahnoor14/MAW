from unicodedata import name
import numpy as np
import pandas as pd
import pubchempy as pcp
import os
import re
import time
import wget
import urllib.parse
import argparse
import ast
import sys

from pybatchclassyfire import *
from pandas import json_normalize
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem import PandasTools

import plotly.express as px

def isNaN(string):
    return string != string

def chemMN(dataframe, naming, name_col):
    #df = pd.read_csv(dataframe)
    df = dataframe
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
    #db_edge.to_csv("/"+ naming+ "_allVSall.csv")

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
    new_df_filename = naming + "_chemMN_Cytoscape.tsv"
    new_df.to_csv(new_df_filename, sep='\t')
    return new_df_filename

def sunburst(dataframe, naming):
    
    #cl = pd.read_csv(input_csv)
    cl = dataframe
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
    name_html = "/"+naming+"_sunburst.html"
    print(name_html)
    fig.write_html(name_html)
    fig.show()
    print("DONE")
    return data

def merge_results(list_of_results):

    all_msi = []
    for result in list_of_results:
        if os.path.exists(result):
            all_msi.append(result)
            
    all_msi_df = pd.concat(map(pd.read_csv, all_msi), ignore_index=True)
    filename = "Final_combined_results.csv"
    all_msi_df.to_csv(filename)
    return all_msi_df, filename

parser = argparse.ArgumentParser(description='MAW-Summary')
parser.add_argument('--naming', type=str, help='can be based on your project name or condition')
# parser.add_argument('--name_col', type=str, default = "IUPAC", help='name of name column in the csv, generally MAW fetches the IUPAC names and not generic names')
# parser.add_argument('--list_of_results', type=str, action = "append", nargs='+', help ="list of paths of merged_results_with_one_Candidates.csv")


# Parse the command-line arguments
args = parser.parse_args()
# naming = args.naming
list_of_results = args.list_of_results

# naming = "coculture"
# list_of_results = ["/Users/mahnoorzulfiqar/OneDriveUNI/GitHub-Repos/MAW/cwl/Linking/exo_neg_mergedResults-with-one-Candidates.csv", 
# "/Users/mahnoorzulfiqar/OneDriveUNI/GitHub-Repos/MAW/cwl/Linking/exo_pos_mergedResults-with-one-Candidates.csv", 
# "/Users/mahnoorzulfiqar/OneDriveUNI/GitHub-Repos/MAW/cwl/Linking/endo_pos_mergedResults-with-one-Candidates.csv", 
# "/Users/mahnoorzulfiqar/OneDriveUNI/GitHub-Repos/MAW/cwl/Linking/endo_neg_mergedResults-with-one-Candidates.csv"]
# print("python script starts")
# print(naming)
print(list_of_results)

# new_results_list = []
# for i in list_of_results:
#     print(i)
#     new_results_list.append(str(i[0]))
df,combineresult = merge_results(list_of_results)


print(f"Combine results written to {combineresult!r}.", file= sys.stderr)

#chemMN(df, naming, name_col = "IUPAC")


# sunburst(df, naming)
