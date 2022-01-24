#!/usr/bin/env python
# coding: utf-8

# # Suspect List Curation

# In[5]:


####import libraries
import pubchempy as pcp
import pandas as pd
import numpy as np
import time
from rdkit import Chem
import re
import wget
import urllib.parse
import pubchempy as pcp
import pandas as pd
#from pandas import ExcelWriter
#from pandas import ExcelFile


# ## CHECK NAMES USING SMILES

# In[8]:


####read suspect list as csv
suspect_list = pd.read_csv("fully_curated_suspect_list_skeletonema.csv")


# In[9]:


#### Take only the entries with SMILES
bool_series = pd.notnull(suspect_list["SMILES"])
S_suspect_list = suspect_list[bool_series]


# In[7]:


####suspect list with all SMILES; Obtain compound with the SMILES 
####and check the name of the compound whether its same or different


db = [] #dataframe with change names
dc= [] #dataframe with same name
xs = [] #dataframe of smiles that have stereochemistry information and hinders pubchempy url

for index, row in S_suspect_list.iterrows():#get inside the suspect list with all SMILES
    try:#try this code
        sm = suspect_list["SMILES"][index]#read smiles from suspect list
        nm = suspect_list["Name"][index]#read names from suspect list
        x = pcp.get_compounds(sm, 'smiles')#pubchempy to obtain a compound using SMILES
        
        for compound in x:
            pac = compound.iupac_name#iupac from pubchempy
            
        for compound in x:
            syn = compound.synonyms#synonyms from pubchempy
            
        synLow = {synL.lower() for synL in syn}#converting syn to low case
        
        if nm.lower() not in synLow :#if they have different names
            
            nmC = syn[0]
            
            df1L = pcp.get_properties(['molecular_weight','molecular_formula', 'inchi'], sm, 'smiles', as_dataframe=True)
            
            for indx, row in df1L.iterrows():
                db.append(
                        {
                            'index':index,
                            'CID':indx,
                            'name':nmC,
                            'isomeric_smiles':sm,
                            'iupac':pac,
                            'synonyms':syn,
                            'molecular_weight':df1L.iloc[0,1],
                            'molecular_formula':df1L.iloc[0,0],
                            'inchi':df1L.iloc[0,2]
                        })
        else:#if the names are same
            df2L = pcp.get_properties(['molecular_weight','molecular_formula', 'inchi'], sm, 'smiles', as_dataframe=True)
            for ind, row in df2L.iterrows():
                dc.append(
                        {
                            'index':index,
                            'CID':ind,
                            'name':nm,
                            'isomeric_smiles':sm,
                            'iupac':pac,
                            'synonyms':syn,
                            'molecular_weight':df2L.iloc[0,1],
                            'molecular_formula':df2L.iloc[0,0],
                            'inchi':df2L.iloc[0,2]
                        })
    except Exception as e:#if there is error, ignore and pass; the error is due to bad request which means either the SMILES mentioned are not correct, not present or the functions couldnt call the url because of SMILES and URL syntax conflict
        xs.append(
            {
                'index':index,
                'isomeric_smiles':sm,
                'error_msg':e
            })
dx = pd.DataFrame(db)#dataframe with different or changed names
dy = pd.DataFrame(dc)#dataframe with same name
ds = pd.DataFrame(xs)#dataframe of smiles that have stereochemistry information and hinders pubchempy url


# In[11]:


#ds.to_csv("wrongSMILES.csv")


# In[13]:


#convert isomeric smiles into non isomeric smiles
ffe= [] #df for non-isomeric SMILES and their index in suspect list
ff4 = [] # df of SMILES with wrong syntax
for i, row in ds.iterrows():
    try:
        cv = ds['isomeric_smiles'][i]
        index = ds['index'][i]
        cvv = Chem.MolFromSmiles(cv)
        cxv = Chem.MolToSmiles(cvv, isomericSmiles = False)
        ff4.append({
            'indexX': index,
            'non_isomeric_smiles':cxv
        })
    except Exception as e:
        ffe.append({
            'indexX': index,
            'smiles':cv,
            'error':e
        })
ff43 = pd.DataFrame(ff4)##dataframe with information on non-isomeric SMILES and their index in suspect list
ffe3 = pd.DataFrame(ffe)##dataframe of SMILES with wrong Syntax (rdkit version)


# In[16]:


#ffe3.to_csv('wrongSMILESrdkit.csv')


# In[17]:


neww = []## for dataframe that contains information on non-isomeric smiles
old = []## for wrong SMILES syntax for PubChem

# take the df with canonical smiles 
for i, row in ff43.iterrows():
    try:
        
        c = ff43['non_isomeric_smiles'][i]#take canonical_smiles
        
        x = pcp.get_compounds(c, 'smiles')# obtain compound using canonical smiles
        
        for compound in x:
            pac = compound.iupac_name#iupac from pubchempy
            
        for compound in x:
            syn = compound.synonyms#synonyms from pubchempy
            
        indexX = ff43['indexX'][i]
        
        # synonyms present, create a df with following list
        if syn:
            neww.append(
                 {
                'indexX':indexX,
                'ID':x,
                'name':syn[0],
                'non_isomeric_smiles':c,
                'iupac':pac,
                'synonyms':syn,
                })
        
        else:
            neww.append(
                {
                'indexX':indexX,
                'ID':x,
                'name':pac,
                'non_isomeric_smiles':c,
                'iupac':pac,
                'synonyms':syn,
                    })
            
    # this exception is made when there are non standardized canonical smiles        
    except Exception as e:
        old.append({
            'indexX':indexX,
            'non_isomeric_smiles':c,
            'error': e
        })
N = pd.DataFrame(neww)#dataframe that contains information on non-isomeric smiles
O = pd.DataFrame(old)#df with wild card or some form of non-standardzied smiles


# In[18]:


O


# In[29]:


#O.to_csv('wrongSMILESPubChem.csv')


# In[19]:


N
#note that some compounds dont have synonyms, so I used iupac for these canonical smiles entries


# In[20]:


# notice that some compounds with correct canonical smiles have no entries in PubChem
NA = N[N.isnull().any(1)]


# In[21]:


NA


# In[30]:


#NA.to_csv('noCompoundPubChem.csv')


# In[22]:


NN = N.dropna()


# In[23]:


NN


# In[ ]:





# In[27]:


NNN=[] # df for metadata for canonical smiles entries in pubchem
for i, row in NN.iterrows():
    print(i)
    splitting0 = str(NN['ID'][i]).split('(') # split the ids into only numberds to psearch in pubchempy
    bx = re.sub('\ |\[|\]|\(|\)', '', splitting0[1])
    # search pubchempy to add metadata
    df4L = pcp.get_properties(['molecular_weight','molecular_formula', 'inchi'], bx, 'cid', as_dataframe=True)
    for ind, row in df4L.iterrows():
        NNN.append({
            'indexX':NN['indexX'][i],
            'ID':NN['ID'][i],
            'name':NN['name'][i],
            'non_isomeric_smiles':NN['non_isomeric_smiles'][i],
            'iupac':NN['iupac'][i],
            'synonyms':NN['synonyms'][i],
            'ID2':ind,
            'molecular_weight':df4L.iloc[0,1],
            'molecular_formula':df4L.iloc[0,0],
            'inchi':df4L.iloc[0,2]
    })
NN2=pd.DataFrame(NNN) # df for metadata for canonical smiles entries in pubchem


# In[28]:


NN2


# In[31]:


#REread the suspect list again with different variable name
suspect_listB = pd.read_csv("fully_curated_suspect_list_skeletonema.csv")


# In[ ]:





# In[234]:


#add a column for non-isomeric SMILES
#suspect_listB['nonIsomeric_SMILES_byRDKit'] = np.nan


# In[32]:


#add these non-isomeric SMILES into the original suspect list
for i, row in suspect_listB.iterrows():
    for j, row in NN2.iterrows():
        if i == NN2["indexX"][j]:
            suspect_listB.loc[i, 'nonIsomeric_SMILES_byRDKit'] = NN2['non_isomeric_smiles'][j]#add correct smiles


# In[237]:


####adding new columns of interest
#suspect_listB['iupac'] = np.nan
#suspect_listB['synonyms'] = np.nan
#suspect_listB['PubChemPY'] = np.nan


# In[33]:


suspect_listB


# In[34]:


NN2.synonyms = NN2.synonyms.astype(str) # synonyms as str


# In[35]:


dy.synonyms = dy.synonyms.astype(str)# synonyms as str


# In[36]:


dx.synonyms = dx.synonyms.astype(str)# synonyms as str


# In[ ]:





# In[262]:


dy


# In[37]:


####adding the above obtained information on NAME CHECKING into the original suspect list
for index, row in suspect_listB.iterrows():#get into the csv suspect list
    for ind, row in dy.iterrows():#compounds with correct names list
        if suspect_listB['SMILES'][index]==dy['isomeric_smiles'][ind] and index==dy['index'][ind]:
            suspect_listB.loc[index, 'iupac'] = dy['iupac'][ind]#add iupac
            suspect_listB.loc[index, 'synonyms'] = dy['synonyms'][ind]#add synonyms
            suspect_listB.loc[index, 'PubChemId'] = dy['CID'][ind]
            suspect_listB.loc[index, 'Monoisotopic_mass'] = dy['molecular_weight'][ind]
            suspect_listB.loc[index, 'Formula'] = dy['molecular_formula'][ind]
            suspect_listB.loc[index, 'InChI'] = dy['inchi'][ind]
            suspect_listB.loc[index, 'PubChemPY'] = 'checked; SAME name'#information on checked by pubchempy; SAME
    for indx, row in dx.iterrows():#compounds with incorrect names list
        if suspect_listB['SMILES'][index]==dx['isomeric_smiles'][indx] and index==dx['index'][indx]:
            suspect_listB.loc[index, 'Name'] = dx['name'][indx]# add correct name
            suspect_listB.loc[index, 'iupac'] = dx['iupac'][indx]#add iupac
            suspect_listB.loc[index, 'synonyms'] = dx['synonyms'][indx]#add synonyms
            suspect_listB.loc[index, 'PubChemId'] = dx['CID'][indx]
            suspect_listB.loc[index, 'Monoisotopic_mass'] = dx['molecular_weight'][indx]
            suspect_listB.loc[index, 'Formula'] = dx['molecular_formula'][indx]
            suspect_listB.loc[index, 'InChI'] = dx['inchi'][indx]
            suspect_listB.loc[index, 'PubChemPY'] = 'checked; NEW name'#information on checked by pubchempy; CHANGED
    for i, row in NN2.iterrows():#compounds with non-isomeric SMILES from rdkit
        if index==NN2['indexX'][i]:
            suspect_listB.loc[index, 'Name'] = NN2['name'][i]
            suspect_listB.loc[index, 'iupac'] = NN2['iupac'][i]
            suspect_listB.loc[index, 'synonyms'] = NN2['synonyms'][i]
            suspect_listB.loc[index, 'PubChemId'] = NN2['ID2'][i]
            suspect_listB.loc[index, 'Monoisotopic_mass'] = NN2['molecular_weight'][i]
            suspect_listB.loc[index, 'Formula'] = NN2['molecular_formula'][i]
            suspect_listB.loc[index, 'InChI'] = NN2['inchi'][i]
            suspect_listB.loc[index, 'PubChemPY'] = 'checked; NEW rdkit SMILES and name'


# In[38]:


suspect_listB


# In[39]:


####can write the suspect list with name checking done
suspect_listB.to_csv('NEWWWsmile_to_name.csv')


#  First Curation Achieved

# ## ADD MISSING STRUCTURAL INFORMATION AND WRONG SYNTAX SMILES

# In[40]:


#### Keep entries with no SMILES
bool_series = pd.notnull(suspect_listB["SMILES"])#check which indices don't have empty SMILES, Just like done above
woS_suspect_list = suspect_listB[-bool_series]#SMILES not present
woS_suspect_list


# In[41]:


####Remove entries with no names, these entries have no names or SMILES so no use in the suspect list
bool_series2 = pd.notnull(woS_suspect_list["Name"])
woS_suspect_list1 = woS_suspect_list[bool_series2]#names present; no empty entries
woS_suspect_list1


# In[42]:


noN_suspect_list = woS_suspect_list[-bool_series2]
noN_suspect_list


# In[82]:


#noN_suspect_list.to_csv('NaNnamesSL.csv')


# In[43]:


####some entries of names that are unknown removed here
woS_suspect_list2 = woS_suspect_list1[woS_suspect_list1.Name != 'unknown']
woS_suspect_list2


# In[44]:


noN_suspect_list1 = woS_suspect_list1[woS_suspect_list1.Name == 'unknown']
noN_suspect_list1


# In[83]:


#noN_suspect_list1.to_csv('UNKNOWNnameEntrySL.csv')


# In[45]:


d =[]#dataframe with names, to find the structural information

als = []#list of names, for which couldn't find any compound

dqs = []#list of compounds with more than one compound entry

####to find the structural information
for index, row in woS_suspect_list2.iterrows():
    nameSL = woS_suspect_list2["Name"][index]#name from suspect list
    y = pcp.get_compounds(nameSL, 'name')#compound based on name
    if y:#if its not empty
        df3 = pcp.get_properties(['isomeric_smiles', 'molecular_weight','molecular_formula', 'inchi','iupac_name'], nameSL, 'name', as_dataframe=True)#get properties relevant to the structure
        for compound in y:
            syn = compound.synonyms#for synonym list because cannot use the above get_properties function for this one
        if len(df3)==1:#if one name corresponds to one compound
            for ind, row in df3.iterrows():#propoerties dataframe
                d.append(
                {
                    'CID':ind,
                    'index':index,
                    'name':nameSL,
                    'isomeric_smiles':df3.iloc[0,2],
                    'molecular_weight':df3.iloc[0,1],
                    'molecular_formula':df3.iloc[0,0],
                    'inchi':df3.iloc[0,3],
                    'iupac':df3.iloc[0,4],
                    'synonyms':syn
                })
        else:
            dqs.append({
                'index':index,
                'name':nameSL
            })
    else:
        als.append({
            'index':index,
            'name':nameSL
        })
dn = pd.DataFrame(d)
qd = pd.DataFrame(dqs)
la = pd.DataFrame(als)


# In[46]:


dn


# In[48]:


la


# In[47]:


qd


# In[ ]:





# In[350]:


###### MAUNAUL ADDITION; CHANGE EVERYTIME ACCORDING TO THE ENTRY!!!!!

dfn = pcp.get_properties(['isomeric_smiles', 'molecular_weight','molecular_formula', 'inchi','iupac_name'], 'Threonic acid', 'name', as_dataframe=True)
dfn1 = dfn.loc[[439535]]
dfn1


# In[54]:


###### MAUNAUL ADDITION; CHANGE EVERYTIME ACCORDING TO THE ENTRY!!!!

y = pcp.get_compounds(439535, 'cid')
for compound in y:
    syn = compound.synonyms
syn = str(syn)
qd.at[1,'synonyms']=syn


# In[49]:


qd.loc[0,'CID']=439535
qd.loc[1,'CID']=74328989
qd.loc[2,'CID']=585998
qd.loc[3,'CID']=604
qd.loc[0,'MolecularFormula']='C4H8O5'
qd.loc[1,'MolecularFormula']='C32H46O16'
qd.loc[2,'MolecularFormula']='C20H18O7'
qd.loc[3,'MolecularFormula']='C6H12O7'
qd.loc[0,'MolecularWeight']=136.1
qd.loc[1,'MolecularWeight']=686.7
qd.loc[2,'MolecularWeight']=370.4
qd.loc[3,'MolecularWeight']=196.16
qd.loc[0,'IsomericSMILES']='C(C(C(C(=O)O)O)O)O'
qd.loc[1,'IsomericSMILES']='COC1=C(C=CC(=C1)CC(COC2C(C(C(C(O2)CO)O)O)O)C(CC3=CC(=C(C=C3)O)OC)COC4C(C(C(C(O4)CO)O)O)O)O'
qd.loc[2,'IsomericSMILES']='C1C2C(COC2OC3=CC4=C(C=C3)OCO4)C(O1)C5=CC6=C(C=C5)OCO6'
qd.loc[3,'IsomericSMILES']='C(C(C(C(C(C(=O)O)O)O)O)O)O'
qd.loc[0,'InChI']='InChI=1S/C4H8O5/c5-1-2(6)3(7)4(8)9/h2-3,5-7H,1H2,(H,8,9)'
qd.loc[1,'InChI']='InChI=1S/C32H46O16/c1-43-21-9-15(3-5-19(21)35)7-17(13-45-31-29(41)27(39)25(37)23(11-33)47-31)18(8-16-4-6-20(36)22(10-16)44-2)14-46-32-30(42)28(40)26(38)24(12-34)48-32/h3-6,9-10,17-18,23-42H,7-8,11-14H2,1-2H3'
qd.loc[2,'InChI']='InChI=1S/C20H18O7/c1-3-15-17(25-9-23-15)5-11(1)19-13-7-22-20(14(13)8-21-19)27-12-2-4-16-18(6-12)26-10-24-16/h1-6,13-14,19-20H,7-10H2'
qd.loc[3,'InChI']='InChI=1S/C6H12O7/c7-1-2(8)3(9)4(10)5(11)6(12)13/h2-5,7-11H,1H2,(H,12,13)'
qd.loc[0,'IUPACName']='2,3,4-trihydroxybutanoic acid'
qd.loc[1,'IUPACName']='2-[2,3-bis[(4-hydroxy-3-methoxyphenyl)methyl]-4-[3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxybutoxy]-6-(hydroxymethyl)oxane-3,4,5-triol'
qd.loc[2,'IUPACName']='5-[[3-(1,3-benzodioxol-5-yl)-1,3,3a,4,6,6a-hexahydrofuro[3,4-c]furan-6-yl]oxy]-1,3-benzodioxole'
qd.loc[3,'IUPACName']='2,3,4,5,6-pentahydroxyhexanoic acid'


# In[55]:


qd


# In[56]:


qd.to_csv('manual_onename_morecids.csv')


# In[385]:


### read the next


# In[386]:


##### I checked all the names for which pubchempy returns no compound, it is either because it doesnt have an entry for that compound OR the name was not properly written in the suspect list e.g: oxalic acid or oxanate etc... 
##### for such compounds, which had weird names and thats why pubchempy didnt detect any compound, I made another small database, which is as below


# In[57]:


la#### names that pubchempy couldnt detect because no entry OR wrong name syntax


# In[436]:


####Manual ADDITION, CHANGE EVERYTIME WITH EACH COMPOUND

#'Hexofuranose or hexose'
#x = pcp.get_compounds('Hexose', 'name')
#x
#dfn = pcp.get_properties(['xlogp', 'rotatable_bond_count','isomeric_smiles', 'molecular_weight','molecular_formula', 'inchi','iupac_name'], 'Hexose', 'name', as_dataframe=True)
#dfn


# In[75]:


y = pcp.get_compounds(206, 'cid')
for compound in y:
    syn = compound.synonyms
syn = str(syn)
la.at[4,'synonyms']=syn


# In[58]:


la.loc[4,'CID']=206
la.loc[5,'CID']=5793
la.loc[13,'CID']=14900
la.loc[15,'CID']=445580
la.loc[16,'CID']=5282743
la.loc[17,'CID']=6506600
la.loc[18,'CID']=21863047
la.loc[20,'CID']=5282367
la.loc[40,'CID']=135911925
la.loc[47,'CID']=971


# In[59]:


la.loc[4,'MolecularFormula']='C6H12O6'
la.loc[5,'MolecularFormula']='C6H12O6'
la.loc[13,'MolecularFormula']='C19H38O4'
la.loc[15,'MolecularFormula']='C22H32O2'
la.loc[16,'MolecularFormula']='C16H30O2'
la.loc[17,'MolecularFormula']='C16H26O2'
la.loc[18,'MolecularFormula']='C16H24O2'
la.loc[20,'MolecularFormula']='C31H40O2'
la.loc[40,'MolecularFormula']='C7H9N5O3'
la.loc[47,'MolecularFormula']='C2H2O4'


# In[60]:


la.loc[4,'MolecularWeight']=180.16
la.loc[5,'MolecularWeight']=180.16
la.loc[13,'MolecularWeight']=330.5
la.loc[15,'MolecularWeight']=328.5
la.loc[16,'MolecularWeight']=254.41
la.loc[17,'MolecularWeight']=250.38
la.loc[18,'MolecularWeight']=248.36
la.loc[20,'MolecularWeight']=444.6
la.loc[40,'MolecularWeight']=211.18
la.loc[47,'MolecularWeight']=90.03


# In[61]:


la.loc[4,'IsomericSMILES']='C(C1C(C(C(C(O1)O)O)O)O)O'
la.loc[5,'IsomericSMILES']='C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O'
la.loc[13,'IsomericSMILES']='CCCCCCCCCCCCCCCC(=O)OCC(CO)O'
la.loc[15,'IsomericSMILES']='CC/C=C\C/C=C\C/C=C\C/C=C\C/C=C\C/C=C\CCC(=O)O'
la.loc[16,'IsomericSMILES']='CCCCCCCCCCCCC/C=C/C(=O)O'
la.loc[17,'IsomericSMILES']='CCCCCCCCC/C=C/C=C/C=C/C(=O)O'
la.loc[18,'IsomericSMILES']='CCCCCCC/C=C/C=C/C=C/C=C/C(=O)O'
la.loc[20,'IsomericSMILES']='CC1=C(C(=O)C2=CC=CC=C2C1=O)C/C=C(\\C)/CC/C=C(\\C)/CC/C=C(\\C)/CCC=C(C)C'
la.loc[40,'IsomericSMILES']='C1C(NC2=C(N1)N=C(NC2=O)N)C(=O)O'
la.loc[47,'IsomericSMILES']='C(=O)(C(=O)O)O'


# In[62]:


la.loc[4,'InChI']='InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2'
la.loc[5,'InChI']='InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6?/m1/s1'
la.loc[13,'InChI']='InChI=1S/C19H38O4/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-19(22)23-17-18(21)16-20/h18,20-21H,2-17H2,1H3'
la.loc[15,'InChI']='InChI=1S/C22H32O2/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16-17-18-19-20-21-22(23)24/h3-4,6-7,9-10,12-13,15-16,18-19H,2,5,8,11,14,17,20-21H2,1H3,(H,23,24)/b4-3-,7-6-,10-9-,13-12-,16-15-,19-18-'
la.loc[16,'InChI']='InChI=1S/C16H30O2/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16(17)18/h14-15H,2-13H2,1H3,(H,17,18)/b15-14+'
la.loc[17,'InChI']='InChI=1S/C16H26O2/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16(17)18/h10-15H,2-9H2,1H3,(H,17,18)/b11-10+,13-12+,15-14+'
la.loc[18,'InChI']='InChI=1S/C16H24O2/c1-2-3-4-5-6-7-8-9-10-11-12-13-14-15-16(17)18/h8-15H,2-7H2,1H3,(H,17,18)/b9-8+,11-10+,13-12+,15-14+'
la.loc[20,'InChI']='InChI=1S/C31H40O2/c1-22(2)12-9-13-23(3)14-10-15-24(4)16-11-17-25(5)20-21-27-26(6)30(32)28-18-7-8-19-29(28)31(27)33/h7-8,12,14,16,18-20H,9-11,13,15,17,21H2,1-6H3/b23-14+,24-16+,25-20+'
la.loc[40,'InChI']='InChI=1S/C7H9N5O3/c8-7-11-4-3(5(13)12-7)10-2(1-9-4)6(14)15/h2,10H,1H2,(H,14,15)(H4,8,9,11,12,13)'
la.loc[47,'InChI']='InChI=1S/C2H2O4/c3-1(4)2(5)6/h(H,3,4)(H,5,6)'


# In[63]:


la.loc[4,'IUPACName']='6-(hydroxymethyl)oxane-2,3,4,5-tetrol'
la.loc[5,'IUPACName']='(3R,4S,5S,6R)-6-(hydroxymethyl)oxane-2,3,4,5-tetrol'
la.loc[13,'IUPACName']='2,3-dihydroxypropyl hexadecanoate'
la.loc[15,'IUPACName']='(4Z,7Z,10Z,13Z,16Z,19Z)-docosa-4,7,10,13,16,19-hexaenoic acid'
la.loc[16,'IUPACName']='(E)-hexadec-2-enoic acid'
la.loc[17,'IUPACName']='2E,4E,6E)-hexadeca-2,4,6-trienoic acid'
la.loc[18,'IUPACName']='(2E,4E,6E,8E)-hexadeca-2,4,6,8-tetraenoic acid'
la.loc[20,'IUPACName']='2-methyl-3-[(2E,6E,10E)-3,7,11,15-tetramethylhexadeca-2,6,10,14-tetraenyl]naphthalene-1,4-dione'
la.loc[40,'IUPACName']='2-amino-4-oxo-5,6,7,8-tetrahydro-3H-pteridine-6-carboxylic acid'
la.loc[47,'IUPACName']='oxalic acid'


# In[64]:


la.loc[4,'CorName']='Hexose'
la.loc[5,'CorName']='Glucose'
la.loc[13,'CorName']='glyceryl palmitate'
la.loc[15,'CorName']='Docosahexaenoic Acid'
la.loc[16,'CorName']='Hexadecenoic acid'
la.loc[17,'CorName']='Hexadecatrienoic acid'
la.loc[18,'CorName']='hexadecatetraenoic acid'
la.loc[20,'CorName']='Vitamin K2'
la.loc[40,'CorName']='6-Carboxy-5,6,7,8-tetrahydropterin'
la.loc[47,'CorName']='oxalic acid'


# In[76]:


la


# In[77]:


laN = la.dropna()


# In[78]:


laN


# In[79]:


laA = la[la.isnull().any(1)]


# In[80]:


laA


# In[81]:


#laA.to_csv('wrongNamesNotInPubchem.csv')


# In[468]:


laN


# In[469]:


dn


# In[84]:


O


# In[537]:


#O['non_isomeric_smiles'][6]


# In[85]:


dex=[]#compound present
dex1=[]
dex2=[]
for i, row in O.iterrows(): 
    for j, row in suspect_listB.iterrows():
        if O['indexX'][i] == j:
            nm = suspect_listB['Name'][j]
            y = pcp.get_compounds(nm, 'name')
            if y:
                df3 = pcp.get_properties(['isomeric_smiles', 'molecular_weight','molecular_formula', 'inchi','iupac_name'], nm, 'name', as_dataframe=True)#get properties relevant to the structure
                for compound in y:
                    syn = compound.synonyms
                if len(df3)==1:
                    for ind, row in df3.iterrows():#propoerties dataframe
                        dex.append(
                        {
                            'CID':ind,
                            'index':j,
                            'name':nm,
                            'isomeric_smiles':df3.iloc[0,2],
                            'molecular_weight':df3.iloc[0,1],
                            'molecular_formula':df3.iloc[0,0],
                            'inchi':df3.iloc[0,3],
                            'iupac':df3.iloc[0,4],
                            'synonyms':syn
                        })
                else:
                    dex1.append({
                        'index':j,
                        'name':nm
                    })
            else:
                dex2.append({
                    'index':j,
                    'name':nm
                })
dndex = pd.DataFrame(dex)#compound present
dndex2 = pd.DataFrame(dex1)#compound with more cios
dndex3 = pd.DataFrame(dex2)#compound not present


# In[86]:


dndex###NOTE:113 is already added earlier in NN2
##take only dndex[2] information (IUPAC and Synonyms)
# directly add to the susect list


# In[102]:


dndex['CID'][2]


# In[105]:


suspect_listB.loc[239,'iupac']= dndex['iupac'][2]
suspect_listB.loc[239,'PubChemId']= dndex['CID'][2]
suspect_listB.loc[239,'PubChemPY']='checked;SAMEname'


# In[106]:


y = pcp.get_compounds(90659181, 'cid')
for compound in y:
    syn = compound.synonyms
syn = str(syn)
suspect_listB.at[239,'synonyms']=syn


# In[108]:


dndex2
##take only dndex2[2] information (IUPAC and Synonyms)
# directly add to the susect list


# In[109]:


x = pcp.get_compounds('Adenosine-GDP-cobinamide', 'name')
x
dfn = pcp.get_properties(['xlogp', 'rotatable_bond_count','isomeric_smiles', 'molecular_weight','molecular_formula', 'inchi','iupac_name'], 'Adenosine-GDP-cobinamide', 'name', as_dataframe=True)
dfn


# In[539]:


dfn1['IsomericSMILES'][135922733]


# In[540]:


dfn1['IsomericSMILES'][135398566]


# In[110]:


suspect_listB.loc[931,'iupac']= dfn['IUPACName'][135922733]
suspect_listB.loc[931,'PubChemId']= 135922733
suspect_listB.loc[931,'PubChemPY']='checked;SAMEname'
y = pcp.get_compounds(135922733, 'cid')
for compound in y:
    syn = compound.synonyms
syn = str(syn)
suspect_listB.at[931,'synonyms']=syn


# In[91]:


dndex3


# In[92]:


#dndex3.to_csv('wrongSMILESNmaesnotinPubchem.csv')


# In[94]:


ffe3


# In[95]:


des=[]
des1=[]
des2=[]
for i, row in ffe3.iterrows(): 
    for j, row in suspect_listB.iterrows():
        if ffe3['indexX'][i] == j:
            nm = suspect_listB['Name'][j]
            y = pcp.get_compounds(nm, 'name')
            if y:
                df3 = pcp.get_properties(['isomeric_smiles', 'molecular_weight','molecular_formula', 'inchi','iupac_name'], nm, 'name', as_dataframe=True)#get properties relevant to the structure
                for compound in y:
                    syn = compound.synonyms
                if len(df3)==1:
                    for ind, row in df3.iterrows():#propoerties dataframe
                        des.append(
                        {
                            'CID':ind,
                            'index':j,
                            'name':nm,
                            'isomeric_smiles':df3.iloc[0,2],
                            'molecular_weight':df3.iloc[0,1],
                            'molecular_formula':df3.iloc[0,0],
                            'inchi':df3.iloc[0,3],
                            'iupac':df3.iloc[0,4],
                            'synonyms':syn
                        })
                else:
                    des1.append({
                        'index':j,
                        'name':nm
                    })
            else:
                des2.append({
                    'index':j,
                    'name':nm
                })
dndes = pd.DataFrame(des)
dndes2 = pd.DataFrame(des1)
dndes3 = pd.DataFrame(des2)


# In[96]:


dndes


# In[111]:


suspect_listB.loc[253,'iupac']= dndes['iupac'][0]
suspect_listB.loc[253,'SMILES']= dndes['isomeric_smiles'][0]
suspect_listB.loc[253,'PubChemId']= 134820202
suspect_listB.loc[253,'PubChemPY']='checked;SAMEname'
y = pcp.get_compounds(134820202, 'cid')
for compound in y:
    syn = compound.synonyms
syn = str(syn)
suspect_listB.at[253,'synonyms']=syn


# In[98]:


dndes2


# In[122]:


dndes3


# In[112]:


x2 = pcp.get_compounds('meso-Tartaric acid', 'name')
x2
dfn2 = pcp.get_properties(['xlogp', 'rotatable_bond_count','isomeric_smiles', 'molecular_weight','molecular_formula', 'inchi','iupac_name'], 'meso-Tartaric acid', 'name', as_dataframe=True)
dfn2


# In[505]:


dfn2['InChI'][447315]


# In[113]:


suspect_listB.loc[84,'PubChemId']=447315
suspect_listB.loc[84,'Formula']='C4H6O6'
suspect_listB.loc[84,'Monoisotopic_mass']=150.09
suspect_listB.loc[84,'SMILES']='C@@H]([C@@H](C(=O)O)O)(C(=O)O)O'
suspect_listB.loc[84,'InChI']='InChI=1S/C4H6O6/c5-1(3(7)8)2(6)4(9)10/h1-2,5-6H,(H,7,8)(H,9,10)/t1-,2+'
suspect_listB.loc[84,'iupac']='(2R,3S)-2,3-dihydroxybutanedioic acid'
suspect_listB.loc[84,'correct_Name']='meso-Tartaric acid'


# In[114]:


y = pcp.get_compounds(447315, 'cid')
for compound in y:
    syn = compound.synonyms
syn = str(syn)
suspect_listB.at[84,'synonyms']=syn


# In[125]:


dndes3NA = dndes3.drop([0])


# In[126]:


dndes3NA


# In[127]:


dndes3NA.to_csv('wrongSMILESNamesnoinpUbchem2.csv')


# In[100]:


for index, row in suspect_listB.iterrows():#get into the csv suspect list
    for ind, row in qd.iterrows():#weird names list
        if index==qd['index'][ind]:
            suspect_listB.loc[index, 'PubChemId'] = qd['CID'][ind]#add pubchemID
            suspect_listB.loc[index, 'SMILES'] = qd['IsomericSMILES'][ind]#add smiles
            suspect_listB.loc[index, 'Monoisotopic_mass'] = qd['MolecularWeight'][ind]#add mass
            suspect_listB.loc[index, 'Formula'] = qd['MolecularFormula'][ind]#add formula
            suspect_listB.loc[index, 'InChI'] = qd['InChI'][ind]#
            suspect_listB.loc[index, 'iupac'] = qd['IUPACName'][ind]#add iupac
            suspect_listB.at[index, 'synonyms'] = qd['synonyms'][ind]#add synonyms
            suspect_listB.loc[index, 'PubChemPY'] = 'checked; structure added'
    for indx, row in dn.iterrows():#normal names list
        if index==dn['index'][indx]:
            suspect_listB.loc[index, 'PubChemId'] = dn['CID'][indx]#add pubchemID
            suspect_listB.loc[index, 'SMILES'] = dn['isomeric_smiles'][indx]#add smiles
            suspect_listB.loc[index, 'Monoisotopic_mass'] = dn['molecular_weight'][indx]#add mass
            suspect_listB.loc[index, 'Formula'] = dn['molecular_formula'][indx]#add formula
            suspect_listB.loc[index, 'InChI'] = dn['inchi'][indx]#
            suspect_listB.loc[index, 'iupac'] = dn['iupac'][indx]#add iupac
            suspect_listB.at[index, 'synonyms'] = dn['synonyms'][indx]#add synonyms
            suspect_listB.loc[index, 'PubChemPY'] = 'checked; structure added'
    for i, row in laN.iterrows():#weird names list
        if index==laN['index'][i]:
            suspect_listB.loc[index, 'PubChemId'] = laN['CID'][i]#add pubchemID
            suspect_listB.loc[index, 'SMILES'] = laN['IsomericSMILES'][i]#add smiles
            suspect_listB.loc[index, 'Monoisotopic_mass'] = laN['MolecularWeight'][i]#add mass
            suspect_listB.loc[index, 'Formula'] = laN['MolecularFormula'][i]#add formula
            suspect_listB.loc[index, 'InChI'] = laN['InChI'][i]#
            suspect_listB.loc[index, 'iupac'] = laN['IUPACName'][i]#add iupac
            suspect_listB.at[index, 'synonyms'] = laN['synonyms'][i]#add synonyms
            suspect_listB.loc[index, 'PubChemPY'] = 'checked; structure added, new name added'
            suspect_listB.loc[index, 'correct_Name'] = laN['CorName'][i]


# In[115]:


suspect_listB


# In[116]:


#convert isotopic masses into float , replace , with .
xn = pd.to_numeric(suspect_listB['Monoisotopic_mass'].astype(str).str.replace(',','.'), errors='coerce')


# In[117]:


#convert isotopic masses into float , replace , with .
xn = pd.to_numeric(suspect_listB['Monoisotopic_mass'].astype(str).str.replace(',','.'), errors='coerce')
suspect_listB['numeric_mass'] = xn


# In[118]:


suspect_listB


# In[119]:


#write the curated suspect list as csv
suspect_listB.to_csv('NEW_suspectlist_Name_And_Structure_IsomericSMILES.csv')


# second curation is done

# ## Check Literature Review, Other Databases

# In[136]:


NA_R = pd.read_csv('noCompoundPubChem.csv')
NA_I =list(NA_R['indexX'])
NA


# In[132]:


noN_suspect_list #### Check the Paper (no other source of information)
#####tell maria about other sterols and lipids present in the paper


# In[133]:


noN_suspect_list1 ##check source MetaboLghts; no other information
##unknowns are unknown but vcan be kept in a separate table, there are few more compounds if they are worth adding?


# In[134]:


laA


# In[137]:


dndex3


# In[138]:


dndes3NA


# In[9]:





# In[10]:


#####read infor for these entries with their respectuve CIDs
datf = {'index': [30, 38, 101, 129, 170, 173, 174, 176, 177, 178, 181, 204, 210, 221, 230, 278, 265, 137, 155, 157, 235, 420, 572, 716, 795, 798, 821, 822, 871, 890, 916, 931], 'CID': [136176976, 136176976, 656504, 2723614, 867, 10316468, 54758663, 46878438, 40490600, 46878427, 65373, 70678568, 25245431, 90657287, 123131567, 25244846, 126961136, 11273547, 45480541, 46931117, 90659074, 135921686, 124201679, 126970647, 124201680, 124201686, 124201684, 124201682, 124201677, 124201675, 70678552, 135922733]}


# In[11]:


dataf = pd.DataFrame(datf)


# In[26]:


S_suspect_list = pd.read_csv('NEW_suspectlist_Name_And_Structure_IsomericSMILES.csv')


# In[13]:


S_suspect_list


# In[14]:


dbf = [] #dataframe with change names
dcf = [] #dataframe with same name
for index, row in S_suspect_list.iterrows():#get inside the suspect list with all SMILES
    for i, row in dataf.iterrows():
        if index == dataf['index'][i]:
            print(index)
            nm = S_suspect_list["Name"][index]
            cid = str(dataf['CID'][i])#read names from suspect list
            x = pcp.get_compounds(cid, 'cid')
            for compound in x:
                pac = compound.iupac_name
            for compound in x:
                syn = compound.synonyms
            synLow = {synL.lower() for synL in syn}#converting syn to low case
            if nm.lower() not in synLow :#if they have different names
                nmC = syn[0]
                df1L = pcp.get_properties(['molecular_weight','molecular_formula', 'inchi', 'isomeric_smiles'], cid, 'cid', as_dataframe=True)
                for indx, row in df1L.iterrows():
                    dbf.append(
                            {
                                'index':index,
                                'CID':indx,
                                'name':nmC,
                                'isomeric_smiles':df1L.iloc[0,2],
                                'iupac':pac,
                                'synonyms':syn,
                                'molecular_weight':df1L.iloc[0,1],
                                'molecular_formula':df1L.iloc[0,0],
                                'inchi':df1L.iloc[0,3]
                            })
            else:#if the names are same
                df2L = pcp.get_properties(['molecular_weight','molecular_formula', 'inchi', 'isomeric_smiles'], cid, 'cid', as_dataframe=True)
                for ind, row in df2L.iterrows():
                    dcf.append(
                            {
                                'index':index,
                                'CID':ind,
                                'name':nm,
                                'isomeric_smiles':df1L.iloc[0,2],
                                'iupac':pac,
                                'synonyms':syn,
                                'molecular_weight':df2L.iloc[0,1],
                                'molecular_formula':df2L.iloc[0,0],
                                'inchi':df2L.iloc[0,3]
                            })
dx = pd.DataFrame(dbf)
dy = pd.DataFrame(dcf)


# In[15]:


dx


# In[16]:


dy


# In[17]:


dx.synonyms = dx.synonyms.astype(str)
dy.synonyms = dy.synonyms.astype(str)


# In[27]:


####adding the above obtained information on NAME CHECKING into the original suspect list
for index, row in S_suspect_list.iterrows():#get into the csv suspect list
    for ind, row in dy.iterrows():#compounds with correct names list
        if index==dy['index'][ind]:
            S_suspect_list.loc[index, 'iupac'] = dy['iupac'][ind]#add iupac
            S_suspect_list.loc[index, 'synonyms'] = dy['synonyms'][ind]#add synonyms
            S_suspect_list.loc[index, 'PubChemId'] = dy['CID'][ind]
            S_suspect_list.loc[index, 'Monoisotopic_mass'] = dy['molecular_weight'][ind]
            S_suspect_list.loc[index, 'Formula'] = dy['molecular_formula'][ind]
            S_suspect_list.loc[index, 'InChI'] = dy['inchi'][ind]
            S_suspect_list.loc[index, 'PubChemPY'] = 'checked; SAME name'#information on checked by pubchempy; SAME
            S_suspect_list.loc[index, 'SMILES'] = dy['isomeric_smiles'][ind]
    for indx, row in dx.iterrows():#compounds with incorrect names list
        if index==dx['index'][indx]:
            S_suspect_list.loc[index, 'Name'] = dx['name'][indx]# add correct name
            S_suspect_list.loc[index, 'iupac'] = dx['iupac'][indx]#add iupac
            S_suspect_list.loc[index, 'synonyms'] = dx['synonyms'][indx]#add synonyms
            S_suspect_list.loc[index, 'PubChemId'] = dx['CID'][indx]
            S_suspect_list.loc[index, 'Monoisotopic_mass'] = dx['molecular_weight'][indx]
            S_suspect_list.loc[index, 'Formula'] = dx['molecular_formula'][indx]
            S_suspect_list.loc[index, 'InChI'] = dx['inchi'][indx]
            S_suspect_list.loc[index, 'PubChemPY'] = 'checked; NEW name'#information on checked by pubchempy; CHANGED
            S_suspect_list.loc[index, 'SMILES'] = dx['isomeric_smiles'][indx]


# In[28]:


S_suspect_list


# In[29]:


del S_suspect_list['Unnamed: 0.1']


# In[30]:


S_suspect_list


# In[32]:


#convert isotopic masses into float , replace , with .
xn = pd.to_numeric(S_suspect_list['Monoisotopic_mass'].astype(str).str.replace(',','.'), errors='coerce')
S_suspect_list['numeric_mass'] = xn


# In[33]:


S_suspect_list.to_csv("sussspext.csv")


# In[1]:





# In[3]:


sl = pd.read_csv("XXXXXXXX.csv")


# In[4]:





# In[52]:


scndlist = S_suspect_list.loc[[65, 68, 92, 98, 99, 100, 103, 105, 106, 108, 129, 143, 144, 146, 150, 188, 200, 218, 298, 300, 301, 371, 387, 391, 403, 147, 148, 242, 246, 249, 254, 255, 260, 267, 348, 411, 413, 457,458, 459, 376, 377, 378, 379, 380, 381, 382, 383, 495, 460, 461, 462, 463, 464, 465]]


# In[53]:


scndlist########list of compounds with missing info


# In[54]:


scndlist.to_csv("SecondSmallList.csv")


# In[36]:


S_suspect_listW = S_suspect_list.drop(S_suspect_list.index[[68, 92, 98, 99, 100, 103, 105, 108, 129, 143, 144, 146, 150, 188, 200, 218, 298, 300, 301, 371, 387, 391, 403, 147, 148, 242, 246, 249, 254, 255, 260, 267, 348, 411, 413, 457, 376, 377, 378, 379, 380, 381, 382, 383, 495, 460, 461, 462, 463, 464, 465]])


# In[37]:


S_suspect_listW


# In[38]:


S_suspect_listX = S_suspect_listW.drop_duplicates()


# In[39]:


S_suspect_listX


# In[44]:


#lst = [68, 92, 98, 99, 100, 103, 105, 108, 129, 143, 144, 146, 150, 188, 200, 218, 298, 300, 301, 371, 387, 391, 403, 147, 148, 242, 246, 249, 254, 255, 260, 267, 348, 411, 413, 457, 376, 377, 378, 379, 380, 381, 382, 383, 495, 460, 461, 462, 463, 464, 465]
#len(lst)


# In[46]:


S_suspect_listX.loc[84, 'SMILES']='[C@@H]([C@@H](C(=O)O)O)(C(=O)O)O'


# In[57]:


S_suspect_listZ = S_suspect_listX.drop(S_suspect_list.index[[65, 106, 458, 459]])


# In[58]:


S_suspect_listZ


# In[91]:


suspectlist = pd.read_csv('fully_curated_suspect_list_skeletonema.csv')


# In[ ]:


# on aciidental repklacement of monoisotopic masses with numeric masses
# thus readding monisoptoc masses 


# In[94]:


for i, row in suspectlist.iterrows():
    for j, row in S_suspect_listZ.iterrows():
        if i == j:
            print(j)
            S_suspect_listZ.loc[j, 'Monoisotopic_mass'] = suspectlist['Monoisotopic_mass'][i]


# In[96]:


#convert isotopic masses into float , replace , with .
mono = pd.to_numeric(S_suspect_listZ['Monoisotopic_mass'].astype(str).str.replace(',','.'), errors='coerce')
S_suspect_listZ['Monoisotopic_mass'] = mono


# In[97]:


S_suspect_listZ


# In[102]:


S_suspect_listZ=S_suspect_listZ.rename({'numeric_mass':'Molecular mass'}, axis=1)


# In[103]:


S_suspect_listZ.to_csv("CURATED_SUSPECT_LIST.csv")


# In[8]:


S_suspect_listZ = pd.read_csv("CURATED_SUSPECT_LIST.csv")


# In[9]:


S_suspect_listZ


# In[4]:


S_suspect_listZ.drop_duplicates()


# In[10]:


#del S_suspect_listZ['Unnamed: 0']
del S_suspect_listZ['Unnamed: 0.1']


# In[11]:


S_suspect_listZ.drop_duplicates()


# In[ ]:





# In[ ]:


#### Master Thesis


# In[43]:


#S_suspect_listNM = sl
#del S_suspect_listNM['Unnamed: 0']
#del S_suspect_listNM['Unnamed: 0.1']
#S_suspect_listNM = S_suspect_listNM.drop(['nonIsomeric_SMILES_byRDKit', 'PubChemPY', 'correct_Name'], axis = 1) 
#S_suspect_listNMX = S_suspect_listNM.drop_duplicates()
# S_suspect_listNMX = S_suspect_listNMX[["Unnamed: 0","Name", "Formula", "Species", "Monoisotopic_mass", "Molecular mass", "SMILES", "InChI", "iupac", "synonyms", "classification1", "classification2", "class", "super_class", "pathway", "ChEBIid", "KEGGid", "PubChemId", "source_database", "Source"]]
#n_sl  = S_suspect_listNMX.drop_duplicates(['Name', 'Formula', 'Species','Molecular mass', 'SMILES', 'InChI', 'iupac', 'classification1', 'classification2', 'class', 'super_class', 'pathway','ChEBIid', 'KEGGid', 'PubChemId', 'source_database', 'Source'])
#S_suspect_listNMX.to_csv("Curated_suspect_list_th.csv")
#S_suspect_listNMX.drop_duplicates(subset='SMILES', keep="last")
#S_suspect_l= S_suspect_listNMX.drop_duplicates(subset='PubChemId', keep="last")
#S_suspect_l.to_csv("checkkkkk.csv")
#df = S_suspect_listNMX
#suspect_list = pd.read_csv("CURATED_SUSPECT_LIST.csv")


# ## Additions from LOTUS

# In[14]:


#####radd infor for these entries with their respectuve CIDs
dats = {'Name': ['ectocarpene', 'Diatoxanthin', '6,9,12-hexadecatrienoic acid', '5-Hydroxyeicosapentaenoic acid', 'Echinenone', '6,9,12,15-hexadecatetraenoic acid', 'Canthaxanthin', '6,9,12,15-hexadecatetraenoic acid', 'Bovinic acid', 'Hexadeca-6,9,12-trienoic acid','13(s)-hpotre'], 'CID': [11744786, 6440986, 5282811, 6439678, 5281236, 5282834, 5281227, 54154481, 5280644, 54296589, 5497123], 'source_database':['CMNPD:CMNPD582', 'LOTUS:LTS0263026', 'LOTUS:LTS0230126', 'LOTUS: LTS0204626', 'LOTUS:LTS0173313', 'LOTUS:LTS0158644', 'LOTUS:LTS0154130', 'LOTUS:LTS0100262','LOTUS:LTS0042897', 'LOTUS:LTS0030494', 'LOTUS:LTS0009634']}
dats = pd.DataFrame(dats)


# In[36]:


sla_d = []
for i, row in dats.iterrows():
    print(i)
    cid = str(dats['CID'][i])
    print(cid)
    x = pcp.get_compounds(cid, 'cid')
    for compound in x:
        pac = compound.iupac_name
    for compound in x:
        syn = compound.synonyms
    df1L = pcp.get_properties(['molecular_weight','molecular_formula', 'inchi', 'isomeric_smiles'], cid, 'cid', as_dataframe=True)
    for indx, row in df1L.iterrows():
        sla_d.append(
            {
                'Name':dats['Name'][i],
                'Formula':df1L.iloc[0,0],
                'Species': 'S.costatum',
                'SMILES':df1L.iloc[0,2],
                'InChI':df1L.iloc[0,3],
                'Monoisotopic_mass':'',
                'classification1':'',
                'classification2':'',
                'ChEBIid':'',
                'KEGGid':'',
                'PubChemId':cid,
                'source_database':dats['source_database'][i],
                'Source':'',
                'nonIsomeric_SMILES_byRDKit':'',
                'iupac':pac,
                'synonyms':syn,
                'PubChemPY':'',
                'correct_Name':'',
                'Molecular mass':df1L.iloc[0,1],
                'class':'',
                'super_class':'',
                'pathway':'',

                            })
sla = pd.DataFrame(sla_d)


# In[37]:


sla


# In[51]:


new_sl = pd.concat([suspect_list, sla], ignore_index = True)


# In[74]:


new_sl


# In[ ]:





# In[72]:


es =[]
for i, row in new_sl.iterrows():
    try:
        print(i)
        c = pcp.Compound.from_cid(int(new_sl['PubChemId'][i]))
        new_sl.loc[i, 'Monoisotopic_mass'] = c.monoisotopic_mass
    except Exception as e:
        es.append([i, e])


# In[83]:


sl = new_sl.loc[new_sl.astype(str).drop_duplicates().index]


# In[84]:


S_suspect_listNM = sl.drop(['nonIsomeric_SMILES_byRDKit', 'PubChemPY', 'correct_Name'], axis = 1) 
S_suspect_listNMX = S_suspect_listNM[["Name", "Formula", "Species", "Monoisotopic_mass", "Molecular mass", "SMILES", "InChI", "iupac", "synonyms", "classification1", "classification2", "class", "super_class", "pathway", "ChEBIid", "KEGGid", "PubChemId", "source_database", "Source"]]


# In[85]:


S_suspect_listNMX


# In[86]:


S_suspect_listNMX.to_csv("Curated_suspect_list_th.csv")


# In[89]:


x = pd.read_csv("Curated_suspect_list_th.csv")
x


# In[90]:


del x['Unnamed: 0']


# In[94]:


x.to_csv("Curated_suspect_list_th.csv")# THESIS


# In[95]:


new_sl.to_csv("CURATED_SUSPECT_LIST.csv")# ACTUAL


# In[34]:


sl = pd.read_csv("XXXXXXXX.csv")


# In[35]:


sl


# In[19]:


sln = pd.read_csv("fully_curated_suspect_list_skeletonema.csv")


# In[20]:


sln


# In[58]:


n_sl


# In[59]:


for i, row in n_sl.iterrows():
    for j, row in sln.iterrows():
        if n_sl['Unnamed: 0'][i] == sln['Unnamed: 0'][j]:
            print(j)
            n_sl.loc[i, 'Name'] = sln['Name'][j]


# In[10]:





# In[23]:


df = pd.read_excel('/Users/mahnoorzulfiqar/Downloads/my_list.xlsx')


# In[24]:


df


# In[16]:


x = []
y = []
for i, row in sl.iterrows():
    for j, row in df.iterrows():
        if sl['Name'][i] == df['Name'][j]:
            print(i)
            x.append(i)
            y.append(j)


# In[60]:


n_sl


# In[62]:


es =[]
for i, row in n_sl.iterrows():
    try:
        print(i)
        c = pcp.Compound.from_cid(int(n_sl['PubChemId'][i]))
        n_sl.loc[i, 'Monoisotopic_mass'] = c.monoisotopic_mass
    except Exception as e:
        es.append([i, e])


# In[64]:


del n_sl['class']
del n_sl['super_class']
del n_sl['pathway']


# In[65]:


n_sl.to_csv('z_sl.csv')


# In[66]:


z = pd.read_csv('z_sl.csv')


# In[67]:


z


# In[68]:


del z['Unnamed: 0']
del z['Unnamed: 0.1']


# In[69]:


z.to_csv("Orig_suspect_list_thesis.csv")


# ## Add nucleosides, nucleotides and missing amino acids 

# In[13]:


cids = {'Names':['Aspartic Acid', 'Glutamine', 'Glycine', 'Methionine', 'Valine', 'cAMP', 'cGMP',
                'ATP', 'ADP', 'AMP', 'GTP', 'GDP', 'GMP', 'UTP', 'UDP', 'UMP', 'dTTP', 'dTDP', 
                 'dTMP', 'CTD', 'CDP', 'CMP', 'Adenosine', 'Cytidine', 'Uridine', 'Guanosine', 
                'Deoxyadenosine', 'Deoxycytidine', 'Deoxythymidine', 'Deoxyguanosine'],
        'CIDS':[5960, 5961, 750, 6137, 6287, 6076, 135398570, 5957, 6022, 6083, 135398633, 
                135398619, 135398631, 6133, 6031, 6030, 64968, 164628, 9700, 6176, 6132, 6131, 
                60961, 6175, 6029, 135398635, 13730, 13711, 5789, 135398592]}
list_remaining = pd.DataFrame(cids)
list_remaining


# In[14]:


suspect_listC = pd.read_csv("CURATED_SUSPECT_LIST.csv")
suspect_listC


# In[39]:


remaining_comp = []
for i, rows in list_remaining.iterrows():
    ids = int(list_remaining['CIDS'][i])
    c = pcp.Compound.from_cid(ids)
    remaining_comp.append(
        {
            'Name' : list_remaining['Names'][i],
            'Formula' : c.molecular_formula,
            'Species' : np.nan,
            'SMILES' : c.isomeric_smiles,
            'InChI' : c.inchi,
            'Monoisotopic_mass' : c.monoisotopic_mass,
            'classification1' : np.nan,
            'classification2' : np.nan,
            'CHEBIid' : np.nan,
            'KEGGid' : np.nan,
            'PubChemId' : ids,
            'source_database' : "PubChem",
            'Source' : np.nan,
            'nonIsomeric_SMILES_byRDKit' : np.nan,
            'iupac' : c.iupac_name,
            'synonyms' : c.synonyms,
            'PubChemPY' : "added with PubChemPY",
            'correct_Name' : np.nan,
            'Molecular mass' : c.molecular_weight,
            'class' : np.nan,
            'super_class' : np.nan,
            'pathway' : np.nan
        }
    )


# In[41]:


rem_comp = pd.DataFrame(remaining_comp)


# In[50]:


rem_comp


# In[52]:


del suspect_listC['Unnamed: 0']


# In[53]:


suspect_listC


# In[54]:


df_all_rows = pd.concat([suspect_listC, rem_comp], ignore_index=True)


# In[55]:


df_all_rows


# In[57]:


pd.DataFrame.to_csv(df_all_rows, "CURATED_SUSPECT_LIST_Remaining.csv")


# ## ADD SMILES FROM PUBCHEMPY

# In[1]:





# In[16]:


sl = pd.read_csv( "CURATED_SUSPECT_LIST_Remaining.csv")


# In[17]:


sl


# In[18]:


cids2 = {'Names':["(E,Z)-2,4-heptadienal", "(2E,4Z)-octa-2,4-dienal"],
        'CIDS':[11788274, 5352875]}
list_remainingsmiles = pd.DataFrame(cids2)
list_remainingsmiles


# In[22]:


remaining_comp2 = []
for i, rows in list_remainingsmiles.iterrows():
    ids = int(list_remainingsmiles['CIDS'][i])
    c = pcp.Compound.from_cid(ids)
    remaining_comp2.append(
        {
            'Name' : list_remainingsmiles['Names'][i],
            'Formula' : c.molecular_formula,
            'Species' : np.nan,
            'SMILES' : c.isomeric_smiles,
            'InChI' : c.inchi,
            'Monoisotopic_mass' : c.monoisotopic_mass,
            'classification1' : np.nan,
            'classification2' : np.nan,
            'CHEBIid' : np.nan,
            'KEGGid' : np.nan,
            'PubChemId' : ids,
            'source_database' : "PubChem",
            'Source' : np.nan,
            'nonIsomeric_SMILES_byRDKit' : np.nan,
            'iupac' : c.iupac_name,
            'synonyms' : c.synonyms,
            'PubChemPY' : "added with PubChemPY",
            'correct_Name' : np.nan,
            'Molecular mass' : c.molecular_weight,
            'class' : np.nan,
            'super_class' : np.nan,
            'pathway' : np.nan
        }
    )


# In[24]:


rem_comp2 = pd.DataFrame(remaining_comp2)
rem_comp2


# In[25]:





# In[28]:


del sl['Unnamed: 0']


# In[30]:


df_all_rows = pd.concat([sl, rem_comp2], ignore_index=True)


# In[31]:


df_all_rows


# In[32]:


pd.DataFrame.to_csv(df_all_rows, "Use_This_CURATED_SUSPECT_LIST.csv")


# In[ ]:





# ## Later Curation

# In[225]:


#read suspect list
sl = pd.read_csv("Use_This_CURATED_SUSPECT_LIST_with_classes.csv")


# In[116]:


# read compounds with no pubchem entry
ck = pd.read_csv("noCompoundPubChem.csv")


# In[117]:


# read compounds with weird SMILES syntax
wl = pd.read_csv("wrongSMILESPubChem.csv")


# In[166]:


# merge two types of SMILES
list_of_smiles = list(ck["non_isomeric_smiles"]) + list(wl['non_isomeric_smiles'])
#list_of_smiles.to_csv("list of smiles.txt")
list_of_smiles = np.unique(list_of_smiles)


# In[ ]:





# In[167]:


kek = []
kek0 = []
smi = []
db = []
for i in list_of_smiles:
    m = Chem.MolFromSmiles(i)
    Chem.Kekulize(m)
    kek_smiles = Chem.MolToSmiles(m,kekuleSmiles=True)
    kek0.append(kek_smiles)


# In[208]:


m = Chem.MolFromSmiles('CC(=CC=CC=C(C)C=CC=C(C)C(=O)CC12C(CC(CC1(O2)C)O)(C)C)C=CC=C(C)C=C=C3C(CC(CC3(C)O)O)(C)C')
Chem.Kekulize(m)
kek_smiles = Chem.MolToSmiles(m,kekuleSmiles=True)


# In[209]:


kek_smiles


# In[154]:


df1L


# In[171]:


x = []
y = []
ids = []
for i in kek0:
    try:
        comp = pcp.get_compounds(i, 'smiles')
        splitting0 = str(comp).split('(') # split the ids into only numberds to psearch in pubchempy
        bx = re.sub('\ |\[|\]|\(|\)', '', splitting0[1])
        if bx :
            print(bx)
            ids.append(bx)
        else:
            y.append(i)
    except:
        x.append(i)


# In[176]:


ids = [607588, 75144734, 122171, 25244683]


# In[239]:


db = []
for i in ids:
    comp = pcp.get_compounds(i, 'cid')
    for compound in comp:
        syn = compound.synonyms#synonyms from pubchempy 
        db.append({
            'Name':compound.iupac_name,
            'Formula':compound.molecular_formula,
            'Species':'',
            'SMILES': compound.isomeric_smiles,
            'InChI':compound.inchi,
            'Monoisotopic_mass':compound.monoisotopic_mass,
            'ChEBIid':'',
            'KEGGid':'',
            'PubChemId': i,
            'source_database':'',
            'Source':'',
            'nonIsomeric_SMILES_byRDKit':'',
            'iupac': compound.iupac_name,
            'synonyms':syn,
            'PubChemPY':'KekulizedSMILES with a PubChemid',
            'correct_Name': '',
            'Molecular mass':compound.molecular_weight,
            'subclass':'',
            'class':'',
            'superclass':''           
            
        }
    )


# In[240]:


data = pd.DataFrame(db)


# In[226]:


del sl['pathway']
del sl['CHEBIid']
del sl['Unnamed: 0']
del sl['Unnamed: 0.1']
del sl['classification1']
del sl['classification2']


# In[241]:


len(sl.columns)


# In[242]:


len(data.columns)


# In[243]:


df_all_rows = pd.concat([sl, data])
df_all_rows


# In[244]:


df_sl = df_all_rows.reset_index(drop=True)


# In[245]:


from pybatchclassyfire import *
import pandas as pd
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
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
import pandas.io.formats.style
from rdkit.Chem import PandasTools
import xlrd
import openpyxl
import statistics
#Import Libraries
import pubchempy as pcp
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import PandasTools
from rdkit.Chem import rdFMCS
import pandas.io.formats.style
import os
import re
import glob
import xlrd
import openpyxl
import statistics
import time
def isNaN(string):
    return string != string
def classification(frame):
    inchis = []
    for i, row in frame.iterrows():
        if not isNaN(frame['SMILES'][i]):
            try:
                InChI = Chem.MolToInchi(Chem.MolFromSmiles(frame["SMILES"][i]))
                InChIKey = Chem.inchi.InchiToInchiKey(InChI)
                inchis.append({
                    'index': i,
                    'smiles':frame["SMILES"][i],
                    'inchi': InChI,
                    'inchikey': InChIKey
                })
            except:
                pass
    inchis = pd.DataFrame(inchis)
    inchis = inchis.loc[-isNaN(inchis['inchikey'])]
    
    ## Retrieve ClassyFire classifications ##
    
    # This first step is done using inchikey and interrogation of the gnps classified structures
    gnps_proxy = True 
    url = "http://classyfire.wishartlab.com"
    proxy_url =  "https://gnps-classyfire.ucsd.edu"
    chunk_size = 1000
    sleep_interval = 12
    
    all_inchi_keys = list(inchis['inchikey'].drop_duplicates())

    resolved_ik_number_list = [0, 0]
    total_inchikey_number = len(all_inchi_keys)

    while True:
    
        start_time = time.time()
    
        print('%s inchikey to resolve' % total_inchikey_number )
        get_classifications_cf_mod(all_inchi_keys, par_level = 6)
    
        cleanse('all_json.json', 'all_json.json')
    
        with open("all_json.json") as tweetfile:
            jsondic = json.loads(tweetfile.read())

        df = json_normalize(jsondic)
        df = df.drop_duplicates( 'inchikey' )
        resolved_ik_number = len( df.drop_duplicates('inchikey').inchikey )
        resolved_ik_number_list.append( resolved_ik_number )
        print('%s resolved inchikeys' % resolved_ik_number )
        print("done in --- %s seconds ---" % (time.time() - start_time))
    
        if resolved_ik_number_list[-1] < resolved_ik_number_list[-2] or resolved_ik_number_list[-1] == resolved_ik_number_list[-3]:
            break
        cleanse('all_json.json', 'all_json_cleaned.json')
        
        with open("all_json_cleaned.json") as tweetfile:
            jsondic = json.loads(tweetfile.read())
            
    flattened_classified_json = json_normalize(jsondic)
    flattened_df = flattened_classified_json.drop_duplicates('inchikey')
    flattened_df['inchikey'] = flattened_df['inchikey'].str.replace(r'InChIKey=', '')
    df_merged = pd.merge(inchis, flattened_df, left_on='inchikey', right_on='inchikey', how='left')
    
    for p, rowp in df_merged.iterrows():
        for q, rowq in frame.iterrows():
            if df_merged["smiles_x"][p] is frame["SMILES"][q]:
                frame.loc[q, 'subclass'] = df_merged["subclass.name"][p]
                frame.loc[q, 'class'] = df_merged["class.name"][p]
                frame.loc[q, 'superclass'] = df_merged["superclass.name"][p]
                #frame.loc[q, 'Classification_Source'] = "ClassyFire"
    #frame.to_csv(input_dir, '/SIRIUS_combined.csv')
    return(frame)


# In[246]:


classification(df_sl)


# In[247]:


df_sl


# In[261]:


df_sl['SMILES'].isna().sum() # 112 entries have no class identified in classifier


# In[255]:


for i, rows in df_sl.iterrows():
    m = Chem.MolFromSmiles(df_sl['SMILES'][i],sanitize=False)
    if m is None:
        print('invalid SMILES')
    else:
        try:
            Chem.SanitizeMol(m)
        except:
            print('invalid chemistry')


# In[ ]:


# all SMILES structures are correct


# In[263]:


sl_new = df_sl.drop_duplicates(['SMILES','InChI','PubChemId', 'iupac', 'Name', 'Monoisotopic_mass'],keep= 'last')


# In[ ]:





# In[264]:


sl_new.to_csv("Use_This_CURATED_SUSPECT_LIST_with_classes_noDups.csv")


# In[ ]:





# In[ ]:


#read suspect list
sl = pd.read_csv("Use_This_CURATED_SUSPECT_LIST_with_classes.csv")


# In[ ]:





# In[175]:


#with open('list_of_smiles.smi', 'w') as f:
   # for line in y:
        #f.write(line)
        #f.write('\n')


# In[ ]:




