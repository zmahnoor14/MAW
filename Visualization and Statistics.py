#!/usr/bin/env python
# coding: utf-8

# # Visualization

# In[157]:


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px
def isNaN(string):
    return string != string
import kaleido
import plotly.graph_objects as go
#Import libraries
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')


# In[202]:


cl = pd.read_csv("Final_Candidate_List.csv")
class_data = cl[['superclass', 'class', 'subclass']]
class_data = class_data.dropna()


# In[159]:


spclass = list(class_data['superclass']) # all superclasses
uniq_spclass = list(np.unique(list(class_data['superclass']))) # only unique super classes
#uniq_spc = [s for s in uniq_spclass if 'nan' not in s ] # only unique super classes with no NA values
len(uniq_spclass)
clss = list(class_data['class'])
uniq_class = list(np.unique(list(class_data['class'])))
#uniq_c = [s for s in uniq_class if 'nan' not in s ]
len(uniq_class)
sbclass = list(class_data['subclass'])
uniq_sbclass = list(np.unique(list(class_data['subclass'])))
#uniq_sbc = [s for s in uniq_sbclass if 'nan' not in s ]
len(uniq_sbclass)

#all characters
Names = ['Organic Compounds'] + uniq_spclass+uniq_class+uniq_sbclass

df = pd.DataFrame(Names)
df['values'] = ''
df['parents'] = ''

df = df.rename(columns={0: 'characters'})
df


# In[175]:


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
        


# In[198]:


data = dict(character = df['characters'], parents = df['parents'], values = df['values'])


# In[200]:


fig = px.sunburst(
    data,
    names='character',
    parents='parents',
    values='values',
)
fig.update_layout(margin = dict(t=0, l=0, r=0, b=0))
fig.show()


# In[247]:


fig.write_image("Sunburst_Chart.png")


# In[ ]:


### for Suspect List


# In[203]:


sl = pd.read_csv("/Users/mahnoorzulfiqar/OneDriveUNI/SuspectList/Use_This_CURATED_SUSPECT_LIST_with_classes_noDups.csv")
class_data = sl[['superclass', 'class', 'subclass']]
class_data = class_data.dropna()


# In[204]:


spclass = list(class_data['superclass']) # all superclasses
uniq_spclass = list(np.unique(list(class_data['superclass']))) # only unique super classes
#uniq_spc = [s for s in uniq_spclass if 'nan' not in s ] # only unique super classes with no NA values
len(uniq_spclass)
clss = list(class_data['class'])
uniq_class = list(np.unique(list(class_data['class'])))
#uniq_c = [s for s in uniq_class if 'nan' not in s ]
len(uniq_class)
sbclass = list(class_data['subclass'])
uniq_sbclass = list(np.unique(list(class_data['subclass'])))
#uniq_sbc = [s for s in uniq_sbclass if 'nan' not in s ]
len(uniq_sbclass)

#all characters
Names = ['Organic Compounds'] + uniq_spclass+uniq_class+uniq_sbclass

df = pd.DataFrame(Names)
df['values'] = ''
df['parents'] = ''

df = df.rename(columns={0: 'characters'})
df


# In[205]:


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
        


# In[206]:


data = dict(character = df['characters'], parents = df['parents'], values = df['values'])


# In[207]:


fig = px.sunburst(
    data,
    names='character',
    parents='parents',
    values='values',
)
fig.update_layout(margin = dict(t=0, l=0, r=0, b=0))
fig.show()


# In[ ]:





# ## Venn Diagram

# In[8]:


cl = pd.read_csv("Candidate_list.csv")
cl.columns


# In[12]:


sl = pd.read_csv("/Users/mahnoorzulfiqar/OneDriveUNI/SuspectList/Use_This_CURATED_SUSPECT_LIST_with_classes_noDups.csv")
sl


# In[18]:


LS22 = cl.loc[-isNaN(cl['SMILES_final'])]


# In[ ]:


cl.loc[cl['SMILES_final'] == 'S']


# In[19]:


# only in SIRIUS_formula
ls2 = list(LS22['Annotation_Source'])
len([x for x in ls2 if 'SuspectList' in x])


# In[98]:


venn2(subsets = (903, 771, 43), set_labels = ('Suspect List', 'Metabolomics Data'), set_colors=('skyblue', 'lightgreen'), alpha = 0.7);


# In[128]:


venn2(subsets = (102, 105, 44), set_labels = ('Suspect List Subclasses', 'Metabolomics Data Subclasses'), set_colors=('pink', 'yellow'), alpha = 0.7);


# In[136]:


venn2(subsets = (62, 75, 35), set_labels = ('Suspect List Classes', 'Metabolomics Data Classes'), set_colors=('pink', 'yellow'), alpha = 0.7);


# In[145]:


venn2(subsets = (12, 16, 10), set_labels = ('Suspect List Superclasses', 'Metabolomics Data Superclasses'), set_colors=('pink', 'yellow'), alpha = 0.7);


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# # Statistics

# In[58]:


import pandas as pd
import numpy as np


# In[70]:


#all results
results = pd.read_csv('/Users/mahnoorzulfiqar/OneDriveUNI/MZML/MetabolomicsResults/Candidate_list.csv')
results


# In[71]:


#remove duplicates based on features
resultsmzrt = results.drop_duplicates(['premz', 'rtmed'],keep= 'last')
resultsmzrt


# In[72]:


# only take entries with some SMILES Entry
resultsmzrt2 = resultsmzrt[-resultsmzrt['SMILES_final'].isna()]
resultsmzrt2


# In[74]:


#np.unique(list(resultsmzrt2['SMILES_final']))


# In[75]:


# only take entries without SMILES Entry
resultsmzrtna = resultsmzrt[resultsmzrt['SMILES_final'].isna()]
resultsmzrtna


# In[76]:


#remove duplicated SMILES entries from resultsmzrt
smiles = resultsmzrt.drop_duplicates(['SMILES_final'],keep= 'last')
smiles


# In[ ]:





# In[77]:


resultsx = pd.concat([resultsmzrtna, smiles])
resultsx# this is unduplicated final list of features


# In[146]:


resultsx.to_csv("Final_Candidate_List.csv")


# In[78]:


#resultsx2 = resultsx[resultsx['Annotation_Source'] != 'SIRIUS_Formula']
#resultsx2# with some database annotation


# In[79]:


# counting number of candidates in different annotation sources


# In[81]:


results2 = resultsx[-resultsx['Annotation_Source'].isna()]
results2#with a source of annotation


# In[82]:


results3 = results2[results2['Annotation_Source'] != 'SIRIUS_Formula']
results3# with some database annotation


# In[86]:


#len(np.unique(list(results2['SMILES_final'])))


# In[87]:


results2[results2['Annotation_Source'] == 'SIRIUS_Formula'] # with Formula annotation


# In[89]:


ls2 = list(results3['Annotation_Source'])
#ls2


# In[ ]:





# In[90]:


# only in spectral databases
len([x for x in ls2 if 'GNPS' in x or 'MassBank' in x or 'HMDB' in x])


# In[91]:


# only in SIRIUS and MetFrag databases
len([x for x in ls2 if 'PubChem' in x or 'KEGG' in x or 'SIRIUS' in x])


# In[92]:


# only in SIRIUS_formula
len([x for x in ls2 if 'SuspectList' in x])


# In[72]:


# counting classification


# In[137]:


classes = resultsx[-resultsx['superclass'].isna()]


# In[138]:


lsc = np.unique(list(classes['superclass']))
lsc
lscc = list(classes['superclass'])


# In[139]:


classesData = []
for i in lsc:
    num = lscc.count(i)
    classesData.append({
        'SuperClass':i,
        'Count':num
    })
classDB = pd.DataFrame(classesData)
classDB


# In[140]:


sl = pd.read_csv("/Users/mahnoorzulfiqar/OneDriveUNI/SuspectList/Use_This_CURATED_SUSPECT_LIST_with_classes_noDups.csv")


# In[141]:


classeSL = sl[-sl['superclass'].isna()]


# In[142]:


lscSL = np.unique(list(classeSL['superclass']))
lscSL
lsccSL = list(classeSL['superclass'])


# In[143]:


classesDataSL = []
for i in lscSL:
    num = lsccSL.count(i)
    classesDataSL.append({
        'SuperClass':i,
        'Count':num
    })
classDBSL = pd.DataFrame(classesDataSL)
classDBSL


# In[144]:


mergedStuff = pd.merge(classDBSL, classDB, on=['SuperClass'], how='inner')
mergedStuff


# In[109]:


classDB.sort_values(by = ['Count'], ascending=False).to_csv("topclasses.csv")


# In[182]:


# names


# In[183]:


namesdf = resultsx[-resultsx['Name'].isna()]
namesdfO = namesdf.sort_values(by = ['Occurence'], ascending=False)


# In[184]:


namesdfO[namesdfO['Occurence'] >= 3].to_csv("top_annotations.csv")

