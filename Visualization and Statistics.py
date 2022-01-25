#!/usr/bin/env python
# coding: utf-8

# # Visualization

# In[181]:


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px
def isNaN(string):
    return string != string
import kaleido
import plotly.graph_objects as go


# In[276]:


cl = pd.read_csv("Candidate_list.csv")


# In[277]:


class_data = cl[['superclass', 'class', 'subclass']]


# In[278]:


class_data = class_data.dropna()
class_data


# In[279]:


spclass = list(class_data['superclass'])
uniq_spclass = list(np.unique(list(class_data['superclass'])))
uniq_spc = [s for s in uniq_spclass if 'nan' not in s ]
clss = list(class_data['class'])
uniq_class = list(np.unique(list(class_data['class'])))
uniq_c = [s for s in uniq_class if 'nan' not in s ]
sbclass = list(class_data['subclass'])
uniq_sbclass = list(np.unique(list(class_data['subclass'])))
uniq_sbc = [s for s in uniq_sbclass if 'nan' not in s ]


# In[293]:


character = ['Organic Compounds'] + uniq_spc + uniq_c + uniq_sbc
values = []
parents = []
for i in character:
    if 'Organic Compounds' in i:
        values.append(0)
        parents.append('')
    elif i in uniq_spclass:
        values.append(spclass.count(i))
        parents.append('Organic Compounds')
    elif i in uniq_class:
        values.append(clss.count(i))
        clsp = class_data['superclass'][class_data[class_data['class'] == i].index.tolist()[0]]
        parents.append(clsp)
    elif i in uniq_sbclass:
        values.append(sbclass.count(i))
        sbclsp = class_data['class'][class_data[class_data['subclass'] == i].index.tolist()[0]]
        parents.append(sbclsp)


# In[294]:


data = dict(character = character, parents = parents, values = values)


# In[295]:


fig = px.sunburst(
    data,
    names='character',
    parents='parents',
    #values='values',
)
fig.show()


# In[296]:


data


# In[298]:


fig =go.Figure(go.Sunburst(
    labels = character, parents = parents #values = values
))
fig.update_layout(margin = dict(t=0, l=0, r=0, b=0))

fig.show()


# In[247]:


fig.write_image("Sunburst_Chart.png")


# # Statistics

# In[2]:


import pandas as pd
import numpy as np


# In[158]:


results = pd.read_csv('/Users/mahnoorzulfiqar/OneDriveUNI/MZML/MetabolomicsResults/Candidate_list.csv')
results


# In[160]:


resultsmzrt = results.drop_duplicates(['premz', 'rtmed'],keep= 'last')
resultsmzrt


# In[161]:


resultsmzrt2 = resultsmzrt[-resultsmzrt['SMILES_final'].isna()]
resultsmzrt2


# In[162]:


resultsmzrtna = resultsmzrt[resultsmzrt['SMILES_final'].isna()]
resultsmzrtna


# In[163]:


#remove duplicated SMILES entries from resultsmzrt
smiles = resultsmzrt.drop_duplicates(['SMILES_final'],keep= 'last')
smiles


# In[ ]:





# In[164]:


resultsx = pd.concat([resultsmzrtna, smiles])
resultsx


# In[145]:


# counting number of candidates in different annotation sources


# In[166]:


results2 = resultsx[-resultsx['Annotation_Source'].isna()]
results2#with a source of annotation


# In[197]:


list(results3[results3['Annotation_Source'].str.contains("SuspectList")==True]['Name'])


# In[169]:


results3 = results2[results2['Annotation_Source'] != 'SIRIUS_Formula']
results3# with some database annotation


# In[ ]:





# In[168]:


results2[results2['Annotation_Source'] == 'SIRIUS_Formula'] # with Formula annotation


# In[171]:


ls2 = list(results3['Annotation_Source'])
ls2


# In[ ]:





# In[172]:


# only in spectral databases
len([x for x in ls2 if 'GNPS' in x or 'MassBank' in x or 'HMDB' in x])


# In[173]:


# only in SIRIUS and MetFrag databases
len([x for x in ls2 if 'PubChem' in x or 'KEGG' in x or 'SIRIUS' in x])


# In[174]:


# only in SIRIUS_formula
len([x for x in ls2 if 'SuspectList' in x])


# In[72]:


# counting classification


# In[176]:


classes = resultsx[-resultsx['superclass'].isna()]


# In[185]:


lsc = np.unique(list(classes['superclass']))
lsc
lscc = list(classes['superclass'])


# In[190]:


classesData = []
for i in lsc:
    num = lscc.count(i)
    classesData.append({
        'SuperClass':i,
        'Count':num
    })
classDB = pd.DataFrame(classesData)
classDB


# In[192]:


classDB.sort_values(by = ['Count'], ascending=False).to_csv("topclasses.csv")


# In[182]:


# names


# In[183]:


namesdf = resultsx[-resultsx['Name'].isna()]
namesdfO = namesdf.sort_values(by = ['Occurence'], ascending=False)


# In[184]:


namesdfO[namesdfO['Occurence'] >= 3].to_csv("top_annotations.csv")


# In[ ]:


#Total features = 1619
#Any annotation from databases = 773
#Spectral database = 165

#Formula = 755
#SIRIUS or Metfrag =  608
#Suspect list = 158

#Superclass = 1027
#Class = 1005
#Subclass = 845

