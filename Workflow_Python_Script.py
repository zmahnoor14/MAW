#!/usr/bin/env python
# coding: utf-8

# ## Load Python function script

# In[1]:


# import the function file
from Workflow_Python_Functions import *


# In[ ]:





# ## Define input directory and input table and Suspect List Path

# In[2]:


#Define input directory, keep all files in same directory and scripts so getwd works
input_dir = os.getcwd()+'/'
input_dir


# In[3]:


# read the suspect list
slistcsv = input_dir + "SkeletonemaSuspectListV1.csv"


# ## Workflow Functions

# ### Results post processing with insilico tools SIRIUS and MetFrag (Compound Databases)

# In[4]:


# SIRIUS Post processing part 2
sirius_postProc2(input_dir, 
                 input_tablecsv = input_dir + "input_table.csv")


# In[5]:


# Metfrag Post processing
metfrag_postproc(input_dir, 
                 input_tablecsv = input_dir + "input_table.csv")


# In[6]:


# combine results from different files; source = "all_insilico"
# note that source can be either all_insilico, SIRIUS, or MetFrag
# the result will be a combined results from different files either with SIRIUS or MetFrag
# with all_insilico as source, two files are genearted: one for SIRIUS and one for MetFrag
combine_insilico(input_dir, 
                 input_tablecsv = input_dir + "input_table.csv",
                Source = "all_insilico")


# In[7]:


# Curation of results from SIRIUS
sirius_curation(input_dir, 
                 siriuscsv = input_dir + "MetabolomicsResults/SIRIUS_combined.csv", 
                 sl = True)


# In[8]:


# Curation of results from MetFrag
metfrag_curation(input_dir, 
                 metfragcsv = input_dir + "MetabolomicsResults/MetFrag_combined.csv", 
                 sl = True)


# In[9]:


combineSM(input_dir, 
          metfragcsv = input_dir + 'MetabolomicsResults/metfrag_curated.csv', 
          siriuscsv = input_dir + 'MetabolomicsResults/sirius_curated.csv')


# ### Results post processing with Spectra package (Spectral Databases)

# In[10]:


# check each mzml file and each database csv result file; perform post processing
spec_postproc(input_dir, 
             Source = "all")


# In[4]:


# combine all spectral databases for each mzml file
combine_specdb(input_dir)


# In[5]:


# combine all spectral databases for all mzml file
combine_allspec(input_dir)


# In[6]:


# only keep good scoring spectral database results
scoring_spec(input_dir, 
             spec_file = input_dir + 'MetabolomicsResults/SD_post_processed_combined_results.csv')


# In[7]:


# suspect list VS spectal databases
# db(str): can be all, gnps, mbank, hmdb, gm(gnps, mbank), hg(hmdb and gnps), hm(hmdb and mbank) 
suspectListScreening(input_dir, 
                     slistcsv, 
                     SpectralDB_Results = input_dir + 'MetabolomicsResults/scoredSpecDB.csv', 
                     db = "all")


# In[8]:


specDB_Curation(input_dir, 
                combinedx = input_dir + 'MetabolomicsResults/SpecDBvsSL.csv', 
                sl = True, 
                db = "all")


# ### Final Results from all Sources

# In[9]:


combine_CuratedR(input_dir, 
                 combinedSDBs = input_dir + 'MetabolomicsResults/curatedSDB.csv', 
                 combinedSMs = input_dir + 'MetabolomicsResults/combinedSM.csv', 
                 data_type = "standards")


# In[10]:


checkSMILES_validity(input_dir, 
                     resultcsv = input_dir + 'MetabolomicsResults/final_curation_without_classes.csv')


# In[11]:


classification(input_dir, 
               resultcsv = input_dir + 'MetabolomicsResults/final_curation_with_validSMILES.csv')


# In[ ]:




