#!/usr/bin/env python
# coding: utf-8

# ## Load Python function script

<<<<<<< HEAD
# In[10]:
=======
# In[1]:
>>>>>>> 25c6491 (cleaned directory)


# import the function file
from Workflow_Python_Functions import *


# In[ ]:





# ## Define input directory and input table

<<<<<<< HEAD
# In[11]:
=======
# In[2]:
>>>>>>> 25c6491 (cleaned directory)


#Define input directory, keep all files in same directory and scripts so getwd works
input_dir = os.getcwd()+'/'
input_dir


<<<<<<< HEAD
# In[12]:


# read input_table.csv generated in the R workflow
input_table = pd.read_csv(input_dir + "input_table.csv")
input_table


# ## Suspect List Path

# In[9]:


# read the suspect list
slistcsv = "/Users/mahnoorzulfiqar/OneDriveUNI/SuspectList/final_SUSPECT_LIST.csv"
=======
# In[3]:


# read input_table.csv generated in the R workflow
input_table = pd.read_csv(input_dir + "input_table.csv")
input_table


# ## Suspect List Path

# In[4]:


# read the suspect list
slistcsv = "/Users/mahnoorzulfiqar/OneDriveUNI/SuspectList/final_SUSPECT_LIST.csv"


# In[ ]:





# ## insilico folder results post processing (SIRIUS and MetFrag Results)

# In[11]:


# SIRIUS result postprocessing part 2 (part 1 in R)
# returns a csv file called input_dir/result_dir/insilico/SiriusResults.csv
sirius_postProc2(input_dir, input_table, slistcsv ,sl = True) 


# In[12]:


# MetFrag result postprocessing
# returns a csv file called input_dir/result_dir/insilico/MetFragResults.csv

# ignore the warnibgs from rdkit
metfrag_postproc(input_dir, input_table, slistcsv ,sl = True)


# In[13]:


# combines the different SiriusResults.csv for different mzml files
# returns input_dir/MetabolomicsResults/SIRIUS_combined.csv
sir_comb_results = combine_insilico(input_dir, input_table, Source = "SIRIUS")


# In[14]:


# combines the different MetFragResults.csv for different mzml files
# returns input_dir/MetabolomicsResults/MetFrag_combined.csv
met_comb_results = combine_insilico(input_dir, input_table, Source = "MetFrag")


# In[15]:


# curation of SIRIUS_combined.csv
# returns input_dir/MetabolomicsResults/sirius_curated.csv
s_cur = sirius_curation(input_dir, 
                siriuscsv = (input_dir+"/MetabolomicsResults/SIRIUS_combined.csv"), 
                sl = True)


# In[16]:


# curation of MetFrag_combined.csv
# returns input_dir/MetabolomicsResults/metfrag_curated.csv
m_cur = metfrag_curation(input_dir, 
                 metfragcsv = (input_dir+"/MetabolomicsResults/MetFrag_combined.csv"), 
                 sl = True)


# In[17]:


# combines the results from metfrag_curated.csv and sirius_curated.csv
# returns input_dir/MetabolomicsResults/curatedSM.csv
com_sm = combineSM(input_dir, m_cur, s_cur)


# In[ ]:





# # Spectral DB Post processing

# In[12]:


# reads individual spectral_dereplication/gnps.csv, spectral_dereplication/hmdb.csv, spectral_dereplication/mbank.csv and processes
# returns processed gnps, hmdb and mbank results in the individual mzml result directory
# ignore rdkit warning/error msgs
spec_postproc(input_dir, Source = "all")


# In[ ]:





# In[13]:


# Combine the gnps, mbank and hmdb into one csv file for each mzml input file
#returns input_dir/result_dir/spectral_dereplication/mergedR.csv
combine_specdb(input_dir)


# In[14]:


# Combine the mergedR.csv from different input files into one csv file
# returns input_dir/MetabolomicsResults/SD_post_processed_combined_results.csv
combine_allspec(input_dir)


# In[ ]:





# In[15]:


# only keeps the high scoring results
# returns input_dir/MetabolomicsResults/SD_post_processed_combined_results.csv
scoring_spec(input_dir, combined = input_dir+"MetabolomicsResults/combinedSpecDB.csv")


# In[ ]:





# In[16]:


# runs the spectral db results against the suspect list entries
# returns input_dir/MetabolomicsResults/SpecDBvsSL.csv
suspectListScreening(input_dir, slistcsv, 
                     SpectralDB_Results = input_dir +"MetabolomicsResults/combinedSpecDB.csv")


# In[ ]:





# In[17]:


# further curates the results from spectral dereplication
# # returns input_dir/MetabolomicsResults/curatedSDB.csv

specDB_Curation(input_dir, 
                combinedx = input_dir +"MetabolomicsResults/SpecDBvsSL.csv", 
                sl = True)


# In[ ]:





# ## Combine results from insilico and spectral_dereplication

# In[ ]:





# In[5]:


# combines the results from curated insilico and dereplication results
# returns input_dir/MetabolomicsResults/final_curation_without_classes.csv
combine_CuratedR(input_dir,
                 combinedSDBs = input_dir + "MetabolomicsResults/curatedSDB.csv", 
                 combinedSMs = input_dir + "/MetabolomicsResults/combinedSM.csv")


# In[2]:


# check validity of smiles
# returns input_dir/MetabolomicsResults/final_curation_with_validSMILES.csv

checkSMILES_validity(input_dir, resultcsv = input_dir + "MetabolomicsResults/final_curation_without_classes.csv")


# In[5]:


# performs calssification for compounds withno classification performed by SIRIUS/Canopus
# returns input_dir/MetabolomicsResults/final_curationList.csv
final_list = classification(input_dir, resultcsv = input_dir + "MetabolomicsResults/final_curation_with_validSMILES.csv")
>>>>>>> 25c6491 (cleaned directory)


# In[ ]:





# ## insilico folder results post processing (SIRIUS and MetFrag Results)

# In[11]:


# SIRIUS result postprocessing part 2 (part 1 in R)
# returns a csv file called input_dir/result_dir/insilico/SiriusResults.csv
sirius_postProc2(input_dir, input_table, slistcsv ,sl = True) 


# In[12]:


# MetFrag result postprocessing
# returns a csv file called input_dir/result_dir/insilico/MetFragResults.csv

# ignore the warnibgs from rdkit
metfrag_postproc(input_dir, input_table, slistcsv ,sl = True)


# In[13]:


# combines the different SiriusResults.csv for different mzml files
# returns input_dir/MetabolomicsResults/SIRIUS_combined.csv
sir_comb_results = combine_insilico(input_dir, input_table, Source = "SIRIUS")


# In[14]:


# combines the different MetFragResults.csv for different mzml files
# returns input_dir/MetabolomicsResults/MetFrag_combined.csv
met_comb_results = combine_insilico(input_dir, input_table, Source = "MetFrag")


# In[15]:


# curation of SIRIUS_combined.csv
# returns input_dir/MetabolomicsResults/sirius_curated.csv
s_cur = sirius_curation(input_dir, 
                siriuscsv = (input_dir+"/MetabolomicsResults/SIRIUS_combined.csv"), 
                sl = True)


# In[16]:


# curation of MetFrag_combined.csv
# returns input_dir/MetabolomicsResults/metfrag_curated.csv
m_cur = metfrag_curation(input_dir, 
                 metfragcsv = (input_dir+"/MetabolomicsResults/MetFrag_combined.csv"), 
                 sl = True)


# In[17]:


# combines the results from metfrag_curated.csv and sirius_curated.csv
# returns input_dir/MetabolomicsResults/curatedSM.csv
com_sm = combineSM(input_dir, m_cur, s_cur)


# In[ ]:





# # Spectral DB Post processing

# In[12]:


# reads individual spectral_dereplication/gnps.csv, spectral_dereplication/hmdb.csv, spectral_dereplication/mbank.csv and processes
# returns processed gnps, hmdb and mbank results in the individual mzml result directory
# ignore rdkit warning/error msgs
spec_postproc(input_dir, Source = "all")


# In[ ]:





# In[13]:


# Combine the gnps, mbank and hmdb into one csv file for each mzml input file
#returns input_dir/result_dir/spectral_dereplication/mergedR.csv
combine_specdb(input_dir)


# In[14]:


# Combine the mergedR.csv from different input files into one csv file
# returns input_dir/MetabolomicsResults/SD_post_processed_combined_results.csv
combine_allspec(input_dir)


# In[ ]:





# In[15]:


# only keeps the high scoring results
# returns input_dir/MetabolomicsResults/SD_post_processed_combined_results.csv
scoring_spec(input_dir, combined = input_dir+"MetabolomicsResults/combinedSpecDB.csv")


# In[ ]:





# In[16]:


# runs the spectral db results against the suspect list entries
# returns input_dir/MetabolomicsResults/SpecDBvsSL.csv
suspectListScreening(input_dir, slistcsv, 
                     SpectralDB_Results = input_dir +"MetabolomicsResults/combinedSpecDB.csv")


# In[ ]:





# In[17]:


# further curates the results from spectral dereplication
# # returns input_dir/MetabolomicsResults/curatedSDB.csv

specDB_Curation(input_dir, 
                combinedx = input_dir +"MetabolomicsResults/SpecDBvsSL.csv", 
                sl = True)


# In[ ]:





# ## Combine results from insilico and spectral_dereplication

# In[ ]:





# In[5]:


# combines the results from curated insilico and dereplication results
# returns input_dir/MetabolomicsResults/final_curation_without_classes.csv
combine_CuratedR(input_dir,
                 combinedSDBs = input_dir + "MetabolomicsResults/curatedSDB.csv", 
                 combinedSMs = input_dir + "/MetabolomicsResults/combinedSM.csv")


# In[5]:


# check validity of smiles
# returns input_dir/MetabolomicsResults/final_curation_with_validSMILES.csv

checkSMILES_validity(input_dir, resultcsv = input_dir + "MetabolomicsResults/final_curation_without_classes.csv")


# In[5]:


# performs calssification for compounds withno classification performed by SIRIUS/Canopus
# returns input_dir/MetabolomicsResults/final_curationList.csv
classification(input_dir, resultcsv = input_dir + "MetabolomicsResults/final_curation_with_validSMILES.csv")

