#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# import the function file
from Workflow_Python_Functions import *


# In[ ]:


#Define input directory, keep all files in same directory and scripts so getwd works
input_dir = os.getwd()+"/data"
input_dir


# In[ ]:


spec_postproc(input_dir, Source = "all")


# In[ ]:


MCSS_for_SpecDB(input_dir, Source)


# In[ ]:


sirius_postproc(input_dir, exp_int = 0.90, csi_score = -150)


# In[ ]:


MCSS_for_SIRIUS(input_dir)    


# In[ ]:


CandidateSelection_SimilarityandIdentity(input_dir)


# In[ ]:


classification(input_dir, resultcsv)

