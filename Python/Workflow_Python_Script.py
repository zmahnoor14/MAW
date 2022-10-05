#!/usr/bin/env python
# coding: utf-8
# import the function file
from Workflow_Python_Functions import *  # Define input directory, keep all files in same directory and scripts so getwd works

input_dir = os.getcwd() + "/data"
input_dir

spec_postproc(input_dir, Source="all")

MCSS_for_SpecDB(input_dir, Source="all")

sirius_postproc(input_dir)

MCSS_for_SIRIUS(input_dir)

CandidateSelection_SimilarityandIdentity(input_dir, standards=False)

merge_all_results(input_dir)

# classification(input_dir, resultcsv = input_dir + "/final_candidates.csv")
