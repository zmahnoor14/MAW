#!/usr/bin/env python
# coding: utf-8


# import the function file
from Workflow_Python_Functions import (
    os,
    spec_postproc,
    MCSS_for_SpecDB,
    sirius_postproc,
    MCSS_for_SIRIUS,
    CandidateSelection_SimilarityandIdentity,
    checkSMILES_validity,
    classification,
    Np_pathways,
    chemMN,
    gnpsMNvsgnpsMAW,
)


# Define input directory, keep all files in same directory and scripts so getwd works
input_dir = os.getwd() + "/data"
input_dir

spec_postproc(input_dir, Source="all")

MCSS_for_SpecDB(input_dir)

sirius_postproc(input_dir, exp_int=0.90, csi_score=-150)

MCSS_for_SIRIUS(input_dir)

CandidateSelection_SimilarityandIdentity(input_dir)

checkSMILES_validity(input_dir)

classification(input_dir)

Np_pathways(input_dir)

chemMN(input_dir)

gnpsMNvsgnpsMAW(input_dir)
