.. Metabolome Annotation Workflow (MAW\) documentation master file, created by
   Kohulan Rajan on Tue Aug  9 13:08:19 2022.

Welcome to Metabolome Annotation Workflow (MAW\)'s documentation!
=================================================================

.. image:: https://github.com/Kohulan/cheminf-jena-logos/blob/main/MAW/MAW.png?raw=true
  :width: 500
  :align: center

This repository hosts Metabolome Annotation Workflow (MAW). The workflow has been developed using the LCMS-2 dataset from a marine diatom Skeletonema marinoi. The workflow takes .mzML format data files as an input in R and performs spectral database dereplication using R Package Spectra and compound database dereplication using SIRIUS and MetFrag (with KEGG and PubChem). The results are saved as .csv files and are post processed in Python using RDKit and PubChemPy.The classification of the tentative candidates from the input data are classified using CANOPUS and ClassyFire, with a python client pybatchclassyfire for ClassyFire.

For comments, bug reports or feature ideas, please use github issues
or send an email to mahnoor.zulfiqar@uni-jena.de

Installation
============

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   modules

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
