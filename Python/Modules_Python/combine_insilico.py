#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python
# coding: utf-8

#!/usr/bin/env python
# make executable in bash chmod +x PyRun

# Libraries

import sys

import pandas as pd
import numpy as np

import os
import glob
import re


def main():
    if len(sys.argv) != 3:
        print(
            "Usage python3 combine_insilico.py input_dir input_tablecsv.csv Source(optional, default= all_insilico)"
        )
    else:
        combine_insilico(sys.argv[1], sys.argv[2], sys.argv[3])


def combine_insilico(input_dir, input_tablecsv, Source="all_insilico"):
    def isNaN(string):
        return string != string

    """combine_insilico function combines the Sirius results from all
    result directories for each input mzml file. It does same for 
    Metfrag.

    Parameters:
    input_dir (str): This is the input directory where all the .mzML 
    files and their respective result directories are stored.
    
    input_table (str): This is the table in csv format (defined in R), 
    which stores a csv table containing columns "mzml_files", which 
    contains liat of all input files with their relative paths, second
    column is "ResultFileName" which is a list of the corresponding
    result relative directories to each mzml files. Lastly, "file_id", 
    contains a file directory. This table will be used to read the 
    Sirius and MetFrag result csv files
    
    Source (str): either "SIRIUS" or "MetFrag"

    Returns:
    
    dataframe: of combined SIRIUS/MetFrag results
    
    csv: stores the dataframe in a csv, named as 
    "input_dir/ResultFileName/MetabolomicsResults/SIRIUS_combined.csv" 
    OR/AND 
    "input_dir/ResultFileName/MetabolomicsResults/MetFrag_combined.csv"
    
    
    Usage:
    combine_insilico(input_dir = "/user/project/", 
    input_table = "/user/project/suspectlist.csv", Source = "SIRIUS")


    """

    input_table = pd.read_csv(input_tablecsv)
    # create a new directory to store all results /MetabolomicsResults/
    path = os.path.join(input_dir, "MetabolomicsResults")
    if not os.path.isdir(path):
        os.mkdir(path)
    # if Sirius results are to be combined
    if Source == "all_insilico" or Source == "SIRIUS":

        # store all files paths here
        all_files = []
        for n, row in input_table.iterrows():
            all_files.append(
                input_dir
                + input_table["ResultFileNames"][n].replace("./", "")
                + "/insilico/SiriusResults.csv"
            )

        # store all dataframes of the results here
        li = []

        for filename in all_files:
            df = pd.read_csv(filename, index_col=None, header=0)
            df["ResultFileNames"] = filename
            li.append(df)

        # join all resulst dataframe
        frame = pd.concat(li, axis=0, ignore_index=True)
        frame.to_csv(input_dir + "/MetabolomicsResults/SIRIUS_combined.csv")

    # if MetFrag results are to be combined
    if Source == "all_insilico" or Source == "MetFrag":

        # store all files paths here
        all_files = []
        for m, row in input_table.iterrows():
            all_files.append(
                input_dir
                + input_table["ResultFileNames"][m].replace("./", "")
                + "/insilico/MetFragResults.csv"
            )
        li = []

        for filename in all_files:
            df = pd.read_csv(filename, index_col=None, header=0)
            df["result_dir"] = filename
            li.append(df)

        frame = pd.concat(li, axis=0, ignore_index=True)
        frame.to_csv(input_dir + "MetabolomicsResults/MetFrag_combined.csv")


if __name__ == "__main__":
    main()
