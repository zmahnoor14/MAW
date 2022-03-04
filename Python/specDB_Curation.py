#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python
# coding: utf-8

#!/usr/bin/env python
#make executable in bash chmod +x PyRun

# Libraries
import os
import glob
import re

import csv 
import time
import json
import sys
import pubchempy as pcp
import numpy as np

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem import PandasTools
import pubchempy as pcp


def specDB_Curation(input_dir, combinedx, sl = True):
    
    """specDB_Curation prioritizes in the following manner: gnps>
    mbank>suspectlist>hmdb

    Parameters:
    input_dir (str): This is the input directory where all the .mzML 
    files and their respective result directories are stored.
    
    combined: dataframe from either suspectListScreening function if
    sl = True OR from scoring_spec if sl = False
    
    Returns:
    dataframe: with curated Spectral DB results
    csv: "MetabolomicsResults/curatedSDB.csv"
    
    Usage:
    specDB_Curation(input_dir = "usr/project/",combinedx, sl = True)

    """
    
    def isNaN(string):
        return string != string

    # define scores again to remove the entries with lower scores
    def HMDB_Scoring(db, i):
        if db['HMDBmax_similarity'][i] >= 0.75 and db['HMDBintScore'][i] >= 0.50 and db['HMDBmzScore'][i] >= 0.50 and db['HQMatchingPeaks'][i]/db['hQueryTotalPeaks'][i] >= 0.50:
            return True
        else:
            return False
    def GNPS_Scoring(db, i):
        if db['GNPSmax_similarity'][i] >= 0.75 and db['GNPSintScore'][i] >= 0.50 and db['GNPSmzScore'][i] >= 0.50 and db['GQMatchingPeaks'][i]/db['gQueryTotalPeaks'][i] >= 0.50:
            return True
        else:
            return False
    def MB_Scoring(db, i):
        if db['MBmax_similarity'][i] >= 0.75 and db['MBintScore'][i] >= 0.50 and db['MBmzScore'][i] >= 0.50 and db['MQMatchingPeaks'][i]/db['mQueryTotalPeaks'][i] >= 0.50:
            return True
        else:
            return False
    
    
    combined = pd.read_csv(combinedx)
    
    
    # remove the similarity scores from low scoring candidates
    for i, row in combined.iterrows():
        if HMDB_Scoring(combined, i) or GNPS_Scoring(combined, i) or MB_Scoring(combined, i):
            print(i)
        else:
            combined['GNPSspectrumID'][i] = np.nan
            combined['MBspectrumID'][i] = np.nan
            combined['HMDBcompoundID'][i] = np.nan
    
    # if sl = True
    if sl:
    
        for i, row in combined.iterrows():
    
        ##### When there is an annotaion from all DBs #####
    
        #all entries with a high scoring annotation in all DBs,
            if not isNaN(combined['GNPSspectrumID'][i]) and not isNaN(combined['MBspectrumID'][i]) and not isNaN(combined['HMDBcompoundID'][i]):
        
                # entries with same candidate from all Spectral DBs
                if combined['tanimotoHG'][i] == 1.0 and combined['tanimotoGM'][i] == 1.0 and combined['tanimotoHM'][i] == 1.0:
                    combined.loc[i, 'Annotation'] = 'GNPS, HMDB, MassBank'
                    #entries with same candidate in suspect list, as in all Spectral DBs
                    if combined['GLname'][i] == combined['HLname'][i]== combined['MLname'][i]:
                        combined.loc[i, 'Annotation'] = 'GNPS, HMDB, MassBank, SuspectList'
                
                # same candidate from GNPS and HMDB        
                if combined['tanimotoHG'][i] == 1.0 and combined['tanimotoGM'][i] != 1.0 and combined['tanimotoHM'][i] != 1.0:
                    combined.loc[i, 'Annotation'] = 'GNPS, HMDB'
                    # if its present in Suspect List
                    if combined['GLname'][i] == combined['HLname'][i]:
                        combined.loc[i, 'Annotation'] = 'GNPS, HMDB, SuspectList'
        
                # same candidate from GNPS and MassBank        
                if combined['tanimotoHG'][i] != 1.0 and combined['tanimotoGM'][i] == 1.0 and combined['tanimotoHM'][i] != 1.0:
                    combined.loc[i, 'Annotation'] = 'GNPS, MassBank'
                    # if its present in Suspect List
                    if combined['GLname'][i] == combined['MLname'][i]:
                        combined.loc[i, 'Annotation'] = 'GNPS, MassBank, SuspectList'
                
                # same candidate from MassBank and HMDB        
                if combined['tanimotoHG'][i] != 1.0 and combined['tanimotoGM'][i] != 1.0 and combined['tanimotoHM'][i] == 1.0:
                    combined.loc[i, 'Annotation'] = 'GNPS, HMDB'
                    # if its present in Suspect List
                    if combined['GLname'][i] == combined['HLname'][i]:
                        combined.loc[i, 'Annotation'] = 'HMDB, MassBank, SuspectList'
                
                # only one database must be selected based on SuspectList annotation
                if combined['tanimotoHG'][i] != 1.0 and combined['tanimotoGM'][i] != 1.0 and combined['tanimotoHM'][i] != 1.0:
            
                    # only GNPS has SuspectList annotation
                    if not isNaN(combined['GLname'][i]):

                        combined.loc[i, 'Annotation'] = 'GNPS, SuspectList'
            
            
                    # only MassBank has SuspectList annotation
                    elif not isNaN(combined['MLname'][i]):
                        combined.loc[i, 'Annotation'] = 'MassBank, SuspectList'
            
            
                    # only HMDB has SuspectList annotation
                    elif not isNaN(combined['HLname'][i]):
                        combined.loc[i, 'Annotation'] = 'HMDB, SuspectList'
            
        
                    # all different annotations, take GNPS
                    else:
                        if not isNaN(combined['GNPSSMILES'][i]):
                            combined.loc[i, 'Annotation'] = 'GNPS'
                        else:
                            combined.loc[i, 'Annotation'] = 'MassBank'
            
        
            ##### When there is an annotation from two DBs #####
    
    
            # only GNPS and HMDB
            if not isNaN(combined['GNPSspectrumID'][i]) and isNaN(combined['MBspectrumID'][i]) and not isNaN(combined['HMDBcompoundID'][i]):

                if not isNaN(combined['GLname'][i]):
                    combined.loc[i, 'Annotation'] = 'GNPS, SuspectList'
                elif not isNaN(combined['HLname'][i]):
                    combined.loc[i, 'Annotation'] = 'HMDB, SuspectList'
                elif not isNaN(combined['GNPSSMILES'][i]):
                    combined.loc[i, 'Annotation'] = 'GNPS'
                else:
                    combined.loc[i, 'Annotation'] = 'HMDB'
        
    
            # only GNPS and MassBank
            if not isNaN(combined['GNPSspectrumID'][i]) and not isNaN(combined['MBspectrumID'][i]) and isNaN(combined['HMDBcompoundID'][i]):

                if not isNaN(combined['GLname'][i]):
                    combined.loc[i, 'Annotation'] = 'GNPS, SuspectList'
                elif not isNaN(combined['MLname'][i]):
                    combined.loc[i, 'Annotation'] = 'MassBank, SuspectList'
                elif not isNaN(combined['GNPSSMILES'][i]):
                    combined.loc[i, 'Annotation'] = 'GNPS'
                else:
                    combined.loc[i, 'Annotation'] = 'MassBank'
    
            # only MassBank and HMDB
    
            if isNaN(combined['GNPSspectrumID'][i]) and not isNaN(combined['MBspectrumID'][i]) and not isNaN(combined['HMDBcompoundID'][i]):

                if not isNaN(combined['MLname'][i]):
                    combined.loc[i, 'Annotation'] = 'MassBank, SuspectList'
                elif not isNaN(combined['HLname'][i]):
                    combined.loc[i, 'Annotation'] = 'HMDB, SuspectList'
                else:
                    combined.loc[i, 'Annotation'] = 'MassBank'
    
    
    
            ##### When there is an annotation from one DBs #####
    
    
            # only GNPS
            if not isNaN(combined['GNPSspectrumID'][i]) and isNaN(combined['MBspectrumID'][i]) and isNaN(combined['HMDBcompoundID'][i]):
        
                #If also SuspectList
                if not isNaN(combined['GLname'][i]):
                    combined.loc[i, 'Annotation'] = 'GNPS, SuspectList'
                elif not isNaN(combined['GNPSSMILES'][i]):
                    combined.loc[i, 'Annotation'] = 'GNPS'
        
            # only MassBank
            if isNaN(combined['GNPSspectrumID'][i]) and not isNaN(combined['MBspectrumID'][i]) and isNaN(combined['HMDBcompoundID'][i]):
                combined.loc[i, 'Annotation'] = 'MassBank'
                #If also SuspectList
                if not isNaN(combined['MLname'][i]):
                    combined.loc[i, 'Annotation'] = 'MassBank, SuspectList'
    
            # only HMDB
            if isNaN(combined['GNPSspectrumID'][i]) and isNaN(combined['MBspectrumID'][i]) and not isNaN(combined['HMDBcompoundID'][i]):
                combined.loc[i, 'Annotation'] = 'HMDB'
                #If also SuspectList
                if not isNaN(combined['HLname'][i]):
                    combined.loc[i, 'Annotation'] = 'HMDB, SuspectList'
            
    else:
        for i, row in combined.iterrows():
    
        ##### When there is an annotaion from all DBs #####
    
        #all entries with a high scoring annotation in all DBs,
            if not isNaN(combined['GNPSspectrumID'][i]) and not isNaN(combined['MBspectrumID'][i]) and not isNaN(combined['HMDBcompoundID'][i]):
        
                # entries with same candidate from all Spectral DBs
                if combined['tanimotoHG'][i] == 1.0 and combined['tanimotoGM'][i] == 1.0 and combined['tanimotoHM'][i] == 1.0:
                    combined.loc[i, 'Annotation'] = 'GNPS, HMDB, MassBank'
                
                # same candidate from GNPS and HMDB        
                if combined['tanimotoHG'][i] == 1.0 and combined['tanimotoGM'][i] != 1.0 and combined['tanimotoHM'][i] != 1.0:
                    combined.loc[i, 'Annotation'] = 'GNPS, HMDB'
        
                # same candidate from GNPS and MassBank        
                if combined['tanimotoHG'][i] != 1.0 and combined['tanimotoGM'][i] == 1.0 and combined['tanimotoHM'][i] != 1.0:
                    combined.loc[i, 'Annotation'] = 'GNPS, MassBank'
                
                # same candidate from MassBank and HMDB        
                if combined['tanimotoHG'][i] != 1.0 and combined['tanimotoGM'][i] != 1.0 and combined['tanimotoHM'][i] == 1.0:
                    combined.loc[i, 'Annotation'] = 'GNPS, HMDB'
                
                
            # all different annotations, take GNPS
            else:
                if not isNaN(combined['GNPSSMILES'][i]):
                    combined.loc[i, 'Annotation'] = 'GNPS'
                else:
                    combined.loc[i, 'Annotation'] = 'MassBank'
            
        
            ##### When there is an annotation from two DBs #####
    
    
            # only GNPS and HMDB
            if not isNaN(combined['GNPSspectrumID'][i]) and isNaN(combined['MBspectrumID'][i]) and not isNaN(combined['HMDBcompoundID'][i]):

                if not isNaN(combined['GNPSSMILES'][i]):
                    combined.loc[i, 'Annotation'] = 'GNPS'
                else:
                    combined.loc[i, 'Annotation'] = 'HMDB'
        
    
            # only GNPS and MassBank
            if not isNaN(combined['GNPSspectrumID'][i]) and not isNaN(combined['MBspectrumID'][i]) and isNaN(combined['HMDBcompoundID'][i]):

                if not isNaN(combined['GNPSSMILES'][i]):
                    combined.loc[i, 'Annotation'] = 'GNPS'
                else:
                    combined.loc[i, 'Annotation'] = 'MassBank'
    
            # only MassBank and HMDB
            if isNaN(combined['GNPSspectrumID'][i]) and not isNaN(combined['MBspectrumID'][i]) and not isNaN(combined['HMDBcompoundID'][i]):
                combined.loc[i, 'Annotation'] = 'MassBank'
    
    
    
            ##### When there is an annotation from one DBs #####
    
    
            # only GNPS
            if not isNaN(combined['GNPSspectrumID'][i]) and isNaN(combined['MBspectrumID'][i]) and isNaN(combined['HMDBcompoundID'][i]):
                if not isNaN(combined['GNPSSMILES'][i]):
                    combined.loc[i, 'Annotation'] = 'GNPS'
        
            # only MassBank
            if isNaN(combined['GNPSspectrumID'][i]) and not isNaN(combined['MBspectrumID'][i]) and isNaN(combined['HMDBcompoundID'][i]):
                combined.loc[i, 'Annotation'] = 'MassBank'
    
            # only HMDB
            if isNaN(combined['GNPSspectrumID'][i]) and isNaN(combined['MBspectrumID'][i]) and not isNaN(combined['HMDBcompoundID'][i]):
                combined.loc[i, 'Annotation'] = 'HMDB'
                
    combined.to_csv(input_dir + "MetabolomicsResults/curatedSDB.csv")
    return(combined)
    
specDB_Curation(sys.argv[1], sys.argv[2], sys.argv[3])

