# Whole Metabolome Annotation and Suspect List of _Skeletonema costatum_

## Suspect List Curation
Initial Suspect List contained compounds from KEGG, PubChem, MetaCyc, and BRENDAenzyme, LOTUS, CMNPD and Literature. The metadata included names, molecular formula, Species names, SMILES, InChI, Monoisotopic mass, different database IDs, and Source/ Reference. For the curation, the structural information from each entry is validated using ```RdKit``` and additional metadata is added via ```pubchempy```. 

First, entries with chemical structure notation such as SMILES were checked for the correct metadata using ```pubchempy```. The additional metadata added includes: IUPAC names, and syonyms. Another column indicates what curation has been done for each entry. This initila step resulted in three types of entries, one with correct metadta and names, second with incorrect name and metadata and third for which ```pubchempy``` returns error. The second type of entries were replaced with correct names and metdata. To handle the third type, the isomeric smiles were converted to canonical smiles uisng ```RdKit``` and then the metadata was added, but some entries had wrong SMILES syntax which caused the error for ```pubchempy``` and thus the any structure information was discarded. Some correct syntax SMILES had no record in PubChem. ## comment, if not in Suspect Lits, add them.


