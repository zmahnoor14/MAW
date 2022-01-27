# Whole Metabolome Annotation and Suspect List of _Skeletonema costatum_

## Suspect List Curation
Initial Suspect List contained compounds from KEGG, PubChem, MetaCyc, and BRENDAenzyme, LOTUS, CMNPD and Literature. The metadata included names, molecular formula, Species names, SMILES, InChI, Monoisotopic mass, different database IDs, and Sources/ References. For the curation, the structural information from each entry is validated using ```RdKit``` and additional metadata is added via ```pubchempy```. 

First, entries with chemical structure notation such as SMILES were checked for the correct metadata using ```pubchempy```. The additional metadata added includes: IUPAC names, and syonyms. Another column indicates what curation has been done for each entry. This initial step resulted in three types of entries, one with correct metadata and names, second with incorrect name and metadata, and third for which ```pubchempy``` returns error. The second type of entries were replaced with correct names and metdata. To handle the third type, the isomeric smiles were converted to canonical smiles uisng ```RdKit``` and then the metadata was added, but some entries had wrong or non-standardized SMILES syntax which caused the error for ```pubchempy```. Data with wrong SMILES syntax were discarded. Some correct syntax SMILES had no record in PubChem. Such SMILES were kekulized and checked in PubChem with metadata.

For entries with names and other metadata but no structure notation, the names were searched in PubChem for isomeric SMILES and other missing metadata. Some enteries had non-conventional names such as Disccharide or 18:1 fatty alcohol. The names were manually checked and changed to conventional names of the compound. for others such as 18:1 fatty alochol, there was not much information to find the right name or strcuture. Such entries were also discarded. For some earlier entries with wrong or non-standardized SMILES syntax, the names were used to extract compounds from PubChem.

Lastly, the SMILES were checked again for syntax errors or invalid chemistry by ```RdKit``` and the duplicates were removed. The final list contains 903 entries with the following columns: Name, Formula, Species, SMILES, InChI, Monoisotopic_mass, PubChemId, source_database, Source, nonIsomeric_SMILES_byRDKit, iupac, synonyms, PubChemPY, correct_Name, Molecular mass, subclass, class', 'superclass.

## Candidate List Curation
Since top candidates for each feature are obtained from different database sources, hence when different sources give different candidates, then a prioritization scheme is requried. The idenification od reference standards has shown that among spectral databases, GNPS and MassBank are prioritized HMDB and among the compound databases, results from SIRIUS are prioritized over MetFrag. In terms of overall prioritizaton, if the candidate has a match in Suspect list, it is given prioritization. In general the scheme is: GNPS > Massbank > SIRIUS > HMDB > MetFrag(PubChem) > MetFrag(KEGG)

### MetFrag Curation
1. for compounds with only Pubchem, or only KEGG, will remain as it is.
2. if there is KEGG and PubChem, calculate a tanimoto score, if it is equal to 1, keep PubChem entry but add KEGG to source add well.
3. if there is less than 1 tanimoto, Check whether the entry has a suspect list compound, add the entry with Suspct list , if no suspectlist, then add PubChem

### SIRIUS Curation
important scores to consider are

SIRIUSscore
CSIFingerIDscore
exp_int >= 0.70
SL_comp (if present, give priority and mention in the annotation)

1. If the explained intensity is greater than 0.70 and there is no suspect list entry
2. If the explained intensity is greater than 0.70 and there is an entry from suspect list
3. if the intensity is less thna 0.70 but it still is similar to an entry in Suspect list

### Candidate Selection from Compound Databases

Main Scheme is: SIRIUS > MetFrag(PubChem) > MetFrag(KEGG) + Suspect List

1. if results has Sirius Structure annotation, and the explained inetnsity is >= 0.70, keep the annotation as is.
2. if annotation has PubChem, by default add SIRIUS And calculate the similarity between pubchem and SIrius results, if similar strcutres, then add Pubchme and sirius, if not then just keep sirius
3. apply same to metfrag results
4. if there is no annotation from SIRIUS, then use Metfrag results


### Spectral Database Dereplication Curation

Main Scheme is: SIRIUS > GNPS > Massbank > HMDB + Suspect List

1. Select only high scoring entries for further curation
2. Calculate the tanimoto similairity score between all entries from three spectral databases and also from suspect list. Add the database to the annotation column which passes the sructural similarity

### Combined Curation
if there are annoatations from both sources:
1. If there are three spectral database sources, check if the SMILES from spectral database and SIRIUS match, if yes, add Sirius to the annotation sources.
2. if there are two spectral databse sources, and HMDB is not in sources, still prioritize spectral database candidate if Sirius produces different structural candidate.
3. if there are two spectral databse sources, and HMDB is on of the sources, prioritize Sirius candidate if there is no strcutural similarity.
4. if there is only one spectral db source prioritize Sirius
5. Add SMILES from one of the sources mentioned in the Annotation_Source Column.


