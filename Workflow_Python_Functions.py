#!/usr/bin/env python
# coding: utf-8

# ### SIRIUS post processing

# In[6]:


def sirius_postProc2(input_table, input_dir):
    for m, row in input_table.iterrows():
        file1  = pd.read_csv(input_table['ResultFileNames'][m] + '/insilico/MS1DATAsirius.csv')
        file1  = pd.read_csv(input_table['ResultFileNames'][m] + '/insilico/MS1DATAsirius.csv')
        file1['MCSSstring'] = np.nan
        file1['MCSS_SMILES'] = np.nan
        file1['Top_can_SL'] = np.nan
        file1['tanimotoSLvsCAN'] = np.nan
        file1['SL_comp'] = np.nan
        for i, row in file1.iterrows():
            if not isNaN(file1['SMILESforMCSS'][i]):
                top_smiles = file1['SMILESforMCSS'][i].split("|")
                mols = []
                for j in top_smiles:
                    for k, row in sl.iterrows():
                        SSms = [Chem.MolFromSmiles(j), Chem.MolFromSmiles(sl['SMILES'][k])]
                        SSfps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=1024) for x in SSms]
                        SStn = DataStructs.FingerprintSimilarity(SSfps[0],SSfps[1])
                        if SStn >= 0.8:
                            file1['Top_can_SL'][i] = j
                            file1['tanimotoSLvsCAN'][i] = SStn
                            file1['SL_comp'][i] = sl['SMILES'][k]
                    sm = Chem.MolFromSmiles(j)
                    mols.append(sm)
                if len(mols) > 1:
                    res = rdFMCS.FindMCS(mols)
                    file1['MCSSstring'][i] = res.smartsString
                    file1['MCSS_SMILES'][i] = Chem.MolToSmiles(Chem.MolFromSmarts(res.smartsString))
        file1.to_csv(input_table['ResultFileNames'][m] + '/insilico/SiriusResults.csv')
        
    all_files = []
    for n, row in input_table.iterrows():
        all_files.append(input_table['ResultFileNames'][n] + '/insilico/SiriusResults.csv')
        
    li = []
    
    for filename in all_files:
        df = pd.read_csv(filename, index_col=None, header=0)
        df["ResultFileNames"] = filename
        li.append(df)

    frame = pd.concat(li, axis=0, ignore_index=True)
    frame.to_csv(input_dir, '/SIRIUS_combined.csv')


# In[28]:


def classification(frame):
    inchis = []
    for l, row in frame.iterrows():
        if frame["Result"][l] == "SIRIUS_FOR":
            sep = 'json/'
            strpd = frame["dir"][l].split(sep, 1)[0] +"json/canopus_summary.tsv"
            if os.path.isfile(strpd):
                canopus = pd.read_csv(strpd, sep='\t')
                if len(canopus) > 0:
                    frame.loc[l, 'most_specific_class'] = canopus["most specific class"][0]
                    frame.loc[l, 'level _5'] = canopus["level 5"][0]
                    frame.loc[l, 'subclass'] = canopus["subclass"][0]
                    frame.loc[l, 'class'] = canopus["class"][0]
                    frame.loc[l, 'superclass'] = canopus["superclass"][0]
                    frame.loc[l, 'all_classifications'] = canopus["all classifications"][0]
        elif frame["Result"][l] == "SIRIUS_STR":
            InChI = Chem.MolToInchi(Chem.MolFromSmiles(frame["SMILES"][l]))
            InChIKey = Chem.inchi.InchiToInchiKey(InChI)
            inchis.append({
                'index': l,
                'smiles':frame["SMILES"][l],
                'inchi': InChI,
                'inchikey': InChIKey
            })
    inchis = pd.DataFrame(inchis)
    
    ## Retrieve ClassyFire classifications ##
    
    # This first step is done using inchikey and interrogation of the gnps classified structures
    gnps_proxy = True 
    url = "http://classyfire.wishartlab.com"
    proxy_url =  "https://gnps-classyfire.ucsd.edu"
    chunk_size = 1000
    sleep_interval = 12
    
    all_inchi_keys = list(inchis['inchikey'].drop_duplicates())

    resolved_ik_number_list = [0, 0]
    total_inchikey_number = len(all_inchi_keys)

    while True:
    
        start_time = time.time()
    
        print('%s inchikey to resolve' % total_inchikey_number )
        get_classifications_cf_mod(all_inchi_keys, par_level = 6)
    
        cleanse('all_json.json', 'all_json.json')
    
        with open("all_json.json") as tweetfile:
            jsondic = json.loads(tweetfile.read())

        df = json_normalize(jsondic)
        df = df.drop_duplicates( 'inchikey' )
        resolved_ik_number = len( df.drop_duplicates('inchikey').inchikey )
        resolved_ik_number_list.append( resolved_ik_number )
        print('%s resolved inchikeys' % resolved_ik_number )
        print("done in --- %s seconds ---" % (time.time() - start_time))
    
        if resolved_ik_number_list[-1] < resolved_ik_number_list[-2] or resolved_ik_number_list[-1] == resolved_ik_number_list[-3]:
            break
        cleanse('all_json.json', 'all_json_cleaned.json')
        
        with open("all_json_cleaned.json") as tweetfile:
            jsondic = json.loads(tweetfile.read())
            
    flattened_classified_json = json_normalize(jsondic)
    flattened_df = flattened_classified_json.drop_duplicates('inchikey')
    flattened_df['inchikey'] = flattened_df['inchikey'].str.replace(r'InChIKey=', '')
    df_merged = pd.merge(inchis, flattened_df, left_on='inchikey', right_on='inchikey', how='left')
    
    return(df_merged)


# In[ ]:





# In[ ]:





# ### MetFrag Result Post Processing

# In[ ]:




