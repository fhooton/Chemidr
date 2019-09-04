import pandas as pd
import numpy as np
import urllib.request as request
import math
import pickle


def clean_term(term, convert_letter = True, w_space = True, is_url=True):
    term = term.lower().strip()
    
    if convert_letter:
        term = greek_letter_converter(term)
    else:
        term = greek_letter_converter(term, convert_letter=False)
    
    if w_space:
        if is_url:
            term = term.replace(' ', '%20') # To replace ' ' in request
        else:
            pass
    else:
        term = term.replace(' ', '')
    return term


def retrieve_info(req):
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" + req + "/synonyms/XML"
    
    with request.urlopen(url) as response:
        xml = response.read()
    
    root = root = etree.fromstring(xml)
    compound_id = float(root.findall(".//{http://pubchem.ncbi.nlm.nih.gov/pug_rest}CID")[0].xpath('.//text()')[0])
    compound_name = str(root.findall(".//{http://pubchem.ncbi.nlm.nih.gov/pug_rest}Synonym")[0].xpath('.//text()')[0])
    
    return compound_id, compound_name


def get_compound_pubchem_info(synonym):
    
    try:
        req = clean_term(synonym, convert_letter=False)
        compound_id, compound_name = retrieve_info(req)
    except:
        req = clean_term(synonym, convert_letter=False, w_space=False)
        compound_id, compound_name = retrieve_info(req)
    
    
    return compound_id, compound_name.lower()


import math
import pickle

def check_key(val, input_dict):
    """
    Evaluates if a cell already contains a key. If not, and if there is a match, insert a key.
    """
    if val is np.nan:
        if val in input_dict:
            input_dict[val]
        else:
            return np.nan
    else:
        return val
    

# Uses the foodb synonym file to create a synonym code for a df, and potentially cross reference the synonyms
# to another df to return what synonyms are unique
def append_synonyms(df, chem_key, use_foodb = True, load_ids=True):
    
    if use_foodb:
        df[chem_key] = df[chem_key].str.strip().str.lower()
        
        fdb_compounds = pd.read_csv('compounds.csv', encoding='latin1')[['id', 'name']]
        fdb_compounds.name = fdb_compounds.name.str.strip().str.lower()
        fdb_compounds.rename(columns={'id' : 'chem_id_f'}, inplace=True)
        
        df = df.merge(fdb_compounds, how = 'left', left_on = chem_key, right_on = 'name')
        
        if load_ids:
            with open('id_dicts/fdb_synonyms.pkl', 'rb') as f:
                input_dict = pickle.load(f)
        else:
            # Need to ensure document is in folder
            compound_synonyms = pd.read_csv('compound_synonymssql.csv')

            # Only keep columns with synonym and synonym id
            syn_reduced = compound_synonyms[['source_id', 'synonym']]
            syn_reduced = syn_reduced.rename(index=str, columns={"source_id": "chem_id_f"})
            syn_reduced.synonym = syn_reduced.synonym.str.lower()
        
            input_dict = {row['synonym'] : row['chem_id_f'] for _, row in syn_reduced.iterrows()}
            
            with open('id_dicts/fdb_synonyms.pkl', 'wb') as f:
                pickle.dump(input_dict, f)

        df.chem_id_f = df.chem_id_f.apply(check_key, input_dict=input_dict)
        
        if load_ids:
            with open('id_dicts/fdb_source_strings.pkl', 'rb') as f:
                input_dict = pickle.load(f)
        else:
            # Dataframe with contents of foodb
            foodb = pd.read_csv('contentssql.csv')

            # Gets the subset of the database pertaining to garlic
            foodb = foodb[['source_id', 'orig_source_name']].drop_duplicates()

            # Transforms all the chemical names to lowercase for syncing
            foodb.orig_source_name = foodb.orig_source_name.str.strip().str.lower()

            # Creates a list of the unique chemicals in garlic from foodb
            input_dict = {row['orig_source_name'] : row['source_id'] for _, row in foodb.iterrows()}
            
            with open('id_dicts/fdb_source_strings.pkl', 'wb') as f:
                pickle.dump(input_dict, f)
        
        df.chem_id_f = df.chem_id_f.apply(check_key, input_dict=input_dict)
        
        if load_ids:
            with open('id_dicts/usda_ids.pkl', 'rb') as f:
                input_dict = pickle.load(f)
        else:
            usda = pd.read_csv('usda_raw_garlic.csv', encoding = 'latin1')
            usda.nut_desc = usda.nut_desc.str.strip().str.lower()
            input_dict = {row['nut_desc'] : row['chem_id'] for _, row in usda.iterrows()}
            
            with open('id_dicts/usda_ids.pkl', 'wb') as f:
                pickle.dump(input_dict, f)
        
        df.chem_id_f = df.chem_id_f.apply(check_key, input_dict=input_dict)
                
    else:
        start = time()

        for idx, row in df.iterrows():
            try:
                ID, name = get_compound_pubchem_info(row[chem_key])
                df.at[idx, 'chem_id_p'] = ID
                df.at[idx, 'pubchem_name'] = name
                
            except:
                #print(row['chemical'])
                pass
            
            if not int(idx) % 20:
                print(idx, 'chems searched in', (time() - start) / 60, "min")

        print("Pubchem ids added in", (time() - start) / 60, "min")
    
    return df

def id_loader(df, chem_key, load, file, fdb = True, pubchem = True):
    
    if not load:
        if pubchem:
            df = append_synonyms(df, chem_key, use_foodb=False)
        if fdb:
            df = append_synonyms(df, chem_key)
        
        total = len(df[chem_key].drop_duplicates())
        
        num_covered = len(df[df.chem_id_p.notnull()].chem_id_p.drop_duplicates())
        print('Pubchem unique compound coverage', num_covered / total, '%')
        
        num_covered = len(df[df.chem_id_f.notnull()].chem_id_f.drop_duplicates())
        print('FooDB unique compound coverage', num_covered / total, '%')
        
        # The maximum pubchem index to make sure id's don't overlap
        max_p_index = 134825000
        
        for idx, row in df.iterrows():
            if not math.isnan(row['chem_id_p']):
                df.at[idx, 'chem_id'] = row['chem_id_p']
            else:
                df.at[idx, 'chem_id'] = row['chem_id_f'] + max_p_index
        
        num_covered = len(df[df.chem_id.notnull()].chem_id.drop_duplicates())
        print('Total unique compound covereage', num_covered / total, '%')

        file = 'misc_save/' + file
        df.to_pickle(file)
    
    return df