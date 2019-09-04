"""
    Author: Forrest Hooton
    Purpose: Label given chemical strings with corresponding id's in popular databases
"""

import pandas as pd
import numpy as np
from time import time
import urllib.request as request
from lxml import etree
import math
import pickle

import os
<<<<<<< HEAD
cwd = os.path.dirname(__file__) # get current location of script
package_path = ''.join(cwd.split('Chemidr')[:-1]) + 'Chemidr'


# Filepath wrapper to make data load-able from Chemidr
def make_fp(fp):
    return f'{package_path}/{fp}'

# Convert text to greet letter using *~greek letter name~*
=======
cwd = os.path.dirname(__file__)
package_path = ''.join(cwd.split('Chemidr')[:-1]) + 'Chemidr'


# Make filepath
def make_fp(fp):
    return f'{package_path}/{fp}'

# Problem inputing letters into csv, so created system to convert them here
>>>>>>> fbc9ca5a8c2a70553e0c6f4b8bb5b75d06f46882
def greek_letter_converter(chem, convert_letter = True):
    if convert_letter:
        chem = chem.replace('*alpha*', 'α')
        chem = chem.replace('*beta*', 'β')
        chem = chem.replace('*gamma*', 'γ')
        chem = chem.replace('*rho*', 'ρ')
        chem = chem.replace('*delta*', 'δ')
    else:
        chem = chem.replace('*alpha*', 'alpha')
        chem = chem.replace('*beta*', 'beta')
        chem = chem.replace('*gamma*', 'gamma')
        chem = chem.replace('*rho*', 'rho')
        chem = chem.replace('*delta*', 'delta')
    return chem


def clean_term(term, convert_letter = True, w_space = True, is_url=True):
    """
        Prepares an input term to be quried in a url

        Input
        ----------------------------------------------------------------
        term : string
            term to clean
        convert_letter : bool (default True)
            Whether or not to convert letter to greek representation
        w_space : bool (default True)
            keep space in term, removes spaces if false
        is_url : bool (default True)
            Replaces spaces with %20 if true

        Returns
        ----------------------------------------------------------------
        term : string
            same term as input after cleaning process

    """
    term = term.lower().strip()
    
    term = greek_letter_converter(term, convert_letter=convert_letter)
    
    # Keeps spaces in string
    if w_space:
        if is_url:
            term = term.replace(' ', '%20') # To replace ' ' in request
        else:
            pass
    else:
        term = term.replace(' ', '')
    return term


def retrieve_info(req):
    """
        retrieves pubchem synonym information and extracts the compound id and primary compound name

        Input
        ----------------------------------------------------------------
        req : string
            chemical string to query in pubmed

        Returns
        ----------------------------------------------------------------
        compound_id : float
            pubchem compound id of req chemical
        compound_name : string
            compound name of primary compound for synonym entry
    """
    # Ex. https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/allicin/synonyms/XML
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{req}/synonyms/XML"
    
    with request.urlopen(url) as response:
        xml = response.read()
    
    # Extracts information from xml string
    root = root = etree.fromstring(xml)
    compound_id = float(root.findall(".//{http://pubchem.ncbi.nlm.nih.gov/pug_rest}CID")[0].xpath('.//text()')[0])
    compound_name = str(root.findall(".//{http://pubchem.ncbi.nlm.nih.gov/pug_rest}Synonym")[0].xpath('.//text()')[0])
    
    return compound_id, compound_name


# try-request wrapper to retrieve pubchem info
def get_compound_pubchem_info(chem):
    
    try:
        req = clean_term(chem, convert_letter=False)
        compound_id, compound_name = retrieve_info(req)
    except:
        req = clean_term(chem, convert_letter=False, w_space=False)
        compound_id, compound_name = retrieve_info(req)
    
    return compound_id, compound_name.lower()


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
def append_foodb_id(df, chem_key, load_ids=True):
    
    df[chem_key] = df[chem_key].str.strip().str.lower()
    
    fdb_compounds = pd.read_csv(make_fp('data/compounds.csv'), encoding='latin1')[['id', 'name']]
    fdb_compounds.name = fdb_compounds.name.str.strip().str.lower()
    fdb_compounds.rename(columns={'id' : 'foodb_id'}, inplace=True)
    
    df = df.merge(fdb_compounds, how = 'left', left_on = chem_key, right_on = 'name')
    
    if load_ids:
        with open(make_fp('intermediate_save/fdb_synonyms.pkl'), 'rb') as f:
            input_dict = pickle.load(f)
    else:
        # Need to ensure document is in folder
        compound_synonyms = pd.read_csv(make_fp('data/compound_synonymssql.csv'))

        # Only keep columns with synonym and synonym id
        syn_reduced = compound_synonyms[['source_id', 'synonym']]
        syn_reduced = syn_reduced.rename(index=str, columns={"source_id": "foodb_id"})
        syn_reduced.synonym = syn_reduced.synonym.str.lower()
    
        input_dict = {row['synonym'] : row['foodb_id'] for _, row in syn_reduced.iterrows()}
        
        with open(make_fp('intermediate_save/fdb_synonyms.pkl'), 'wb') as f:
            pickle.dump(input_dict, f)

    df.foodb_id = df.foodb_id.apply(check_key, input_dict=input_dict)
    
    if load_ids:
        with open(make_fp('intermediate_save/fdb_source_strings.pkl', 'rb')) as f:
            input_dict = pickle.load(f)
    else:
        # Dataframe with contents of foodb
        foodb = pd.read_csv(make_fp('data/contentssql.csv'))

        # Gets the subset of the database pertaining to garlic
        foodb = foodb[['source_id', 'orig_source_name']].drop_duplicates()

        # Transforms all the chemical names to lowercase for syncing
        foodb.orig_source_name = foodb.orig_source_name.str.strip().str.lower()

        # Creates a list of the unique chemicals in garlic from foodb
        input_dict = {row['orig_source_name'] : row['source_id'] for _, row in foodb.iterrows()}
        
        with open(make_fp('intermediate_save/fdb_source_strings.pkl'), 'wb') as f:
            pickle.dump(input_dict, f)
    
    df.foodb_id = df.foodb_id.apply(check_key, input_dict=input_dict)
    
    if load_ids:
        with open(make_fp('intermediate_save/usda_ids.pkl'), 'rb') as f:
            input_dict = pickle.load(f)
    else:
        usda = pd.read_csv(make_fp('data/usda_raw_garlic.csv'), encoding = 'latin1')
        usda.nut_desc = usda.nut_desc.str.strip().str.lower()
        input_dict = {row['nut_desc'] : row['chem_id'] for _, row in usda.iterrows()}
        
        with open(make_fp('intermediate_save/usda_ids.pkl'), 'wb') as f:
            pickle.dump(input_dict, f)
    
    df.foodb_id = df.foodb_id.apply(check_key, input_dict=input_dict)
                
    return df


def append_pubchem_id(df, chem_key):
    start = time()

    i = 0
    for idx, row in df.iterrows():
<<<<<<< HEAD
        try:
=======
        # try:
>>>>>>> fbc9ca5a8c2a70553e0c6f4b8bb5b75d06f46882
        ID, name = get_compound_pubchem_info(row[chem_key])
        df.at[idx, 'pubchem_id'] = ID
        df.at[idx, 'pubchem_name'] = name
            
<<<<<<< HEAD
        except:
            # print(row[chem_key])
            pass
=======
        # except:
            # print(row[chem_key])
            # pass
>>>>>>> fbc9ca5a8c2a70553e0c6f4b8bb5b75d06f46882
        
        if not i % 20:
            print(idx, 'chems searched in', (time() - start) / 60, "min")

        i += 1

    print("Pubchem ids added in", (time() - start) / 60, "min")

    return df


def id_searcher(df, chem_key, fdb = True, pubchem = True):
    
    if pubchem:
        df = append_pubchem_id(df, chem_key)
    if fdb:
        df = append_foodb_id(df, chem_key)
    
    total = len(df[chem_key].drop_duplicates())
    
<<<<<<< HEAD
    # num_covered = len(df[df.pubchem_id.notnull()].pubchem_id.drop_duplicates())
    # print('Pubchem unique compound coverage', num_covered / total, '%')
    
    # num_covered = len(df[df.foodb_id.notnull()].foodb_id.drop_duplicates())
    # print('FooDB unique compound coverage', num_covered / total, '%')
    
    # Manually looked up the maximum pubchem index to make sure id's don't overlap
=======
    num_covered = len(df[df.pubchem_id.notnull()].pubchem_id.drop_duplicates())
    print('Pubchem unique compound coverage', num_covered / total, '%')
    
    num_covered = len(df[df.foodb_id.notnull()].foodb_id.drop_duplicates())
    print('FooDB unique compound coverage', num_covered / total, '%')
    
    # The maximum pubchem index to make sure id's don't overlap
>>>>>>> fbc9ca5a8c2a70553e0c6f4b8bb5b75d06f46882
    max_p_index = 134825000
    
    for idx, row in df.iterrows():
        if not math.isnan(row['pubchem_id']):
            df.at[idx, 'chem_id'] = row['pubchem_id']
        else:
            df.at[idx, 'chem_id'] = row['foodb_id'] + max_p_index
    
<<<<<<< HEAD
    # num_covered = len(df[df.chem_id.notnull()].chem_id.drop_duplicates())
    # print('Total unique compound covereage', num_covered / total, '%')
=======
    num_covered = len(df[df.chem_id.notnull()].chem_id.drop_duplicates())
    print('Total unique compound covereage', num_covered / total, '%')
>>>>>>> fbc9ca5a8c2a70553e0c6f4b8bb5b75d06f46882

    # file = 'intermediate_save/' + file
    # df.to_pickle(file)
    
    return df