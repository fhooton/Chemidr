"""
    Author: Forrest Hooton
    Purpose: Label given chemical strings with corresponding id's in popular databases (pubchem and foodb)

    Primary Functions:
        get_compound_pubchem_info: returns pubchem id and name given a string
        append_foodb_id: adds column with foodb id's given a dataframe with chemical strings
        append_pubchem_id: adds column with pubchem ids's given a dataframe with chemical strings
        id_searcher: adds columns with foodb and pubchem id's, along with a composite column

    Example Usage:
        test = pd.DataFrame({'chem' : ['allicin', 'diallyl disulphide']})
        id_searcher(test, 'chem')
"""

import pandas as pd
import numpy as np
from time import time
import urllib.request as request
from lxml import etree
import math
import pickle

import os
cwd = os.path.dirname(__file__) # get current location of script
package_path = 'Chemidr'.join(cwd.split('Chemidr')[:-1]) + 'Chemidr' # get path to head of Chemidr package


# Filepath wrapper to make data load-able from Chemidr
def __make_fp__(fp):
    return f'{package_path}/{fp}'

# Convert text to greet letter using *~greek letter name~*
def __greek_letter_converter__(chem, convert_letter = True):
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


def __clean_term__(term, convert_letter = True, w_space = True, is_url=True):
    """
        Prepares an input term to be quried in a url

        Input
        ----------------------------------------------------------------
        term : str
            term to clean
        convert_letter : bool (default True)
            Whether or not to convert letter to greek representation
        w_space : bool (default True)
            keep space in term, removes spaces if false
        is_url : bool (default True)
            Replaces spaces with %20 if true

        Returns
        ----------------------------------------------------------------
        term : str
            same term as input after cleaning process

    """
    term = term.lower().strip()
    
    term = __greek_letter_converter__(term, convert_letter=convert_letter)
    
    # Keeps spaces in string
    if w_space:
        if is_url:
            term = term.replace(' ', '%20') # To replace ' ' in request
        else:
            pass
    else:
        term = term.replace(' ', '')
    return term


def __retrieve_info__(req):
    """
        retrieves pubchem synonym information and extracts the compound id and primary compound name

        Input
        ----------------------------------------------------------------
        req : str
            chemical string to query in pubmed

        Returns
        ----------------------------------------------------------------
        compound_id : float
            pubchem compound id of req chemical
        compound_name : str
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
        req = __clean_term__(chem, convert_letter=False)
        compound_id, compound_name = __retrieve_info__(req)
    except:
        req = __clean_term__(chem, convert_letter=False, w_space=False)
        compound_id, compound_name = __retrieve_info__(req)
    
    return compound_id, compound_name.lower()


# Evaluates if a cell already contains a key. If not, and if there is a match, insert a key
def __check_key__(row, id_col, str_col, input_dict):
    if np.isnan(row[id_col]):
        if row[str_col] in input_dict:
            return input_dict[row[str_col]]
        else:
            return np.nan
    else:
        return row[id_col]
    

def append_foodb_id(df, chem_key, load_ids=True):
    """
        label chemicals with foodb id where there is a match to chemicals in foodb files under the data directory

        Input
        ----------------------------------------------------------------
        df : pd.DataFrame
            dataframe containing chemical names
        chem_key : str
            column name with the chemical strings
        load_ids : bool (default True)
            Whether or not to use pre-altered pubchem data. Want to use true after the program has run once to reduce runtime

        Returns
        ----------------------------------------------------------------
        df : pd.DataFrame
            input dataframe with a foodb_id column
    """
    df[chem_key] = df[chem_key].str.strip().str.lower()
    
    # Read in foodb primary compound names
    fdb_compounds = pd.read_csv(__make_fp__('data/compounds.csv'), encoding='latin1')[['id', 'name']]
    fdb_compounds.name = fdb_compounds.name.str.strip().str.lower()
    fdb_compounds.rename(columns={'id' : 'foodb_id'}, inplace=True)
    
    # Merge on matching names
    df = df.merge(fdb_compounds, how = 'left', left_on = chem_key, right_on = 'name')
    
    if load_ids:
        with open(__make_fp__('intermediate_save/fdb_synonyms.pkl'), 'rb') as f:
            input_dict = pickle.load(f)
    else:
        # Need to ensure document is in folder
        compound_synonyms = pd.read_csv(__make_fp__('data/compound_synonymssql.csv'))

        # Only keep columns with synonym and synonym id
        syn_reduced = compound_synonyms[['source_id', 'synonym']]
        syn_reduced = syn_reduced.rename(index=str, columns={"source_id": "foodb_id"})
        syn_reduced.synonym = syn_reduced.synonym.str.strip().str.lower()
    
        input_dict = {row['synonym'] : row['foodb_id'] for _, row in syn_reduced.iterrows()}
        
        with open(__make_fp__('intermediate_save/fdb_synonyms.pkl'), 'wb') as f:
            pickle.dump(input_dict, f)

    df.foodb_id = df.apply(__check_key__, id_col='foodb_id', str_col=chem_key, input_dict=input_dict, axis=1)
    
    if load_ids:
        with open(__make_fp__('intermediate_save/fdb_source_strings.pkl'), 'rb') as f:
            input_dict = pickle.load(f)
    else:
        # Dataframe with contents of foodb
        foodb = pd.read_csv(__make_fp__('data/contentssql.csv'))

        # Gets the subset of the database pertaining to garlic
        foodb = foodb[['source_id', 'orig_source_name']].drop_duplicates()

        # Transforms all the chemical names to lowercase for syncing
        foodb.orig_source_name = foodb.orig_source_name.str.strip().str.lower()

        # Creates a list of the unique chemicals in garlic from foodb
        input_dict = {row['orig_source_name'] : row['source_id'] for _, row in foodb.iterrows()}
        
        with open(__make_fp__('intermediate_save/fdb_source_strings.pkl'), 'wb') as f:
            pickle.dump(input_dict, f)
    
    df.foodb_id = df.apply(__check_key__, id_col='foodb_id', str_col=chem_key, input_dict=input_dict, axis=1)
    
    if load_ids:
        with open(__make_fp__('intermediate_save/usda_ids.pkl'), 'rb') as f:
            input_dict = pickle.load(f)
    else:
        usda = pd.read_csv(__make_fp__('data/usda_raw_garlic.csv'), encoding = 'latin1')
        usda.nut_desc = usda.nut_desc.str.strip().str.lower()
        input_dict = {row['nut_desc'] : row['chem_id'] for _, row in usda.iterrows()}
        
        with open(__make_fp__('intermediate_save/usda_ids.pkl'), 'wb') as f:
            pickle.dump(input_dict, f)
    
    df.foodb_id = df.apply(__check_key__, id_col='foodb_id', str_col=chem_key, input_dict=input_dict, axis=1)
                
    return df


def append_pubchem_id(df, chem_key):
    """
        label chemicals with pubchem id by searching the pubchem synonyms

        Input
        ----------------------------------------------------------------
        df : pd.DataFrame
            dataframe containing chemical names
        chem_key : str
            column name with the chemical strings

        Returns
        ----------------------------------------------------------------
        df : pd.DataFrame
            input dataframe with a pubchem_id column
    """
    start = time()

    i = 0
    for idx, row in df.iterrows():
        # get_compound_pubchem_info() creates an error if the chemical does not exist
        try:
            ID, name = get_compound_pubchem_info(row[chem_key])
            df.at[idx, 'pubchem_id'] = ID
            df.at[idx, 'pubchem_name'] = name
            
        except:
            # print(row[chem_key])
            pass
        
        if not i % 20:
            print(idx, 'chems searched in', (time() - start) / 60, "min")

        i += 1

    print("Pubchem ids added in", (time() - start) / 60, "min")

    return df


def id_searcher(df, chem_key, fdb = True, pubchem = True):
    """
        main function to assign chemical keys in pubchem and foodb

        Input
        ----------------------------------------------------------------
        df : pd.DataFrame
            dataframe containing chemical names
        chem_key : str
            column name with the chemical strings
        fdb : bool (default True)
            assign foodb ids
        pubchem : bool (default True)
            assign pubchem ids

        Returns
        ----------------------------------------------------------------
        df : pd.DataFrame
            input dataframe with a pubchem_id + foodb_id columns, and composite column chem_id
    """
    if pubchem:
        df = append_pubchem_id(df, chem_key)
    if fdb:
        df = append_foodb_id(df, chem_key)
    
    total = len(df[chem_key].drop_duplicates())
    
    # num_covered = len(df[df.pubchem_id.notnull()].pubchem_id.drop_duplicates())
    # print('Pubchem unique compound coverage', num_covered / total, '%')
    
    # num_covered = len(df[df.foodb_id.notnull()].foodb_id.drop_duplicates())
    # print('FooDB unique compound coverage', num_covered / total, '%')
    
    # Manually looked up the maximum pubchem index to make sure id's don't overlap
    max_p_index = 134825000
    
    # If there is no pubchem_id, make chem_id the foodb_id + the maximum pubchem id
    for idx, row in df.iterrows():
        if not math.isnan(row['pubchem_id']):
            df.at[idx, 'chem_id'] = row['pubchem_id']
        else:
            df.at[idx, 'chem_id'] = row['foodb_id'] + max_p_index
    
    # num_covered = len(df[df.chem_id.notnull()].chem_id.drop_duplicates())
    # print('Total unique compound covereage', num_covered / total, '%')

    # file = 'intermediate_save/' + file
    # df.to_pickle(file)
    
    return df