import pandas as pd
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs


# Convert SMILE's to chemical fingerprints
def get_fingerprint_string(SMILE):
    # Gets dictionary of subcompountnets and their counts
    sub_dict = AllChem.GetMorganFingerprint(Chem.MolFromSmiles(SMILE),1).GetNonzeroElements()

    fingerprint = []
    for key, value in sub_dict.items():
        fingerprint = fingerprint + [str(key)] * value

    return ' '.join(fingerprint)


# Organize fingerprints into correct molecular order
def get_ordered_fingerprint_string(SMILE):
    bi = {}
    AllChem.GetMorganFingerprint(Chem.MolFromSmiles(SMILE), radius=1, bitInfo=bi)

    mol = pd.DataFrame()
    for key, value in bi.items():
        for i in range(len(value)):
            sub = pd.Series()
            sub['val'] = str(key)

            sub['order'] = value[i][0]
            sub['radius'] = value[i][1]

            mol = mol.append(sub, ignore_index=True)

    mol = mol.sort_values(by=['order', 'radius']).reset_index(drop=True)
    #display(mol)
    
    mol_string = " ".join(mol.val.tolist())
        
    return mol_string.strip()


# def UNK_replacement(string):


def train_mol2vec(fingerprints, dim_embedding=100):
	fingerprints = [f.split() for f in chems.fingerprint.tolist()]

	model = Word2Vec(fingerprints, size=dim_embedding, sg=1, window=5, min_count=1, workers=4)
	model.train(fingerprints, total_examples=model.corpus_count, epochs=model.epochs)


def calc_fingerprint_vector(fingerprint, model=None):
    fingerprints = fingerprint.split()
    
    fingerprint_substructures = [model.wv[f] for f in fingerprints]
    
    fingerprint_vector = np.mean(fingerprint_substructures, axis=0)
    
    fingerprint_vector = ' '.join([str(i) for i in fingerprint_vector])
    
    return fingerprint_vector