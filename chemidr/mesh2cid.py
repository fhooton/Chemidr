'''
Author: Jiali Cheng

Purpose: Replace MeshID with CID / SID

Use case: check the main function
'''

PREFIX = '/scratch/cheng.jial/foodome'
import pandas as pd
import pickle
import requests
from lxml import etree
import requests
import mpi4py.MPI as MPI
import argparse

def mesh2cid(mesh):
    '''
    Look up the chemical CID/SID based on the Mesh ID
    Input
    ----------------------------------------------------------------
    mesh : str
        Mesh ID of the chemical

    Returns
    ----------------------------------------------------------------
    output : dict
        dictionary of Mesh ID, 
        tag (whether it's CID or SID) and 
        value (CID value or SID value)
    '''
    try:
        url = 'http://www.ncbi.nlm.nih.gov/pcsubstance?term=%22{}%22%5bSourceID%5d'.format(mesh)
        response = requests.get(url).text
        html = etree.HTML(response)
        element = html[0].find("meta[@property='og:url']")
        info = element.get('content').split('/')
        tag = info[-2]
        value = info[-1]
        return {'mesh': mesh, 'tag': tag, 'value': value}
    except:
        return "None"


if __name__ == "__main__":

    comm = MPI.COMM_WORLD
    comm_rank = comm.Get_rank()
    comm_size = comm.Get_size()
    
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', dest='PATH', type=str, help='path of the result')
    parser.set_defaults(flag=False)
    args = parser.parse_args()

    cd = pd.read_csv('%s/data/pubmed-chemical-disease.csv'%PREFIX)
    mesh_list = cd.ChemicalID.drop_duplicates().values

    size = len(mesh_list) / comm_size
    local_mesh_list = mesh_list[comm_rank*size : comm_rank*size + size]
    local_result = [mesh2cid(i) for i in local_mesh_list]
    