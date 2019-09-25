
import numpy as np
import math
import urllib.request as request
import requests
import time
import json
from lxml import etree


def cids2inchikeys(cids, as_dict=False, use_prefix=False):
	"""
        Retrieves InChIKeys from PubChem using Pubchem CIDS

        Input
        ----------------------------------------------------------------
        cids : list
            list of pubchem cid's for InChIKeys (needs to be ints, but also included int typecast)
        as_dict : bool (default False)
            returns dictionary of info if true, list otherwise
       	use_prefix : bool (default False)
       		only return the prefix of inchikeys (before the first -), which contains the structural information
       		(to find out more see https://www.inchi-trust.org/technical-faq-2/)

        Returns
        ----------------------------------------------------------------
        inchikeys : dict or list
            dictionary with CID's as keys and InChIKeys as values if as_dict is True, otherwise list
            of InChIKeys to preserve order
    """
	cids = __divide_list__([str(int(i)) for i in cids])

	if as_dict: inchikeys = {}
	else: inchikeys = []

	# Loop over divisions of ids to avoid overloading query
	for ids in cids:

		# Create url for InChIKey query
		url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{','.join(ids)}/property/InChIKey/JSON"

		
		r = __safe_urlopen__(url)

		# option to return InChIKey's as list or as dict (dict has certainty in case some cids aren't
		# retrieved, list preserves order)
		if as_dict:
			new_dict = {
				p['CID'] : p['InChIKey'] for p in json.loads(r)['PropertyTable']['Properties']
			}
			inchikeys.update(new_dict)

		else:
			new_list = [
				p['InChIKey'] for p in json.loads(r)['PropertyTable']['Properties']
			]
			inchikeys += new_list

	if use_prefix:
		if isinstance(inchikeys, dict): inchikeys = {cid : ikey.split('-')[0] for cid, ikey in inchikeys.items()}
		else: inchikeys = [ikey.split('-')[0] for ikey in inchikeys]

	return inchikeys


# Divides list into even divisions with a maximum of 100 elements
def __divide_list__(ids):
	num_divisions = int(math.ceil(len(ids) / 100))

	split_ids = np.array_split(np.asarray(ids), num_divisions)
	split_ids = [np.ndarray.tolist(split_ids[i]) for i in range(len(split_ids))]

	return split_ids


def __safe_urlopen__(url):
    response = requests.get(url)

    if response.status_code == 200: # Successful
        return response.content

    elif response.status_code == 429: # Too many requests
        # print('Retrying...')
        time.sleep(.5)
        return __safe_urlopen__(url)

    elif response.status_code == 503: # PUGREST.ServerBusy
        # print('Retrying...')
        time.sleep(1)
        return __safe_urlopen__(url)

    elif response.status_code == 404: # PUGREST.NotFound (aka doesn't exist)
        return None


def mesh2pid(mesh):
	url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pcsubstance&term={mesh}&retmode=json'

	r = __safe_urlopen__(url)

	if r is not None:
		j = json.loads(r)

		if j['esearchresult']['count'] != 0:

			sid = j['esearchresult']['idlist'][0] # get first sid result

			url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/{sid}/xml'

			xml = __safe_urlopen__(url)

			try:
				root = etree.fromstring(xml)

			except:
				print(sid)
				print(type(xml))
				print(xml)
				return

			cids = root.findall(".//{http://www.ncbi.nlm.nih.gov}PC-CompoundType_id_cid")

			if len(cids) > 0:
				cid = cids[0].xpath('./text()')[0]
			else:
				cid = np.nan

			return {mesh : {'mesh' : mesh, 'sid' : sid, 'cid' : cid}}

		else:
			return {mesh : {'mesh' : mesh, 'sid' : np.nan, 'cid' : np.nan}}

	else:
		return {mesh : {'mesh' : mesh, 'sid' : np.nan, 'cid' : np.nan}}

	# url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pccompound&term={mesh}&retmode=json'

	# r = __safe_urlopen__(url)

	# if r is not None:
	# 	j = json.reads(r)

	# 	if j['esearchresult']['count'] != 0:

	# 		sid = j['esearchresult']['idlist'][0] # get first sid result

	# 		url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/{sid}/xml'

	# 		xml = __safe_urlopen__(url)

	# 		root = etree.from_string(xml)

	# 		cids = root.findall(".//{http://www.ncbi.nlm.nih.gov}PC-CompoundType_id_cid")

	# 		if len(cids) > 0:
	# 			cid = cids[0].xpath('./text()')[0]
	# 		else:
	# 			cid = np.nan

	# 		return {'mesh' : mesh, 'sid' : sid, 'cid' : cid}