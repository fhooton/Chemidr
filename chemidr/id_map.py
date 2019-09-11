
import urllib.request as request
import json

def inchikeys_from_cids(cids, as_dict=False):
	"""
        Retrieves InChIKeys from PubChem using Pubchem CIDS

        Input
        ----------------------------------------------------------------
        cids : list
            list of pubchem cid's for InChIKeys (needs to be ints)

        Returns
        ----------------------------------------------------------------
        inchikeys : dict or list
            dictionary with CID's as keys and InChIKeys as values if as_dict is True, otherwise list
            of InChIKeys to preserve order
    """
	cids = [str(i) for i in cids]

	# Create url for InChIKey query
	url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{','.join(cids)}/property/InChIKey/JSON"

	with request.urlopen(url) as response:
		r = response.read()

	# option to return InChIKey's as list or as dict (dict has certainty in case some cids aren't
	# retrieved, list preserves order)
	if as_dict:
		inchikeys = {
			p['CID'] : p['InChIKey'] for p in json.loads(r)['PropertyTable']['Properties']
		}

	else:
		inchikeys = [
			p['InChIKey'] for p in json.loads(r)['PropertyTable']['Properties']
		]

	return inchikeys
	