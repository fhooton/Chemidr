
import urllib.request as request
import json

def inchikeys_from_cids(cids):
	cids = [str(i) for i in cids]

	url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{','.join(cids)}/property/InChIKey/JSON"

	with request.urlopen(url) as response:
		r = response.read()

	inchikeys = [
		p['InChIKey'] for p in json.loads(r)['PropertyTable']['Properties']
	]

	return inchikeys
	