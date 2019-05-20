# Functions to request data using the Interpro RESTful API
# cmdoret, 20190520
import os
import requests
import json


def get_domain_organisms(domain, taxid=None):
    """
    Returns the number of species in target taxonomic group for a given interpro domain.
    
    Parameters
    ----------
    taxid : int
        Uniprot taxid of interest. If None, the total number of hits for the 
        input domain is returned, otherwise only the number of hits in the given
        taxid is returned.
    domain : str
        Interpro domain ID to query.

    Returns
    -------
    int:
        Number of species in taxonomic group that contain target domain.

    Examples
    --------
    # Get hits for p53 tumour suppressor family in human
    >>> get_domain_organisms(9606, 'ipr002117')
    3
    # Same for all eukaryotes
    >>> get_domain_organisms(2759, 'ipr002117)
    40
    """
    base_url = "https://www.ebi.ac.uk/interpro/beta/api/protein/reviewed/entry/interpro/{ipr_id}"

    # Add taxid filter if requested.
    if taxid is not None:
        base_url += "/taxonomy/uniprot/{taxid}"

    response = requests.get(url.format(ipr_id=domain, taxid=str(taxid)))
    contents = json.loads(response.content)
    return contents["count"]
