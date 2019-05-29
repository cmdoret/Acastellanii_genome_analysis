# Functions to request data using the Interpro RESTful API
# cmdoret, 20190520
import os
import time
import requests
import json
from omadb import Client as Oma
from omadb.OMARestAPI import ClientException, ClientTimeout


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

    response = requests.get(base_url.format(ipr_id=domain, taxid=str(taxid)))
    try:
        contents = json.loads(response.content)
        count = contents["count"]
    except (json.JSONDecodeError, KeyError):
        count = 0
    return count


def get_oma_hog(seq):
    """
    Given an input protein sequence, retrieve the finest taxonomic node for the
    best matching hierarchical ortholog group using the OMA api.
    
    Parameters
    ----------
    seq : str
        The amino-acid sequence of the query protein.
    
    Returns
    -------
    list of str :
        The list of taxonomic levels, from finest to largest, at which the HOG
        of the best hit is defined.
    """
    o = Oma()
    fail_countdown = 10
    while fail_countdown > 0:
        try:
            prots = o.proteins.search(seq)
            # First target is the closest hit
            targets = prots.targets
            fail_countdown = -1
        except (ClientException, ClientTimeout):
            print("Request failed, %i tries left" % fail_countdown)
            fail_countdown -= 1

    # If request did not work after 10 tries
    if fail_countdown == 0:
        print("All 10 requests failed, giving up.")
        targets = []
    if len(targets) == 0:
        organisms = [None]
    else:
        best_target = targets[0]
        # HOG taxonomic levels of the target go from finest to largest
        organisms = best_target.hog_levels
    return organisms
