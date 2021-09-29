# Download proteome of input organism from NCBI
import time
import urllib
from Bio import Entrez


def safe_request(fun):
    """
    Wraps function requesting data to allow safe errors and retry.
    
    Parameters
    ----------
    fun : python function
        The python function that queries a server
    
    Returns
    -------
    wrapped_f : python function
        A wrapped version of the input function which will call itself recursively
        every 5 seconds if the server is overloaded and will return None if the
        query record does not exist
    """

    def wrapped_f(*args, **kwargs):

        try:
            a = fun(*args, **kwargs)
            return a
        except urllib.error.HTTPError as e:
            if e.code in (429, 502):
                time.sleep(5)
                print(
                    "Sending too many requests, sleeping 5sec and retrying..."
                )
                wrapped_f(*args, **kwargs)
            else:
                print("No sequence found for %s" % args[0])
                return None

    return wrapped_f


@safe_request
def name_to_proteins(
    name, db="protein", email="someone@email.com", filters=""
):
    """
    Given an organism name, fetch all available protein sequences.

    Parameters
    ----------
    name : str
        Name of an organism or group of organism (e.g. Legionella).
    db : str
        Entrez db name to use for the search.
    email : str
        User email for identification.
    filters : str
        Additional filters to use for the query. E.g. " AND refseq[filter]" to
        filter results from refseq

    Returns
    -------
    seq_record : TextIOWrapper
        iterable object containing the query Fasta records.
    
    """
    Entrez.email = email
    query_string = name + "[Orgn]"
    time.sleep(0.1)
    # First search to see the number of hits (returns values for 10 by default)
    query = Entrez.read(
        Entrez.esearch(term=query_string, db=db, email=email, retmax=10 ** 9)
    )
    seqs = ""
    s = 0
    chunk_size = 1000
    # Fetch proteins by batch of 100 to avoid maximum number of queries
    for e in range(chunk_size, len(query["IdList"]), chunk_size):
        print(f"fetching entries {s} to {e}")
        fetched = False
        while not fetched:
            try:
                seqs += Entrez.efetch(
                    id=query["IdList"][s:e], rettype="fasta", db=db
                ).read()
                fetched = True
            except:
                print("Fetching failed, trying again.")
                time.sleep(5)
        s = e
        time.sleep(0.1)
    e = len(query["IdList"])
    print(f"fetching entries {s} to {e}")
    seqs += Entrez.efetch(
        id=query["IdList"][s:e], rettype="fasta", db=db
    ).read()
    return seqs


org_df = snakemake.params["org"]
organism = org_df.loc[
    org_df["clean_name"] == f"{snakemake.wildcards.organism}", "name"
]
time.sleep(0.1)  # Do not spam NCBI :)
proteome = name_to_proteins(
    organism, email="demo@email.com", filters=" AND refseq[filter]"
)
try:
    print(f"Writing {proteome.count('>')} proteins for {organism}.")
    with open(snakemake.output[0], "w") as outf:
        outf.write(proteome)
except TypeError:
    print(f"No proteome found for {organism[0]}")
    pass
