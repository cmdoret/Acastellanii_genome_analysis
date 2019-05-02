# Misc functions to operate on fasta files in the frame of the A. castellanii
# genome analysis project
# cmdoret, 20190502
from Bio import SeqIO
from Bio import Entrez
import re


def retrieve_refseq_ids(in_ids, db, out_fa):
    """
    Given a refseq db in fasta format and a list of incomplete query refseq ID,
    extract the queried genomes into a new fasta file by matching the IDs.

    Parameters
    ----------
    in_ids : str
        Path to a file containing one refseq ID per line.
    db : str
        Path to the refseq database in fasta format.
    out_fa : str
        Path to the output fasta file containing query genomes.
    """
    query_ids = open(in_ids).read().splitlines()
    found = []
    with open(out_fa, "w") as genomes:
        for query_rec in SeqIO.parse(db, "fasta"):
            if re.search("|".join(query_ids), query_rec.id):
                query_rec.id = re.search(r"[^\.]*", query_rec.id).group()
                found.append(query_rec.id)
                SeqIO.write(query_rec, genomes, "fasta")
    print("%d genomes found among the %d queries." % (len(found), len(query_ids)))


def fetch_refseq_genome(seq_id, email="someone@email.com"):
    """
    Downloads a genome corresponding to input sequence ID.
    
    Parameters
    ----------
    seq_id : str
        A refseq sequence accession ID (e.g. NC_19130).
    email : str
        User email address to download from refseq.

    Returns
    -------
    seq_record : Bio.Seq
        Seq object containing the query Fasta record.
    """
    Entrez.email = email
    try:
        with Entrez.efetch(
            db="nuccore", rettype="fasta", retmode="full", id=seq_id
        ) as handle:
            seq_record = SeqIO.read(handle, "fasta")
    except:
        print("No sequence found for %s" % seq_id)
        seq_record = None
    return seq_record

