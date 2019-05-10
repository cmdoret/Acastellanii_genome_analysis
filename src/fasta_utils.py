# Misc functions to operate on fasta files in the frame of the A. castellanii
# genome analysis project
# cmdoret, 20190502
from Bio import SeqIO
from Bio import Entrez
from BCBio import GFF
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


def fetch_fasta(seq_id, db="nucleotide", email="someone@email.com"):
    """
    Downloads a genome corresponding to input sequence ID.
    
    Parameters
    ----------
    seq_id : str
        A refseq sequence accession ID (e.g. NC_19130).
    db : str
        A valid Entrez database name. Some possible values are: cdd, gap, dbvar,
        epigenomics, nucest, gene, genome, gds, geoprofiles, nucgss, homologene,
        mesh, nuccore, protein, snp, sra, taxonomy, unigene
    email : str
        User email address to download from refseq.

    Returns
    -------
    seq_record : Bio.Seq
        Seq object containing the query Fasta record.
    """
    Entrez.email = email
    try:
        with Entrez.efetch(db=db, rettype="fasta", retmode="text", id=seq_id) as handle:
            seq_record = SeqIO.read(handle, "fasta")
    except:
        print("No sequence found for %s" % seq_id)
        seq_record = None
    return seq_record


def name_to_proteins(name, db="protein", email="someone@email.com"):
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

    Returns
    -------
    seq_record : TextIOWrapper
        iterable object containing the query Fasta records.
    
    """
    Entrez.email = email
    try:
        # First search to see the number of hits (returns values for 10 by default)
        query = Entrez.read(Entrez.esearch(term=name + "[Orgn]", db=db, email=email))
        # Get N hits
        count = query["Count"]
        # Real search, specifying max number of values
        query = Entrez.read(Entrez.esearch(term=name + "[Orgn]", db=db, retmax=count))
        seqs = Entrez.efetch(id=query["IdList"], rettype="fasta", db=db)
    except:
        print("No sequence found for %s" % name)
        seqs = None
    return seqs


def retrieve_id_annot(id, out_gff, mode="w", email="someone@email.com"):
    """
    Queries genbank record for an input ID and retrieves the genome annotations
    in GFF format. Amino acid sequences are included in the GFF.
    
    Parameters
    ----------
    id : str
        Sequence accession ID to query via Entrez.
    out_gff : str
        Path to the output GFF file.
    mode : str
        Mode in which to open the output GFF file. Should be 'w' or 'a'.
    email : str
        Personal email to provide when querying Entrez.

    """
    handle = Entrez.efetch(id=id, db="nucleotide", email=email, rettype="gb")
    record = SeqIO.parse(handle, "genbank")
    with open(out_gff, mode) as gff_handle:
        GFF.write(record, gff_handle, include_fasta=True)


def gff_seq_extract(gff, fa):
    """
    Extracts sequence from the attributes of CDS in a GFF into a fasta file

    Parameters
    ----------
    gff_in : str
        Path to the input GFF file containing "translation" and "protein_id" attributes.
    fa_out : str
        Path to the fasta file where the protein sequences should be written.
    """
    with open(gff, "r") as gff_in, open(fa, "w") as fa_out:
        for line in gff_in:
            seq_ok, id_ok = False, False
            fields = line.split("\t")
            if fields[2] == "CDS" and not fields[0].startswith("#>"):
                desc = fields[-1].split(";")
                for attr in desc:
                    if re.search("protein_id=", attr):
                        prot_id = attr.split("=")[1]
                        id_ok = True
                    elif re.search("translation=", attr):
                        seq = attr.split("=")[1]
                        seq_ok = True
            if seq_ok and id_ok:
                fa_out.writelines([prot_id, seq])

