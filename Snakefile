# Genome analysis of Acanthamoeba castellanii Neff and "old" strains. Includes
# * Viral and bacterial gene content
# * Synteny analysis
# * Gene annotation stats
# * comparative analysis between both strains
# Input amobea genomes assembled using long reads, shotgun Illumina and Hi-C

# cmdoret, 20190410

shell.executable("/bin/bash")
from Bio import SeqIO
import re
from os.path import join
from src import fasta_utils as fu
from src import misc_utils as mu

# Strains Neff and "old" of Acanthamoeba castellanii
Acastellanii_strains = ["Ac_neff", "Ac_old"]
# Bacterial groups of interest
bacteria_groups = "Legionella,Rickettsia,Parachlamydia," + \
                  "Protochlamydia,Amoebophilus,Procabacter"
email = 'cmatthey@pasteur.fr'
# email = input("Please enter your email to download refseq genomes.")

# ===========================
# Set file paths if needed
DATA = 'data/'
GENOMES = join(DATA, "input", 'genomes')
ANNOT = join(DATA, "input", 'annotations')
OUT = join(DATA, 'out')
# ===========================

rule all: 
    #input: expand("blast/{amoeba}.txt", amoeba=strains)
    input : expand(join(GENOMES, 'others/{group}_euk_assoc.fasta'), group=["bacteria", 'virus'])

# Download bacterial and viral genomes from interesting groups from Refseq.
rule download_refseq:
    input: join(GENOMES, 'others/{group}_accession.txt')
    output: join(GENOMES, 'others/{group}_euk_assoc.fasta')
    run:
        num_ids = sum(1 for line in open(input[0], 'r'))
        with open(output[0], 'w') as out_fa, open(input[0], 'r') as in_id:
            id_number = 0
            for seq_id in in_id:
                mu.progbar(id_number, num_ids, "Downloading genomes")
                genome = fu.fetch_refseq_genome(seq_id, email=email)
                SeqIO.write(genome, out_fa, 'fasta')
                id_number += 1
    

rule make_blast_db:
    input: "fa/amoeba/{amoeba}.fa"
    output: "fa/amoeba/{amoeba}.blastdb"
    shell: "makeblastdb -in {input} -input_type fasta -dbtype nucl -out {output} && sleep 5"

rule query_viruses:
    input:
      query = "fa/virus/{amoeba}_viruses.fa",
      db_file = "fa/amoeba/{amoeba}.blastdb.nhr"
    output: "blast/{amoeba}.txt"
    threads: 4
    params:
        db = "fa/amoeba/{amoeba}.blastdb",
        of="6",
        e_val="1e-3"
    shell: "blastn -num_threads {threads} -db {params.db} -query {input.query} -evalue {params.e_val} -outfmt {params.of} -out {output}"

rule hist:
    input: "lengths/{id}.data"
    output: "hists/{id}.png"
    shell: "scripts/hist.R {input} {output}"

rule mcscanx_virus:
    input:
    output:
