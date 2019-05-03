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
Acastellanii_strains = ["Neff", "Old"]
# Bacterial and viral groups of interest
taxo_groups = {
    'bacteria': [
        "Legionella", "Rickettsia", "Parachlamydia", "Protochlamydia", 
        "Amoebophilus", "Procabacter"
        ],
                  
    'virus': [
        "Pandoravirus", "Marseillevirus", "Megavirus", 
        "Acanthamoeba polyphaga mimivirus", "Human herpesvirus", "Megavirus", 
        "Aureococcus anophagefferens virus", "Pithovirus sibericum", 
        "Lausannevirus"
        ]
}

email = 'cmatthey@pasteur.fr'
# email = input("Please enter your email to download refseq genomes.")

# ===========================
# Set file paths if needed
DATA = 'data/'
IN = join(DATA, 'input')
OUT = join(DATA, 'out')
TMP = join(DATA, 'tmp')
GENOMES = join(IN, 'genomes')
ANNOT = join(IN, 'annotations')
DB = join(IN, 'db')

# ===========================

rule all:
    input: 
        expand(join(OUT, "blast", "{groups}.{amoeba}.txt"), 
                   groups=list(taxo_groups.keys()), amoeba=Acastellanii_strains)


# 01 Download bacterial and viral genomes from interesting species
rule fetch_genomes:
    input:  join(GENOMES, 'others', '{group}_accession.txt')
    output: join(GENOMES, 'others', '{group}_euk_assoc.fa')
    run:
        num_ids = sum(1 for line in open(input[0], 'r'))
        with open(output[0], 'w') as out_fa, open(input[0], 'r') as in_id:
            id_number = 0
            for seq_id in in_id:
                mu.progbar(id_number, num_ids, "Downloading genomes")
                genome = fu.fetch_fasta(seq_id, email=email)
                SeqIO.write(genome, out_fa, 'fasta')
                id_number += 1

# 02 Download all protein sequences from bacterial and viral groups of interest
rule fetch_proteins:
    output: expand(join(ANNOT, 'others', '{group}_proteins.fa'),group=taxo_groups.keys()) 
    params:
    run:
        num_ids = len(taxo_groups['{group}'])
        with open(output[0], 'w') as out_fa:
            id_number = 0
            for organism in taxo_groups['{group}']:
                mu.progbar(id_number, num_ids, "Downloading proteins")
                fu.name_to_proteins(organism, email=email)
                id_number += 1

# 03 Make a nucleotide blast database of amoeba
rule make_blastn_db:
    input: join(GENOMES, "amoeba", "{amoeba}.fa")
    output: touch(join(TMP, '{amoeba}.blastdb'))
    params:
        db = join(DB, 'blast', '{amoeba}.blastdb')
    shell: "makeblastdb -in {input} -input_type fasta -dbtype nucl -out {params.db} && sleep 5"

# 04 Find other genomes which have hits in the amoeba
rule query_genomes:
    input:
      query = join(GENOMES, 'others', '{group}_euk_assoc.fa'),
      db_touch = join(TMP, '{amoeba}.blastdb')
    output: join(OUT, "blast", "{group}.{amoeba}.txt")
    threads: 4
    params:
        db = join(DB, 'blast', '{amoeba}.blastdb'),
        of="6",
        e_val="1e-10"
    shell:
        """
        blastn -num_threads {threads} \
               -db {params.db} \
               -query {input.query} \
               -evalue {params.e_val} \
               -outfmt {params.of} \
               -out {output}
        """

rule hist:
    input: "lengths/{id}.data"
    output: "hists/{id}.png"
    shell: "scripts/hist.R {input} {output}"

rule mcscanx_virus:
    input:
    output:
