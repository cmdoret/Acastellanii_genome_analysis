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
import time
from os.path import join
import numpy as np
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

email = input("enter email adress for genome queries: ")
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
                   groups=list(taxo_groups.keys()), amoeba=Acastellanii_strains),
        join(OUT, "MCScanX", "MCScanX_genomes.fa"), 
        expand(join(TMP, "{group}_genomes_filtered.fa"), group=taxo_groups.keys()),
        expand(join(OUT, 'plots', '{amoeba}_annot_stats.pdf'), amoeba=["Neff", "NEFF_v1.43"])

# 00 General annotations stats from GFF files
rule amoeba_annot_stats:
    input: join(ANNOT, 'amoeba', '{amoeba}.gff')
    output: join(OUT, 'plots', '{amoeba}_annot_stats.pdf')
    shell: "Rscript scripts/annot_stats.R {input} {output}"

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
                time.sleep(1) # Do not spam NCBI :)
                genome = fu.fetch_fasta(seq_id, email=email)
                try:
                    SeqIO.write(genome, out_fa, 'fasta')
                    id_number += 1
                except TypeError:
                    pass
                


# 02 Download all protein sequences from bacterial and viral groups of interest
rule fetch_annotations:
    input:  join(GENOMES, 'others', '{group}_accession.txt')
    output: join(ANNOT, 'others', '{group}_annot.gff')
    run:
        num_ids = sum(1 for line in open(input[0], 'r'))
        with open(input[0], 'r') as in_id:
            id_number = 0
            for seq_id in in_id:
                mu.progbar(id_number, num_ids, "Downloading annotations")
                time.sleep(1)
                fu.retrieve_id_annot(seq_id, output[0], mode='a', email=email)
                id_number += 1


# 03 Make a nucleotide blast database of amoeba
rule make_amoeba_blastn_db:
    input: join(GENOMES, "amoeba", "{amoeba}.fa")
    output: touch(join(TMP, '{amoeba}.blastdb'))
    params:
        db = join(DB, 'blast', '{amoeba}.blastdb')
    shell: "makeblastdb -in {input} -input_type fasta -dbtype nucl -out {params.db} && sleep 5"


# 04 Find viral/bacterial genomes which have hits in the amoeba
rule blast_genomes_vs_amoeba:
    input:
      query = join(GENOMES, 'others', '{group}_euk_assoc.fa'),
      db_touch = join(TMP, '{amoeba}.blastdb')
    output: join(OUT, "blast", "{group}.{amoeba}.txt")
    threads: 4
    params:
        db = join(DB, 'blast', '{amoeba}.blastdb'),
        of="6",
        e_val="1e-11"
    shell:
        """
        blastn -num_threads {threads} \
               -db {params.db} \
               -query {input.query} \
               -evalue {params.e_val} \
               -outfmt {params.of} \
               -out {output}
        """


# 05 Filter annotations and genomes from hosts associated with amoeba
rule filter_genomes:
    input:
        genomes = join(GENOMES, 'others', '{group}_euk_assoc.fa'),
        annot = join(ANNOT, "others", "{group}_annot.gff"),
        select_id = expand(join(OUT, 'blast', '{group}.{amoeba}.txt'), group=taxo_groups.keys(), amoeba=Acastellanii_strains)
    output: 
        join(TMP, "{group}_genomes_filtered.fa"),
        join(TMP, "{group}_annot_filtered.gff")
    run:
        # Get unique virus / bact IDs that blasted agains amoeba genomes
        load_ids = lambda x: np.loadtxt(x, usecols=(0,), delimiter='\t', dtype=str)
        ids = [load_ids(arr) for arr in input['select_id'] if re.search(wildcards["group"], arr)]
        ids = np.concatenate(ids)
        ids = set(ids)
        # Write fasta with only genomes of those organisms
        with open(output[0], 'w') as out_fa:
            for rec in SeqIO.parse(input['genomes'], 'fasta'):
                print(rec.id)
                if rec.id in ids:
                    SeqIO.write(rec, out_fa, 'fasta')

        with open(input['annot'], 'r') as in_gff, open(output[1], 'w') as out_gff:
            for record in in_gff:
                if record.split('\t')[0] in ids:
                    out_gff.write(record)

# 06 Extract protein sequence from bact /vir GFF
rule extract_prot_seq:
    input: join(TMP, "{group}_annot_filtered.gff")
    output: join(TMP, "{group}_prot_filtered.fa")
    run:
        fu.gff_seq_extract(input[0], output[0])


# 07 Make blastn db for filtered bact/virus genomes
rule make_other_blastn_db:
    input: join(TMP, "{group}_genomes_filtered.fa")
    output: touch(join(TMP, '{group}.blastdb'))
    params:
        db = join(DB, 'blast', '{group}.blastdb')
    shell: "makeblastdb -in {input} -input_type fasta -dbtype nucl -out {params.db} && sleep 5"


# 06 Combine filtered genomes from all amoeba-host group combo and rm duplicates entries
rule merge_genomes:
    input:
        expand(
            join(TMP, "{group}_genomes_filtered.fa"), 
                 group=taxo_groups.keys(), 
                 amoeba=Acastellanii_strains
               )
    output: join(OUT, "MCScanX", "MCScanX_genomes.fa")
    run:
        consumed_ids = []
        with open(output[0], 'w') as out_fa:
            for genome in input[:]:
                for rec in SeqIO.parse(genome, 'fasta'):
                    if rec.id not in consumed_ids:
                        consumed_ids.append(rec.id)
                        SeqIO.write(rec, out_fa, 'fasta')

# 06 Get collinearity blocks between amoeba and viruses or bacteria
rule mcscanx_virus:
    input:
        combined_genomes = join(OUT, "MCScanX", "MCScanX_genomes_{group}.fa"),
        combined_annot = join(OUT, "MSCanX", "MCScanX_annot_{group},gff")
    output: touch(join(OUT, "MCScanX", "MCScanX.done"))
    params:
        out_dir = join(OUT, "MCScanX")
    shell:
        """
        bash scripts/gen_mcscanx_input.sh -g {input.combined_annot} \
                                          -o {params.out_dir} \
                                          -r {input.combined_genomes}
        """


# Compares new Neff assembly with old one
rule comp_amoeba:
    input:
    output:
    shell: "quast "
