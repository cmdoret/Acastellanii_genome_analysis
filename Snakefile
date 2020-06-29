# Genome analysis of Acanthamoeba castellanii Neff and "old" strains. Includes
# * HGT analysis
# * Synteny analysis
# * Gene annotation stats
# * comparative analysis between both strains (TODO)
# Input amobea genomes assembled using long reads, shotgun Illumina and Hi-C

# cmdoret, 20190410

shell.executable("/bin/bash")
from Bio import SeqIO
import re
import time
from os.path import join
import numpy as np
import pandas as pd
from src import fasta_utils as fu
from src import misc_utils as mu
from src import orthology_utils as ou
# library
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn2
mpl.use('TkAgg')


# Strains Neff and "old" of Acanthamoeba castellanii
#Acastellanii_strains = ["Neff", "C3"]
# OrthoMCL compatible taxon names for both strains 
#abbr = {"Neff": "Acn", "C3": "Acc"}
# Bacterial and viral groups of interest
#taxo_groups = {
#    'bacteria': [
#        "Legionella", "Rickettsia", "Parachlamydia", "Protochlamydia", 
#        "Amoebophilus", "Procabacter"
#        ],
#                  
#    'virus': [
#        "Pandoravirus", "Marseillevirus", "Megavirus", 
#        "Acanthamoeba polyphaga mimivirus", "Human herpesvirus", "Megavirus", 
#        "Aureococcus anophagefferens virus", "Pithovirus sibericum", 
#        "Lausannevirus"
#        ]
#}

configfile: "config.yaml"
#validate(config, schema='schemas/config.schema.yaml')
NCPUS = config['n_cpus']
samples = pd.read_csv(config['samples'], sep='\t', dtype=str, comment='#').set_index(['strain'], drop=False)
#validate(samples, schema="schemas/samples.schema.yaml")

# Load the list of organisms to be used for comparisons and use clean names for files (no spaces, no capitals)
organisms = pd.read_csv(config['compare_species'], sep='\t', comment='#', dtype=str)
organisms['clean_name'] = organisms['name'].apply(lambda n: n.lower().replace(" ", "_"))


email = 'cmatthey@pasteur.fr'
# email = input("Please enter your email to download refseq genomes.")

# ===========================
# Set file paths if needed
OUT = config['out_dir']
TMP = config['tmp_dir']
IN = config['in_dir']
GENOMES = join(IN, 'genomes')
DB = join(IN, 'db')
CIRCOS = join(IN, 'misc', 'circos_conf')

## WILDCARD CONSTRAINTS
wildcard_constraints:
  strain="|".join(samples.index),
  organism = '|'.join(organisms['clean_name'])

# ===========================

rule all:
    input:
        join(OUT, 'GLOOME'),
        join(OUT, 'plots', 'annot_stats.svg'),
        expand(join(OUT, 'plots', 'rdna_mat_{strain}.svg'), strain=samples.strain),
        join(OUT, "MCScanX", "MCScanX.done"),
        join(OUT, 'plots', 'circos.svg'),
        #expand(join(OUT, 'go_enrich', '{amoeba}_enrich.txt'), amoeba="Neff"),
        join(OUT, 'plots', 'assembly_radars.svg'),
        join(OUT, 'plots', 'gene_families_venn.svg'),
        #join(OUT, 'orthofinder_blast', 'similarity_profile_bact.svg'),
        join(OUT, 'plots', 'acastellanii_quast_report'),
        join(OUT, 'stats', 'annot_stats.tsv'),
        join(OUT, 'plots', 'hgt_stats.svg'),
        #join(OUT, 'figures', 'hgt_similarity.svg'),
        join(OUT, 'plots', 'busco_comparison.svg'),
        join(OUT, 'go_enrich', 'hgt_go_enrich.tsv')




include: 'rules/00_annot_stats.smk'
include: 'rules/01_downloaders.smk'
include: 'rules/02_horizontal_gene_transfer.smk'
include: 'rules/03_rdna.smk'
include: 'rules/04_synteny.smk'


# 10 GO enrichment test for HGT candidates
rule GO_enrich:
    input:
        annot = lambda w: samples.annotations[f'{w.amoeba}'],
        candidates = join(OUT, '{amoeba}_HGT_candidates.txt')
    output: join(OUT, 'go_enrich', '{amoeba}_enrich.txt')
    shell: "Rscript scripts/go_enrich.R {input.annot} {input.candidates} {output}"
