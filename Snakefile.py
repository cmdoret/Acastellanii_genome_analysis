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
from matplotlib_venn import venn3
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

amoeba = pd.read_csv(config['amoeba_species'], sep='\t', comment='#', dtype=str)


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

vir_df = pd.read_csv(join(IN, 'misc', 'virus_names.tsv'), sep='\t', header=None)
bact_df = pd.read_csv(join(IN, 'misc', 'bacteria_names.tsv'), sep='\t', header=None)

# ===========================

rule all:
    input:
        join(OUT, 'GLOOME'),
        expand(join(OUT, 'plots', '{amoeba}_annot_stats.svg'), amoeba=["Neff", "C3", "NEFF_v1"]),
        expand(join(OUT, 'plots', 'rdna_mat_{amoeba}.svg'), amoeba=samples.strain),
        join(OUT, "MCScanX", "MCScanX.done"),
        join(OUT, 'plots', 'circos.svg'),
        #expand(join(OUT, 'go_enrich', '{amoeba}_enrich.txt'), amoeba="Neff"),
        join(OUT, 'plots', 'assembly_radars.svg'),
        join(OUT, 'specific_genes', 'acastellanii.svg'),
        join(OUT, 'orthofinder', 'blast', 'similarity_profile_bact.svg')


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
