# Genome analysis of Acanthamoeba castellanii Neff and "old" strains. Includes
# * HGT analysis
# * Synteny analysis
# * Gene annotation stats
# * comparative analysis between both strains (TODO)
# Input amobea genomes assembled using long reads, shotgun Illumina and Hi-C

# cmdoret, 20190410

shell.executable("/bin/bash")
from collections import defaultdict
import re
import time
from os.path import join, basename
import numpy as np
import pandas as pd
# library


configfile: "config.yaml"
#validate(config, schema='schemas/config.schema.yaml')
NCPUS = config['n_cpus']
samples = pd.read_csv(config['samples'], sep='\t', dtype=str, comment='#').set_index(['strain'], drop=False)
#validate(samples, schema="schemas/samples.schema.yaml")

# Load the list of organisms to be used for comparisons and use clean names for files (no spaces, no capitals)
organisms = pd.read_csv(config['compare_species'], sep='\t', comment='#', dtype=str)
organisms['clean_name'] = organisms['name'].apply(lambda n: n.lower().replace(" ", "_"))


# email = input("Please enter your email to download refseq genomes.")

# ===========================
# Set file paths if needed
OUT = config['out_dir']
TMP = config['tmp_dir']
IN = config['in_dir']
SHARED = join(IN, 'shared_assets')
GENOMES = join(SHARED, 'genomes')
DB = join(IN, 'db')

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
        join(OUT, 'plots', 'assembly_radars.svg'),
        join(OUT, 'plots', 'gene_families_venn.svg'),
        join(OUT, 'plots', 'acastellanii_quast_report'),
        join(OUT, 'stats', 'annot_stats.tsv'),
        join(OUT, 'plots', 'busco_comparison.svg'),
        expand(join(OUT, 'plots', 'virus_{strain}.svg'), strain=samples.strain),
        expand(join(OUT, 'virus', 'spatial', '{strain}_regions_pileup.txt'), strain=samples.strain),
        expand(join(OUT, 'virus', 'spatial', '{strain}_borders.tsv'), strain=samples.strain),
        expand(join(OUT, 'hgt', '{strain}_windows_hgt.tsv'), strain=samples.strain),
        expand(join(OUT, 'plots', '{strain}_scaffolds.svg'), strain=samples.strain),
        join(OUT, 'plots', 'c3_vs_neff_gap_compressed_div.svg')
        #join(OUT, 'figures', 'hgt_similarity.svg'),
        #join(OUT, 'orthofinder_blast', 'similarity_profile_bact.svg'),
        #expand(join(OUT, 'go_enrich', '{amoeba}_enrich.txt'), amoeba="Neff"),





include: 'rules/00_downloaders.smk'
include: 'rules/01_annot_stats.smk'
include: 'rules/02_horizontal_gene_transfer.smk'
include: 'rules/03_rdna.smk'
include: 'rules/04_synteny.smk'
include: 'rules/05_viruses.smk'
include: 'rules/06_phylogeny.smk'


# 10 GO enrichment test for HGT candidates
rule GO_enrich:
    input:
        annot = lambda w: samples.annotations[f'{w.amoeba}'],
        candidates = join(OUT, '{amoeba}_HGT_candidates.txt')
    output: join(OUT, 'go_enrich', '{amoeba}_enrich.txt')
    shell: "Rscript scripts/go_enrich.R {input.annot} {input.candidates} {output}"
