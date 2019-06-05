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

# Strains Neff and "old" of Acanthamoeba castellanii
Acastellanii_strains = ["Neff", "C3"]
# OrthoMCL compatible taxon names for both strains 
abbr = {"Neff": "Acn", "C3": "Acc"}
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
CIRCOS = join(IN, 'misc', 'circos_conf')

vir_df = pd.read_csv(join(IN, 'misc', 'virus_names.tsv'), sep='\t', header=None)
bact_df = pd.read_csv(join(IN, 'misc', 'bacteria_names.tsv'), sep='\t', header=None)

# ===========================

rule all:
    input:
        expand(join(OUT, "blast", "{groups}.{amoeba}.txt"), 
                   groups=list(taxo_groups.keys()), amoeba=Acastellanii_strains),
        join(OUT, 'orthoMCL', 'amoeba_groups.txt'),
        #expand(join(TMP, "{group}_genomes_filtered.fa"), group=taxo_groups.keys()),
        expand(join(OUT, 'plots', '{amoeba}_annot_stats.pdf'), amoeba=["Neff", "NEFF_v1.43"]),
        join(TMP, 'Neff_hog_taxon.txt'),
        expand(join(OUT, '{amoeba}_sighunt.bed'), amoeba=Acastellanii_strains),
        join(OUT, "MCScanX", "MCScanX.done"),
        expand(join(OUT, 'plots', 'circos_{amoeba}.svg'), amoeba="Neff"),
        expand(join(OUT, 'go_enrich', '{amoeba}_enrich.txt'), amoeba="Neff")

include: 'workflows/downloaders.smk'
include: 'workflows/orthomcl.smk'

# 00 General annotations stats from amoeba GFF files
rule amoeba_annot_stats:
    input: join(ANNOT, 'amoeba', '{amoeba}.gff')
    output: join(OUT, 'plots', '{amoeba}_annot_stats.pdf')
    shell: "Rscript scripts/annot_stats.R {input} {output}"

# 01 Scan amoeba genomes to identify regions of different 4-mer signatures
rule sighunt_scan:
    input: join(GENOMES, 'amoeba', '{amoeba}.fa')
    output:
      sig = join(OUT, '{amoeba}_sighunt.bed'),
      agg = join(OUT, '{amoeba}_sighunt.bed.agg')
    shell: "Rscript scripts/sighunt_analysis.R {input} {output.sig}"


# 10: Filter amoeba proteins containing prokaryotic interpro domains
rule interpro_filter:
    input: join(ANNOT, 'amoeba', '{amoeba}_annotations.txt')
    output: join(TMP, '{amoeba}_domain_groups.txt')
    params:
        euk = 2759,
        bac = 2,
        vir = 10239
    run:

        
        def sum_domains_hits(domains, taxid=None):
            """
            Wraps interpro_utils.get_domain_organisms to count hits in input
            taxid over multiple interpro domains separated by a semicolon.
            
            Parameters
            ----------
            domain : str
                InterPro domain to query, or multiple domains separated by ";".
            taxid : int or str
                Taxonomic group in which to count hits for the input domains.
            
            Returns
            -------
            int :
                The cumulative number of hits for input domains in the given taxid
            """

            hits = 0
            try:
                for domain in domains.split(";"):
                    try:
                        # Check if domain has been queried previously
                        hits = sum_domains_hits.memo[(domain, taxid)]
                    except KeyError:
                        # First time seeing this domain: send query
                        result = ou.get_domain_organisms(domain, taxid)
                        # Add result to the known combinations
                        sum_domains_hits.memo[(domain, taxid)] = result
                        hits += result
            # No domain available
            except AttributeError:
                pass
            
            # Show progress    
            sum_domains_hits.n_done += 1
            mu.progbar(
                sum_domains_hits.n_done, 
                sum_domains_hits.tot_prot, 
                'Querying interpro domains ({0}/{1})'.format(
                    sum_domains_hits.n_done, 
                    sum_domains_hits.tot_prot))
            return hits
        
        prot_tbl = pd.read_csv(input[0], sep='\t', header=0)
        # Count domain hits in eukaryotes, bacteria and viruses for each prot
        groups = ["euk", "bac", "vir"]
        sum_domains_hits.tot_prot = prot_tbl.shape[0] * len(groups)
        # Track progress
        sum_domains_hits.n_done = 0
        for group in groups:
            # Used to store the hit counts for each taxon-domain combination, in order
            # to only compute each combination once
            sum_domains_hits.memo = {} # Flush dictionary for new group
            prot_tbl["n_%s" % group] = prot_tbl.InterPro.apply(
                    lambda x: sum_domains_hits(x, params[group])
                    )
        
        # Get total hit count for each protein
        n_cols = [x for x in prot_tbl.columns if x.startswith('n_')]
        prot_tbl["n_all"] = prot_tbl.loc[:, n_cols].sum(axis=1)
        
        # Compute proportion of hits in each group
        for group in groups:
            prot_tbl["prop_%s" % group] = round(prot_tbl["n_%s" % group] / prot_tbl.n_all, 4)
        
        # Output newly generated columns to file
        out_cols = ["GeneID", "Contig", "Start", "Stop"] + n_cols + ["n_all"] + \
                   [x for x in prot_tbl.columns if x.startswith('prop_')]
        prot_tbl.loc[:, out_cols].to_csv(output[0], sep='\t', index=False)


# 06 Get collinearity blocks between amoeba and viruses or bacteria
rule mcscanx_amoeba:
    input:
        prot = join(ANNOT, "amoeba", "Neff_proteins.fa"),
        annot = join(ANNOT, "amoeba", "Neff.gff")
    output: touch(join(OUT, "MCScanX", "MCScanX.done"))
    threads: 12
    params:
        out_dir = join(OUT, "MCScanX")
    shell:
        """
        bash scripts/gen_mcscanx_input.sh -g {input.annot} \
                                          -o {params.out_dir} \
                                          -f {input.prot} \
                                          -c {threads}
        MCScanX -s 3 {params.out_dir}/MCScanX_in
        """

# 07 For each candidate prokaryotic gene, retrieve the closest taxonomic group at
# which the HOG is defined.
rule closest_hog:
    input:
        prot = join(ANNOT, 'amoeba', '{amoeba}_proteins.fa'),
        domains = join(TMP, '{amoeba}_domain_groups.txt')
    output: join(TMP, '{amoeba}_hog_taxon.txt')
    run:
        dom = pd.read_csv(input['domains'], sep='\t')
        dom['finest_taxon'] = ""
        dom['largest_taxon'] = ""
        dom = dom.loc[dom.prop_bac > 0.8, :]
        prok_genes = dom.GeneID.tolist()
        done_prok, tot_prok = 0, len(prok_genes)
        for prot in SeqIO.parse(input['prot'], 'fasta'):
            gene = prot.id.split("-")[0]
            if gene in prok_genes:
                orgs = ou.get_oma_hog(str(prot.seq))
                try:
                    dom.loc[dom.GeneID == gene, "finest_taxon"] = orgs[0]
                    dom.loc[dom.GeneID == gene, 'largest_taxon'] = orgs[-1]
                except IndexError:
                    dom.loc[dom.GeneID == gene, 'finest_taxon'] = None
                    dom.loc[dom.GeneID == gene, 'largest_taxon'] = None
                done_prok += 1
                mu.progbar(done_prok, tot_prok, "Fetching HOG taxons")
        dom.to_csv(output[0], sep='\t', index=False)


rule call_hgt_candidates:
    input:
        sighunt = join(OUT, '{amoeba}_sighunt.bed'),
        interpro = join(TMP, '{amoeba}_domain_groups.txt')
    output: join(OUT, '{amoeba}_HGT_candidates.txt')
    shell: "Rscript scripts/call_hgt.R {input.sighunt} {input.interpro} {output}"



# 08 Compute proportion of prokaryotic genes along genomic bins for visualisation
rule bin_interpro_genes:
    input: join(TMP, '{amoeba}_domain_groups.txt')
    output: join(OUT, '{amoeba}_interpro_bac_hist.txt')
    params:
        bin_size = 100000,
        threshold = 0.85
    run:
        mu.bin_tsv_file(
            input[0],
            output[0],
            params["bin_size"], 
            chrom=1, 
            start=2, 
            end=3, 
            val=9, 
            fun=lambda x: len(x[(x>params['threshold']) & (~np.isnan(x))]) / 
                          max(1, len(x[~np.isnan(x)]))
        )

# 09 Generate circos plot and required inputs
rule circos:
    input:
        ref = join(GENOMES, 'amoeba', '{amoeba}.fa'),
        sighunt = join(OUT, '{amoeba}_sighunt.bed.agg'),
        interpro = join(OUT, '{amoeba}_interpro_bac_hist.txt'),
        candidates = join(OUT, '{amoeba}_HGT_candidates.txt'),
        mcscx_flag = join(OUT, 'MCScanX', 'MCScanX.done')
    output:
        karyotype = join(OUT, 'plots', 'circos_{amoeba}.svg')
    params:
        mcsx_prefix = join(OUT, 'MCScanX', 'MCScanX_in')
    shell:
        """
        bash scripts/gen_circos_files.sh {input.ref} {input.sighunt} \
                                         {input.interpro} {input.candidates} \
                                         {params.mcsx_prefix}
        circos -conf {CIRCOS}/circos.conf -outputfile {output}
        """

# 10 GO enrichment test for HGT candidates
rule GO_enrich:
    input:
        annot = join(ANNOT, 'amoeba', '{amoeba}_annotations.txt'),
        candidates = join(OUT, '{amoeba}_HGT_candidates.txt')
    output: join(OUT, 'go_enrich', '{amoeba}_enrich.txt')
    shell: "Rscript scripts/go_enrich.R {input.annot} {input.candidates} {output}"