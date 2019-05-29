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

vir_df = pd.read_csv(join(IN, 'misc', 'virus_names.tsv'), sep='\t', header=None)
bact_df = pd.read_csv(join(IN, 'misc', 'bacteria_names.tsv'), sep='\t', header=None)

# ===========================

rule all:
    input: 
        expand(join(OUT, "blast", "{groups}.{amoeba}.txt"), 
                   groups=list(taxo_groups.keys()), amoeba=Acastellanii_strains),
        join(OUT, 'orthoMCL', 'amoeba_groups.txt'),
        expand(join(TMP, "{group}_genomes_filtered.fa"), group=taxo_groups.keys()),
        expand(join(OUT, 'plots', '{amoeba}_annot_stats.pdf'), amoeba=["Neff", "NEFF_v1.43"]),
        join(TMP, 'Neff_hog_taxon.txt'),
        expand(join(OUT, '{amoeba}_sighunt.bed'), amoeba=Acastellanii_strains)


# 00 General annotations stats from amoeba GFF files
rule amoeba_annot_stats:
    input: join(ANNOT, 'amoeba', '{amoeba}.gff')
    output: join(OUT, 'plots', '{amoeba}_annot_stats.pdf')
    shell: "Rscript scripts/annot_stats.R {input} {output}"

# 00b Scan amoeba genomes to identify regions of different 4-mer signatures
rule sighunt_scan:
    input: join(GENOMES, 'amoeba', '{amoeba}.fa')
    output: join(OUT, '{amoeba}_sighunt.bed')
    shell: "Rscript scripts/sighunt_analysis.R {input} {output}"


# 01 Download bacterial and viral genomes from interesting species
rule fetch_genomes:
    input:  join(IN, 'misc', '{group}_names.tsv')
    output: join(TMP, 'genomes', '{group}_euk_assoc.fa')
    run:
        num_ids = sum(1 for line in open(input[0], 'r'))
        with open(output[0], 'w') as out_fa, open(input[0], 'r') as in_id:
            id_number = 0
            for organism in in_id:
                seq_id = organism.split("\t")[0]
                mu.progbar(id_number, num_ids, "Downloading %s genomes      " % wildcards['group'])
                time.sleep(0.1) # Do not spam NCBI :)
                genome = fu.fetch_fasta(seq_id, email=email)
                try:
                    SeqIO.write(genome, out_fa, 'fasta')
                    id_number += 1
                except TypeError:
                    pass


# 02 Download all protein sequences from bacterial and viral groups of interest
rule fetch_annotations:
    input: join(IN, 'misc', '{group}_names.tsv'),
    output: join(TMP, 'annot', '{group}_annot.gff')
    run:
        num_ids = sum(1 for line in open(input[0], 'r'))
        with open(input[0], 'r') as in_id:
            id_number = 0
            for organism in in_id:
                seq_id = organism.split('\t')[0]
                mu.progbar(id_number, num_ids, "Downloading %s annotations     " % wildcards['group'])
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
      query = join(TMP, 'genomes', '{group}_euk_assoc.fa'),
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
        genomes = join(TMP, 'genomes', '{group}_euk_assoc.fa'),
        annot = join(TMP, "annot", "{group}_annot.gff"),
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
                if rec.id in ids:
                    SeqIO.write(rec, out_fa, 'fasta')

        with open(input['annot'], 'r') as in_gff, open(output[1], 'w') as out_gff:
            for record in in_gff:
                if record.split('\t')[0] in ids:
                    out_gff.write(record)

# 06 Extract protein sequence from bact / vir GFF
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

            
# 08 Adjust fasta headers of amoeba proteomes for orthoMCL
rule prepare_orthoMCL_proteins:
    input:
        join(ANNOT, 'amoeba', '{amoeba}_proteins.fa')
    output:
        join(OUT, 'orthoMCL', 'fasta', '{amoeba}.fasta')
    run:
        taxon  = abbr["{wildcards.amoeba}"]
        with open(output[0], 'w') as mcl_fa:
            # All A. castellanii Neff headers prefixed with Acn|
            for prot in SeqIO.parse(input[0], 'fasta'):
                prot.id = taxon + "|" + prot.id
                prot.description = ""
                SeqIO.write(prot, mcl_fa, 'fasta')

# 09: Generate bed files corresponding to the orthoMCL fasta to use them in MCScanX
rule gen_orthoMCL_bed:
    input:
        join(ANNOT, 'amoeba', '{amoeba}_annotations.txt')
    output:
        join(OUT, 'orthoMCL', 'bed', '{amoeba}.bed')
    run:
        taxon = abbr["{wildcards.amoeba}"]
        with open(input[0]) as in_gff, open(output[0], 'w') as mcl_bed:
            for prot in in_gff:
                fields = prot.split('\t')
                bed_fields = {
                    'chr': fields[2], 
                    'start': fields[3], 
                    'end': fields[4], 
                    'id': taxon + "|" + fields[0]
                }
                mdl_bed.writeline("{chr}\t{start}\t{end}\t{id}".format(**bed_fields))
                

# 10 Pull and setup orthoMCL docker container, then run the ortholog group analysis
rule orthoMCL:
    input: expand(join(OUT, 'orthoMCL', 'fasta', '{amoeba}.fasta'), amoeba=Acastellanii_strains)
    output: 
        groups = join(OUT, 'orthoMCL', 'amoeba_groups.txt'),
        fasta = join(OUT, 'orthoMCL', 'amoeba_proteins_orthoMCL.fasta')
    threads: 12
    params:
        config = join('scripts', 'orthoMCL.config'),
        
    run:
        """
        cat {input} > {output.fasta}
        outdir=$PWD/$(dirname {output.groups})
        cp {params.config} $outdir
        bash scripts/orthoMCL_setup.sh {params.fasta} {threads} $outdir 
        """


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

# For each candidate prokaryotic gene, retrieve the closest taxonomic group at
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


# Compares new Neff and C3 assemblies
rule comp_amoeba:
    input:
    output:
    shell: "quast "
