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
from src import interpro_utils as iu

# Strains Neff and "old" of Acanthamoeba castellanii
Acastellanii_strains = ["Neff", "C3"]
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
        expand(join(OUT, 'plots', '{amoeba}_annot_stats.pdf'), amoeba=["Neff", "NEFF_v1.43"])


# 00 General annotations stats from amoeba GFF files
rule amoeba_annot_stats:
    input: join(ANNOT, 'amoeba', '{amoeba}.gff')
    output: join(OUT, 'plots', '{amoeba}_annot_stats.pdf')
    shell: "Rscript scripts/annot_stats.R {input} {output}"


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

            
# 08 Adjust fasta headers of amoeba proteomes and combine all into a single fasta
rule prepare_orthoMCL_proteins:
    input:
        join(ANNOT, 'amoeba', 'Neff_proteins.fa'),
    output:
        join(OUT, 'orthoMCL', 'fasta', 'amoeba_proteins_orthoMCL.fasta')
    run:
        with open(output[0], 'w') as mcl_fa:
            # All A. castellanii Neff headers prefixed with Acn|
            for neff_prot in SeqIO.parse(input[0], 'fasta'):
                neff_prot.id = "Acn|" + neff_prot.id
                neff_prot.description = ""
                SeqIO.write(neff_prot, mcl_fa, 'fasta')
       

# 09 Pull and setup orthoMCL docker container, then run the ortholog group analysis
rule orthoMCL:
    input: join(OUT, 'orthoMCL', 'fasta', 'amoeba_proteins_orthoMCL.fasta')
    output: join(OUT, 'orthoMCL', 'amoeba_groups.txt')
    threads: 12
    params:
        config = join('scripts', 'orthoMCL.config')
    shell:
        """
        outdir=$PWD/$(dirname {output})
        cp {params.config} $outdir
        bash scripts/orthoMCL_setup.sh {input} {threads} $outdir 
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
            for domain in domains.split(";"):
                hits += iu.get_domain_organisms(domain, taxid)
            return hits

        prot_tbl = pd.read_csv(input[0], sep='\t', header=0)
        # Count domain hits in eukaryotes, bacteria and viruses for each prot
        groups = ["euk", "bac", "vir"]
        for group in groups:
            prot_tbl["n_%s" % group] = prot_tbl.InterPro.apply(
                    lambda x: sum_domains_hits(x, params[group])
                    )
        n_cols = [x for x in prot_tbl.columns if x.startswith('n_')]
        prot_tbl["n_all"] = prot_tbl.loc[:, n_cols].sum(axis=1)
        
        # Compute proportion of hits in each group
        for group in groups:
            prot_tbl["prop_%s" % group] = prot_tbl["n_%s" % group] / prot_tbl.n_all
        
        # Output newly generated columns to file
        out_cols = n_cols + ["n_all"] + [x for x in prot_tbl.columns if x.startswith('prop_')]
        prot_tbl.loc[:, out_cols].to_csv(output[0], sep='\t')

### KEEP FOR LATER ###
# 10 Combine filtered genomes from all amoeba-host group combo and rm duplicates entries
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


# Compares new Neff and C3 assemblies
rule comp_amoeba:
    input:
    output:
    shell: "quast "
