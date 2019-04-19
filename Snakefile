# Analysis of the virus content in strains of amoeba based on a predetermined list of viruses.
# cmdoret, 20190410
shell.executable("/bin/bash")
from Bio import SeqIO
import re

# Strains Neff and "old" of Acanthamoeba castellanii
strains = ["Ac_neff", "Ac_old"]

rule all: 
    input: expand("blast/{amoeba}.txt", amoeba=strains)

rule get_virus_seq:
    input: 
        refseq = '/data/resources/db/refseq/refseq_viral_database_20190116.fasta',
        vir_id = 'virus/accession.txt'
    output: 
        ac_vir = "fa/virus/{amoeba}_viruses.fa"
    run:
      viruses = open(input.vir_id).read().splitlines()
      found = []
      with open(output.ac_vir, 'w') as ac_viruses:
        for virus in SeqIO.parse(input.refseq, 'fasta'):
          if re.search("|".join(viruses), virus.id):
            virus.id = re.search(r'[^\.]*', virus.id).group()
            found.append(virus.id)
            SeqIO.write(virus, ac_viruses, 'fasta')
      print("%d viruses found among the %d queries." % (len(found), len(viruses)))


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

