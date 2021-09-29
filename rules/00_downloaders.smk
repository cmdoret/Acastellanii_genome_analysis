# Rules to fetch resources and data from online sources

# 00 Download all protein sequences from groups of interest
# Note it takes for threads, but only uses 1 in practice. This is to limit the
# number of instances running (e.g. 3 on 12 threads and avoid spamming ncbi
rule fetch_proteomes:
    output: join(OUT, 'proteomes', '{organism}_raw.fa')
    params:
        org = organisms
    threads: 4
    conda: '../envs/viz.yaml'
    script: '../scripts/00_fetch_proteomes.py'

# Remove redundant (99% identical) proteins from each proteome to discard duplicates
rule cdhit_proteomes:
  input: join(OUT, 'proteomes', '{organism}_raw.fa')
  output: join(OUT, 'filtered_proteomes', '{organism}.fa')
  params:
    sim = 0.99
  conda: '../envs/cdhit.yaml'
  shell: 'cd-hit -i {input} -o {output} -c {params.sim} -M 0'

# Download data assets from zenodo record
rule get_zenodo_assets:
  output:
    join(IN, 'genomes', 'NEFF_v1.fa'),
    join(IN, 'annotations', 'NEFF_v1.43.gff'),
    join(IN, 'cds', 'NEFF_v1.43.fa'),
    expand(join(IN, 'genomes', '{strain}_assembly.fa'), strain=['Neff', 'C3']),
    expand(join(IN, 'annotations', '{strain}_annotations.gff'), strain=['Neff', 'C3']),
    expand(join(IN, 'proteomes', '{strain}_proteins.fa'), strain=['Neff', 'C3']),
    expand(join(IN, 'annotations', 'rnammer', '{strain}.gff'), strain=['Neff', 'C3']),
    expand(join(IN, 'cds', '{strain}_cds.fa'), strain=['Neff', 'C3']),
    expand(join(IN, 'cool', '{strain}.mcool'), strain=['Neff', 'C3']),
    url_tbl = temp(join(TMP, 'zenodo_urls.tsv'))
  conda: '../envs/zenodo_get.yaml'
  params:
    in_dir = IN
  shell:
    """
    zenodo_get -d https://doi.org/10.5281/zenodo.5507417 -w {output.url_tbl}
    wget $(grep "shared_assets" {output.url_tbl}) -o shared_assets.tar.gz
    tar xzvf shared_assets.tar.gz --directory={params.in_dir}
    """

# Download bacterial and viral genomes from interesting species
#rule fetch_genomes:
#    input:  join(IN, 'misc', '{group}_names.tsv')
#    output: join(TMP, 'genomes', '{group}_euk_assoc.fa')
#    run:
#        num_ids = sum(1 for line in open(input[0], 'r'))
#        with open(output[0], 'w') as out_fa, open(input[0], 'r') as in_id:
#            id_number = 0
#            for organism in in_id:
#                seq_id = organism.split("\t")[0]
#                mu.progbar(id_number, num_ids, "Downloading %s genomes      " % wildcards['group'])
#                time.sleep(0.1) # Do not spam NCBI :)
#                genome = fu.fetch_fasta(seq_id, email=email)
#                try:
#                    SeqIO.write(genome, out_fa, 'fasta')
#                    id_number += 1
#                except TypeError:
#                    pass

#rule fetch_annotations:
#    input: join(IN, 'misc', '{group}_names.tsv'),
#    output: join(TMP, 'annot', '{group}_annot.gff')
#    run:
#        num_ids = sum(1 for line in open(input[0], 'r'))
#        with open(input[0], 'r') as in_id:
#            id_number = 0
#            for organism in in_id:
#                seq_id = organism.split('\t')[0]
#                mu.progbar(id_number, num_ids, "Downloading %s annotations     " % wildcards['group'])
#                time.sleep(1)
#                fu.retrieve_id_annot(seq_id, output[0], mode='a', email=email)
#                id_number += 1

