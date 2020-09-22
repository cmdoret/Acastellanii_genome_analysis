# Rules to fetch resources and data from online sources

# 00 Download all protein sequences from groups of interest
# Note it takes for threads, but only uses 1 in practice. This is to limit the
# number of instances running (e.g. 3 on 12 threads and avoid spamming ncbi
rule fetch_proteomes:
    output: temp(join(OUT, 'proteomes', '{organism}_raw.fa'))
    params:
        org = organisms
    threads: 4
    run:
        org_df = params['org']
        organism = org_df.loc[org_df['clean_name'] == f'{wildcards.organism}', 'name']
        time.sleep(0.1) # Do not spam NCBI :)
        proteome = fu.name_to_proteins(organism, email=email, filters=" AND refseq[filter]")
        try:
            print(f"Writing {proteome.count('>')} proteins for {organism}.")
            with open(output[0], 'w') as outf:
                outf.write(proteome)
        except TypeError:
            print(f"No proteome found for {organism[0]}")
            pass

# Remove redundant (95% identical) proteins from each proteome
rule cdhit_proteomes:
  input: join(OUT, 'proteomes', '{organism}_raw.fa')
  output: join(OUT, 'filtered_proteomes', '{organism}.fa')
  params:
    sim = 0.95
  conda: '../envs/cdhit.yaml'
  shell: 'cd-hit -i {input} -o {output} -c {params.sim} -M 0'

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

