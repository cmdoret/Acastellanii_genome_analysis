# Rules to fetch resources and data from online sources

# 00 Download all protein sequences from groups of interest
rule fetch_proteomes:
    output: directory(join(OUT, 'proteomes'))
    params:
        org = organisms
    run:
        os.makedirs(output[0])
        n_orgs = params['org'].shape[0]
        for num_dl, org_row in enumerate(params['org'].iterrows()):
            organism = org_row[1][1]
            print(organism)
            mu.progbar(num_dl, n_orgs, "Downloading proteomes")
            time.sleep(0.1) # Do not spam NCBI :)
            proteome = fu.name_to_proteins(organism, email=email, filters=" AND refseq[filter]")
            fname = organism.lower().replace(" ", "_")
            try:
                print(f"Writing {proteome.count('>')} proteins for {organism}.")
                with open(join(output[0], fname + ".fa"), 'w') as outf:
                    outf.write(proteome)
            except TypeError:
                print(f"No proteome found for {organism[0]}")
                pass
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

