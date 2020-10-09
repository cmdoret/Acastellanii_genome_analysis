
# Retrieve HGT sequences from Neff v1 publication and translate them
rule get_v1_hgt_fa:
    output: join(OUT, 'hgt', 'NEFF_v1_hgt_cds.fa')
    params:
        hgt_ids = join(IN, 'misc', 'NEFF_v1_HGT.tsv'),
        cds = join(IN, 'cds', 'NEFF_v1.43.fa')
    conda: '../envs/seqkit.yaml'
    shell: "seqkit grep -r -f {params.hgt_ids} {params.cds} | seqkit translate > {output}"

# Same for the gff3 file...
rule get_v1_hgt_gff:
    output: join(OUT, 'hgt', 'NEFF_v1_hgt.tsv')
    params:
        gff = join(IN, 'annotations', 'NEFF_v1.43.gff'),
        hgt = join(IN, 'misc', 'NEFF_v1_HGT.tsv')
    run:
        hgt_ids = {i.strip() for i in open(params['hgt'])}
        outf = open(output[0], 'w')
        with open(params['gff'], 'r') as gff:
            hgt = False
            n_exons = 0
            for line in gff:
                fields = line.split('\t')
                if re.match('^#', line):
                    continue
                elif fields[2] == 'gene':
                    # ID chrom start end type strand n_exon gene_len mrna_len attributes
                    if hgt:
                        outf.write(f"{curr_hgt.strip()}\t{n_exons}\n")
                    gene_id = re.search('gene:([^;]+);', fields[8]).group(1)
                    if gene_id in hgt_ids:
                        hgt = True 
                        n_exons = 0
                        curr_hgt = line
                    else:
                        hgt = False
                elif fields[2] == 'exon':
                    n_exons += 1


# Compute genome statistics in sliding windows on both genomes
rule genomic_windows:
    input: lambda w: samples.genome[w.strain]
    output: join(OUT, 'dnaglider', '{strain}_windows.tsv')
    params:
        win = 1000,
        step = 1000
    threads: 6
    shell:
        """
        dnaglider -threads {threads} \
                  -fields "GC,GCSKEW,ATSKEW,ENTRO,KMER" \
                  -kmers 2,3,4 \
                  -window {params.win} \
                  -stride {params.step} \
                  -fasta {input} \
                  -out {output}
        """

# Get corresponding v2 proteins to HGT detected in v1 (liftover)
rule blast_v1_hgt_vs_ac:
    input:
        v1_hgt = join(OUT, 'hgt', 'NEFF_v1_hgt_cds.fa')
    output:
        blast = join(OUT, 'blast', 'hgt_v1_vs_{strain}.blast')
    params:
        v2_prots = lambda w: samples.proteome[f"{w.strain}"],
        v2_db = lambda w: temp(join(TMP, 'blast', f'{w.strain}_db'))
    threads: 12
    shell:
        """
        makeblastdb -dbtype prot -in {params.v2_prots} -title neff -out {params.v2_db}
        blastp -evalue 1e-30 \
               -outfmt "6 qseqid sseqid pident evalue" \
               -db {params.v2_db} \
               -num_threads {threads} \
               -query {input.v1_hgt} > {output.blast}
        """

# Filter the best hit in v2 proteome for every HGT
rule filter_blast_hits:
    input: join(OUT, 'blast', 'hgt_v1_vs_{strain}.blast')
    output: join(OUT, 'blast', 'hgt_v1_vs_{strain}_filtered.blast')
    params:
        pident = 95
    run:
        df = pd.read_csv(
            input[0],
            sep='\t',
            header=None,
            names=['query', 'target', 'pident', 'evalue']
        )
        keep_idx = df.groupby('query')['pident'].transform(max) == df['pident']
        # For each query HGT, retain the match with highest sequence identity
        df = df[keep_idx]
        # Filter to keep only matches with > 95% identity
        df.loc[df.pident > params['pident']].to_csv(output[0], sep='\t', index=False)

# retrieve coordinates of lifted hgt and put a hgt (1/0) label on each gene annotation
rule get_hgt_coords_v2:
    input:
        blast = join(OUT, 'blast', 'hgt_v1_vs_{strain}_filtered.blast'),
        annot = join(OUT, 'stats', 'annot_stats.tsv')
    output: join(OUT, 'hgt', '{strain}_genes_hgt.bed')
    run:
        blast = pd.read_csv(input['blast'], sep='\t')
        annot = pd.read_csv( input['annot'], sep='\t')
        strain = wildcards['strain']
        annot = annot.loc[annot.ID.str.startswith(strain), :]
        annot.ID = annot.ID.str.replace(f"{strain}_", "")
        annot.chrom= annot.chrom.str.replace(f"{strain}_", "")
        blast['hgt'] = 1
        blast.target = blast.target.str.replace('-T1', '')
        annot = blast.merge(annot, how='right', left_on='target', right_on='ID')
        annot.hgt = annot.hgt.fillna(0).astype(int)
        annot = annot.loc[:, ['chrom', 'start', 'end', 'ID', 'hgt', 'n_exon']]
        annot.to_csv(output[0], sep='\t', header=None, index=False)


# Intersect genome windows and hgt to get average stats for each gene
rule windows_inter_hgt:
    input:
        win = join(OUT, 'dnaglider', '{strain}_windows.tsv'),
        hgt = join(OUT, 'hgt', '{strain}_genes_hgt.bed')
    output: join(OUT, 'hgt', '{strain}_windows_hgt.tsv')
    threads: 6
    shell:
        """
        echo -e "chrom\tstart\tend\tgeneID\tHGT\tNEXON\tGC\tGCSKEW\tATSKEW\tENTRO\t2MER\t3MER\t4MER" > {output}
        bedtools intersect -a {input.win} -b {input.hgt} -wb \
            | sort --parallel={threads} -k11,11 -k12,12n \
            | bedtools groupby -g 11,12,13,14,15,16 -c 4,5,6,7,8,9,10 -o mean \
            | sort --parallel={threads} -k1,1 -k2,2n >> {output}
        """
