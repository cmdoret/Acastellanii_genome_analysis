# Rules involved in setting up and running orthoMCL


# 01 Adjust fasta headers of amoeba proteomes for orthoMCL
rule prepare_orthoMCL_proteins:
    input:
        join(ANNOT, 'amoeba', '{amoeba}_proteins.fa')
    output:
        join(OUT, 'orthoMCL', 'fasta', '{amoeba}.fasta')
    run:
        taxon = abbr["{wildcards.amoeba}"]
        with open(output[0], 'w') as mcl_fa:
            # All A. castellanii Neff headers prefixed with Acn|
            for prot in SeqIO.parse(input[0], 'fasta'):
                prot.id = taxon + "|" + prot.id
                prot.description = ""
                SeqIO.write(prot, mcl_fa, 'fasta')

# 02: Generate bed files corresponding to the orthoMCL fasta to use them in MCScanX
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
                

# 03 Pull and setup orthoMCL docker container, then run the ortholog group analysis
rule orthoMCL:
    input: expand(join(OUT, 'orthoMCL', 'fasta', '{amoeba}.fasta'), amoeba=Acastellanii_strains)
    output: 
        groups = join(OUT, 'orthoMCL', 'amoeba_groups.txt'),
        fasta = join(OUT, 'orthoMCL', 'amoeba_proteins_orthoMCL.fasta')
    threads: 12
    params:
        config = join(IN, 'misc', 'orthoMCL.config'),
        
    run:
        """
        cat {input} > {output.fasta}
        outdir=$PWD/$(dirname {output.groups})
        cp {params.config} $outdir
        bash scripts/orthoMCL_setup.sh {params.fasta} {threads} $outdir 
        """