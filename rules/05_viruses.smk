VIR_GENOMES = join(TMP, 'virus_genomes')

# Download viruses genomes from NCBI/Refseq
rule get_viruse:
    output: directory(join(VIR_GENOMES, '{virus}'))
    message:
        """
        Downloading assembly {params.assembly} for {wildcards.virus}.
        """
    params:
        assembly = lambda w: organisms.loc[organisms.clean_name == w.virus, 'assembly'].values[0],
        genomedir = VIR_GENOMES
    conda: "../envs/genomepy.yaml"
    threads: 4
    shell: "genomepy install -g '{params.genomedir}/{wildcards.virus}' {params.assembly}"

# Align virus to amoeba genome
rule map_virus_amoeba:
    input:
        virus = join(VIR_GENOMES, '{virus}'),
        amoeba = lambda w: samples.genome[w.strain]
    output: join(OUT, 'virus', '{strain}', '{virus}.paf')
    threads: 6
    conda: '../envs/genomepy.yaml'
    shell: "minimap2 -xasm20 -t {threads} {params.amoeba} {input.virus}/*/*fa > {output}"


# Generate a bed file for each strain. Each virus-aligned segment is an interval with 
# similarity as score and virus species as name
rule parse_virus_matches:
    input: expand( join(OUT, 'virus', '{{strain}}', '{virus}.paf'), virus=organisms.loc[organisms['type'] == 'virus', 'clean_name'])
    output: join(OUT, 'virus', '{strain}_summary.tsv')
    params:
        pafdir = lambda w: join(OUT, 'virus', w.strain)
    conda: '../envs/genomepy.yaml'
    shell:
        """
        fd ".*paf" {params.pafdir} -x awk -vOFS='\t' -vvir={{/.}} '{{print $6,$8,$9,vir,$10/$11}}' {{}} > {output}
        """

rule aggregate_virus_matches:
    input: join(OUT, 'virus', '{strain}_summary.tsv')
    output: join(OUT, 'plots', 'virus_{strain}.svg')
    run:
        import matplotlib.pyplot as plt
        import seaborn as sns
        import matplotlib as mpl
        mpl.use("Agg")
        vir = pd.read_csv(
            input[0], sep="\t", names=["chrom", "start", "end", "virus", "ident"]
            )
        vir["seqlen"] = vir.end - vir.start
        # Amount of each virus inserted
        sns.scatterplot(data=vir, x="seqlen", y="ident", hue="virus")
        plt.savefig(output[0])

# Merge overlapping and neighbouring viral segments into viral "regions"
rule merge_virus_segments:
    input: join(OUT, 'virus', '{strain}_summary.tsv')
    output: join(OUT, 'virus', '{strain}_regions.tsv')
    params:
        neigh_dist = 10000
    conda: '../envs/genomepy.yaml'
    shell:
        """
        sort -k1,1 -k2,2n {input} \
            | bedtools merge -o mean -c 5 -i - -d {params.neigh_dist} \
            > {output}
        """

# Visualise the 3D structure around viral regions
rule viral_regions_pileup:
    input: join(OUT, 'virus', '{strain}_regions.tsv')
    output:
        tbl = join(OUT, 'virus', 'spatial', '{strain}_regions_pileup.txt'),
        fig = join(OUT, 'virus', 'spatial', '{strain}_regions_pileup.svg')
    params:
        cool = lambda w: samples.cool[w.strain] + "::/resolutions/8000"
    conda: "../envs/hic.yaml"
    shell:
        """
        # Discard coordinates with chromosomes absent from the cool file
        # and get the bed file into 2d coordinates
        cut -f1-3 {input} \
            | grep -f <(cooler dump -t chroms {params.cool} | awk '{{print $1"\t"}}') \
            | awk -vOFS='\t' '{{print $0,$1,$2+10000,$3+10000}}' \
            | coolpup.py {params.cool} \
                   - \
                   --nshifts 10 \
                   --mindist 0 \
                   --minshift 100000 \
                   --outname {output.tbl} \
                   --log WARNING
        plotpup.py {output.tbl} \
                   --col_names {wildcards.strain} \
                   --enrichment 0 \
                   --vmin 1 \
                   --output {output.fig}                  
        """
rule viral_regions_borders:
    input: join(OUT, 'virus', '{strain}_regions.tsv')
    output: multiext(join(OUT, 'virus', 'spatial', '{strain}_borders'), '.tsv', '.json', '.pdf')
    params:
        cool = lambda w: samples.cool[w.strain] + "::/resolutions/8000",
        base = join(OUT, 'virus', 'spatial', '{strain}_borders'),
        tpos = temp(join(TMP, '{strain}_virpos.tmp'))
    threads: 12
    conda: "../envs/hic.yaml"
    shell:
        """
        # Convert the bed file int 2d coordinates.
        cut -f1-3 {input} | awk -vOFS='\t' '{{print $0,$0}}' > {params.tpos}
        chromosight quantify --pattern=borders \
                             --threads {threads} \
                             --perc-zero=60 \
                             --perc-undetected=60 \
                             --win-size=31 \
                             {params.tpos} \
                             {params.cool} \
                             {params.base}  
        """