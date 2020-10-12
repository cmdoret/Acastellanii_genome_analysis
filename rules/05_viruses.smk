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
        virus = directory(join(VIR_GENOMES, '{virus}'))
    output: join(OUT, 'virus', '{strain}', '{virus}.paf')
    params:
        amoeba = lambda w: samples.genome[w.strain]
    threads: 6
    shell: "minimap2 -xasm20 -t {threads} {params.amoeba} {input.virus}/*/*fa > {output}"


# Generate a bed file for each strain. Each virus-aligned segment is an interval with 
# similarity as score and virus species as name
rule parse_virus_matches:
    input: expand( join(OUT, 'virus', '{{strain}}', '{virus}.paf'), virus=organisms.loc[organisms['type'] == 'virus', 'clean_name'])
    output: join(OUT, 'virus', '{strain}_summary.tsv')
    params:
        pafdir = lambda w: join(OUT, 'virus', w.strain)
    shell:
        """
        fd ".*paf" {params.pafdir} -x awk -vOFS='\t' -vvir={{/.}} '{{print $6,$8,$9,vir,$10/$11}}' {{}} > {output}
        """

rule aggregate_virus_matches:
    input:
        virus = join(OUT, 'virus', '{strain}_summary.tsv'),
        chroms = join(TMP, '{strain}_chroms.sizes'),
        hgt = join(OUT, 'hgt', '{strain}_genes_hgt.bed')
    output: join(OUT, 'plots', 'virus_{strain}.svg')
    shell: "python ./scripts/05_aggregate_viruses.py {input.virus} {input.chroms} {input.hgt} {output}"
