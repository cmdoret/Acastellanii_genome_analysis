# Prepend strain name to scaffolds, genes and proteins to avoid name 
# collisions between strains.
rule rename_entries:
    input:
        genome = lambda w: samples.genome[w.strain],
        cds = lambda w: samples.cds[w.strain],
        proteome = lambda w: samples.proteome[w.strain],
        annotations = lambda w: samples.annotations[w.strain]
    output:
        genome = join(TMP, 'renamed', '{strain}_genome.fa'),
        cds = join(TMP, 'renamed', '{strain}_cds.fa'),
        proteome = join(TMP, 'renamed', '{strain}_proteome.fa'),
        annotations = join(TMP, 'renamed', '{strain}_annotations.gff')
    shell:
        """
        st={wildcards.strain}
        sed 's/^>\(.*\)/>'"$st"'_\\1/' {input.genome} > {output.genome}
        sed 's/^>\(.*\)/>'"$st"'_\\1/' {input.cds} > {output.cds}
        sed 's/^>\(.*\)/>'"$st"'_\\1/' {input.proteome} > {output.proteome}
        sed 's/^\([^#]\)/'"$st"'_\\1/' {input.annotations} |
            sed 's/ID=/ID='"$st"'_/' | 
            sed 's/Parent=/Parent='"$st"'_/' \
            > {output.annotations}
        """

# Just copying v1 assembly so that filename pattern matches other strains
rule rename_v1:
    input:
        genome = join(IN, 'genomes', 'NEFF_v1.fa'),
        annot = join(IN, 'annotations', 'NEFF_v1.43.gff')
    output:
        genome = join(TMP, 'renamed', 'NEFF_v1_genome.fa'),
        annot = join(TMP, 'renamed', 'NEFF_v1_annotations.gff')
    shell:
        """
        cp {input.genome} {output.genome}
        cp {input.annot} {output.annot}
        """


# 00 General annotations stats from amoeba GFF files
rule amoeba_annot_stats:
    input:
        v1 = join(TMP, 'renamed', 'NEFF_v1_annotations.gff'),
        neff = join(TMP, 'renamed', 'Neff_annotations.gff'),
        c3 = join(TMP, 'renamed', 'C3_annotations.gff')
    output:
        tbl = join(OUT, 'stats', 'annot_stats.tsv'),
        plt = join(OUT, 'plots', 'annot_stats.svg')
    conda: '../envs/r.yaml'
    shell: "Rscript scripts/01_annot_stats.R {input.v1} {input.neff} {input.c3} {output.tbl} {output.plt}"


# 00b Visualise assembly stats
rule plot_assembly:
    input: expand(join(TMP, 'renamed', '{strain}_genome.fa'), strain=samples.strain)
    output: join(OUT, 'plots', 'assembly_radars.svg')
    params:
        neff_v1 = join(IN, 'genomes', 'NEFF_v1.fa'),
        assembly_tbl = join(TMP, 'assembly_tbl.tsv')
    conda: '../envs/r.yaml'
    shell:
        """
        assembly-stats -t {input} {params.neff_v1} > {params.assembly_tbl}
        Rscript scripts/01_assembly_stats.R {params.assembly_tbl} {output}
        """

rule quast_report:
  input: expand(join(TMP, 'renamed', '{strain}_genome.fa'), strain=samples.strain)
  output: directory(join(OUT, 'plots', 'acastellanii_quast_report'))
  params:
    ref_fa = join(IN, 'genomes', 'NEFF_v1.fa'),
    ref_gff = join(IN, 'annotations', 'NEFF_v1.43.gff')
  conda: '../envs/quast.yaml'
  threads: NCPUS
  shell: 'quast -t {threads} -e -o {output} {input} {params.ref_fa}'

rule busco:
    input: join(TMP, 'renamed', '{assembly}_genome.fa')
    output: directory(join(OUT, 'busco', '{assembly}'))
    conda: '../envs/busco.yaml'
    threads: NCPUS
    shell:
        """
        busco -f --long --mode genome -c {threads} --lineage eukaryota -i {input} -o $(basename {output})
        mv "$(basename {output})" "$(dirname {output})/"
        """

# Clean busco summary to generate a usable tsv file
rule format_busco:
    input: join(OUT, 'busco', '{assembly}')
    output: join(OUT, 'busco', '{assembly}_summary.tsv')
    shell:
        """
        sed 's/^[\t ]*//' {input}/run_eukaryota_odb10/short_summary.txt \
          | grep -v '^[#\*]' \
          | sed '/^$/d' \
          | sed 's/[\t ]*$//' \
          | grep -v '\%' \
          > {output}
        """

# Generate BUSCO-style stacked barplots for all three assemblies
rule plot_busco:
    input:
        v1 = join(OUT, 'busco', 'NEFF_v1_summary.tsv'),
        neff = join(OUT, 'busco', 'Neff_summary.tsv'),
        c3 = join(OUT, 'busco', 'C3_summary.tsv')
    output: join(OUT, 'plots', 'busco_comparison.svg')
    conda: '../envs/viz.yaml'
    script: '../scripts/01_plot_busco.py'


# Compute chromosome sizes in new assembly
rule get_chrom_sizes:
    params:
        genome = lambda w: samples.genome[w.strain],
    output: join(TMP, '{strain}_chroms.sizes')
    shell:
        """
        awk -vOFS='\t' '
            /^>/ {{if (seqlen){{print seqlen}}; printf "%s\t",substr($0, 2) ;seqlen=0;next; }}
            {{ seqlen += length($0)}}END{{print seqlen}}
            ' {params.genome} > {output}
        """

rule agg_features:
    input:
        hgt = join(OUT, 'hgt', '{strain}_genes_hgt.bed'),
        vir = join(OUT, 'virus', '{strain}_summary.tsv')
    output: join(TMP, '{strain}_features.bed')
    params:
        rdna = lambda w: samples.rdna[w.strain]
    shell:
        """
        grep -v "^#" {params.rdna} | awk -vOFS='\t' '{{print $1,$4,$5,$9}}' > {output}
        awk -vOFS='\t' '$5 == 1 {{print $1,$2,$3,"HGT"}}' {input.hgt} >> {output}
        awk -vOFS='\t' '{{print $1,$2,$3,"virus"}}' {input.vir} >> {output}
        """

# Plot scaffolds
rule viz_scaffolds:
    log: join(OUT, 'log', '{strain}_viz_scaffolds.log')
    input:
        features = join(TMP, '{strain}_features.bed')
    output: join(OUT, 'plots', '{strain}_scaffolds.svg')
    params:
        genome = lambda w: samples.genome[w.strain],
        title =  lambda w: f"A. castellanii {w.strain}"
    conda: '../envs/viz.yaml'
    script: "../scripts/00_plot_karyo.py"