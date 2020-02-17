# Prepend strain name to scaffolds, genes and proteins to avoid name 
# collisions between strains.
rule rename_entries:
    input:
        genome = lambda w: samples.genome[f'{w.strain}'],
        proteome = lambda w: samples.proteome[f'{w.strain}'],
        annotations = lambda w: samples.annotations[f'{w.strain}']
    output:
        genome = join(TMP, 'renamed', '{strain}_genome.fa'),
        proteome = join(TMP, 'renamed', '{strain}_proteome.fa'),
        annotations = join(TMP, 'renamed', '{strain}_annotations.gff')
    shell:
        """
        st={wildcards.strain}
        sed 's/^>\(.*\)/>'"$st"'_\\1/' {input.genome} > {output.genome}
        sed 's/^>\(.*\)/>'"$st"'_\\1/' {input.proteome} > {output.proteome}
        sed 's/^\([^#]\)/'"$st"'_\\1/' {input.annotations} |
            sed 's/ID=/ID='"$st"'_/' | 
            sed 's/Parent=/Parent='"$st"'_/' \
            > {output.annotations}
        """

# 00 General annotations stats from amoeba GFF files
rule amoeba_annot_stats:
    input: join(TMP, 'renamed', '{amoeba}_annotations.gff')
    output:
        tbl = join(OUT, 'stats', '{amoeba}_annot_stats.tsv'),
        plt = join(OUT, 'plots', '{amoeba}_annot_stats.svg')
    conda: '../envs/r.yaml'
    shell: "Rscript scripts/00_annot_stats.R {input} {output.tbl} {output.plt}"


rule Neff_v1_annot_stats:
    input: join(IN, 'annotations', 'NEFF_v1.43.gff')
    output:
        tbl = join(OUT, 'stats', 'NEFF_v1_annot_stats.tsv'),
        plt = join(OUT, 'plots', 'NEFF_v1_annot_stats.svg')
    conda: '../envs/r.yaml'
    shell: "Rscript scripts/00_annot_stats.R {input} {output.tbl} {output.plt}"


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
        Rscript scripts/00_assembly_stats.R {params.assembly_tbl} {output}
        """

rule quast_report:
  input: expand(join(TMP, 'renamed', '{strain}_genome.fa'), strain=samples.strain)
  output: directory(join(OUT, 'plots', 'acastellanii_quast_report'))
  params:
    ref_fa = join(IN, 'genomes', 'NEFF_v1.fa'),
    ref_gff = join(IN, 'annotations', 'NEFF_v1.43.gff')
  conda: '../envs/quast.yaml'
  threads: NCPUS
  shell: 'quast -t {threads} -e -g {params.ref_gff} -r {params.ref_fa} -o {output} {input}'
