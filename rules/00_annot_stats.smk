# 00 General annotations stats from amoeba GFF files
rule amoeba_annot_stats:
    input: lambda w: samples.annotations[f'{w.amoeba}']
    output: join(OUT, 'plots', '{amoeba}_annot_stats.svg')
    shell: "Rscript scripts/00_annot_stats.R {input} {output}"

rule Neff_v1_annot_stats:
    input: join(IN, 'annotations', 'NEFF_v1.43.gff')
    output: join(OUT, 'plots', 'NEFF_v1_annot_stats.svg')
    shell: "Rscript scripts/00_annot_stats.R {input} {output}"

# 00b Visualise assembly stats
rule radar_plot_assembly:
    input: samples.genome
    output: join(OUT, 'plots', 'assembly_radars.svg')
    params:
        neff_v1 = join(IN, 'genomes', 'NEFF_v1.fa'),
        assembly_tbl = join(TMP, 'assembly_tbl.tsv')
    shell:
        """
        assembly-stats -t {input} {params.neff_v1} > {params.assembly_tbl}
        Rscript scripts/00_radar_assembly_stats.R {params.assembly_tbl} {output}
        """