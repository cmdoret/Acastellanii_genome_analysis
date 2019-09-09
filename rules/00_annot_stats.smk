# 00 General annotations stats from amoeba GFF files
rule amoeba_annot_stats:
    input: lambda w: samples.annotations[f'{w.amoeba}']
    output: join(OUT, 'plots', '{amoeba}_annot_stats.svg')
    shell: "Rscript scripts/annot_stats.R {input} {output}"

# 00b Visualise assembly stats
rule radar_plot_assembly:
    input: expand(join(GENOMES, 'amoeba', '{amoeba}.fa'), amoeba=samples.strain + ["NEFF_v1"])
    output: join(OUT, 'plots', 'assembly_radars.svg')
    params:
        assembly_tbl = join(TMP, 'assembly_tbl.tsv')
    shell:
        """
        assembly-stats -t {input} > {params.assembly_tbl}
        Rscript scripts/radar_assembly_stats.R {params.assembly_tbl} {output}
        """