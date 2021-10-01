
# 00 merge genomes and annotations
rule combine_strains:
    input:
        prot = expand(join(TMP, 'renamed', '{strain}_proteome.fa'), strain=samples.strain),
        annot = expand(join(TMP, 'renamed', '{strain}_annotations.gff'), strain=samples.strain),
        genome = expand(join(TMP, 'renamed', '{strain}_genome.fa'), strain=samples.strain)
    output:
        prot = join(TMP, 'merged', 'proteins.fa'),
        genome = join(TMP, 'merged', 'genome.fa'),
        annot = join(TMP, 'merged', 'annot.gff')
    shell:
        """
        echo -n "" > {output.prot}
        for i in {input.prot}; do
            cat $i >> {output.prot}
        done

        echo -n "" > {output.genome}
        for i in {input.genome}; do
            cat $i >> {output.genome}
        done

        echo -n "" > {output.annot}
        for i in {input.annot}; do
            cat $i >> {output.annot}
        done
        """

# 06 Get collinearity blocks between amoeba and viruses or bacteria
# Requires: github.com/wyp1125/MCScanX
rule mcscanx_amoeba:
    input:
        prot = join(TMP, 'merged', 'proteins.fa'),
        annot = join(TMP, 'merged', 'annot.gff')
    output: touch(join(OUT, "MCScanX", "MCScanX.done"))
    threads: 12
    params:
        out_dir = join(OUT, "MCScanX")
    shell:
        """
        bash scripts/04_gen_mcscanx_input.sh -g {input.annot} \
                                             -o {params.out_dir} \
                                             -f {input.prot} \
                                             -c {threads}
        MCScanX -s 3 {params.out_dir}/MCScanX_in
        """

# 09 Generate circos plot and required inputs
# Note this requires circos-tools, which is distributed separately from circos
rule prep_circos:
    input:
        ref = join(TMP, 'merged', 'genome.fa'),
        mcscx_flag = join(OUT, 'MCScanX', 'MCScanX.done')
    output:
        mcsx = join(TMP, 'misc', 'circos_conf', 'mcsx.txt'),
        cfg = join(TMP, 'misc', 'circos_conf', 'circos.conf')
    params:
        mcsx_prefix = join(OUT, 'MCScanX', 'MCScanX_in'),
        dir = join(TMP, 'misc', 'circos_conf'),
        cfg = join(IN, 'misc', 'circos.conf')
    shell:
        """
        cp {params.cfg} {output.cfg}
        bash scripts/04_gen_circos_files.sh {input.ref} {params.mcsx_prefix} {params.dir}
        """

# Filter out links that involve chromosomes absent from karyotype (e.g. too short)
rule filter_links:
    input: join(TMP, 'misc', 'circos_conf', 'mcsx.txt')
    output: join(TMP, 'misc', 'circos_conf', 'mcsx_filtered.txt')
    params:
        circos_dir = join(TMP, 'misc', 'circos_conf')
    run:
        karyo = pd.read_csv(join(params['circos_dir'], 'karyotype.txt'), sep=' ', header=None)
        chroms = list(karyo.iloc[:, 2])
        links = pd.read_csv(input[0], sep=' ', header=None,
            names=['c1', 's1', 'e1', 'c2', 's2', 'e2', 'col']
        )
        filtered = links.loc[(np.isin(links.c1,chroms)) & (np.isin(links.c2, chroms)), :]
        filtered.to_csv(output[0], sep=' ',header=None, index=False)



rule circos:
    input:
        ref = join(TMP, 'merged', 'genome.fa'),
        mcsx = join(TMP, 'misc', 'circos_conf', 'mcsx_filtered.txt'),
        cfg = join(TMP, 'misc', 'circos_conf', 'circos.conf')
        #candidates = join(OUT, 'HGT_candidates.txt'),
    output:
        karyotype = join(OUT, 'plots', 'circos.svg')
    params:
        circos_dir = join(TMP, 'misc', 'circos_conf'),
    conda: '../envs/circos.yaml'
    shell:
        """
        bundlelinks -strict -min_bundle_size 1000 -max_gap 50000 \
                    -links {input.mcsx} | sed 's/lgrey=$//' > {params.circos_dir}/bundles.txt
        circos -conf {input.cfg} -outputfile {output}
        """

# Measure sequence divergence between both strains
# We use gap-excluded nucleotide divergence (i.e.
# Proportion of mismatch in aligned blocks)
rule get_divergence:
    input:
        c3 = join(TMP, 'renamed', 'C3_genome.fa'),
        neff = join(TMP, 'renamed', 'Neff_genome.fa'),
    output:
        paf = join(TMP, 'c3_vs_neff.paf'),
        div = join(OUT, 'div', 'c3_vs_neff_gap_excluded_div.bedgraph')
    conda: '../envs/genomepy.yaml'
    shell:
        """
        minimap2 -c {input.neff} {input.c3} -x map-ont \
        > {output.paf}
        python ./scripts/04_compute_seq_divergence.py \
            {output.paf} \
            {output.div}
        """

rule plot_divergence:
    input: join(OUT, 'div', 'c3_vs_neff_gap_excluded_div.bedgraph')
    output: join(OUT, 'plots', 'c3_vs_neff_gap_excluded_div.svg')
    conda: '../envs/viz.yaml'
    script: "../scripts/04_plot_seq_divergence.py"