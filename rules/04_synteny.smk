
# 00 merge genomes and annotations
rule combine_strains:
    input:
        prot = samples.proteome,
        annot = samples.annotations,
        genome = samples.genome
    output:
        prot = join(TMP, 'merged', 'proteins.fa'),
        genome = join(TMP, 'merged', 'genome.fa'),
        annot = join(TMP, 'merged', 'annot.gff')
    shell:
        """
        echo -n "" > {output.prot}
        for i in {input.prot}; do
            # Prepend strain-specific id to each scaffold
            sp=$(basename $i)
            sp=${{sp%%_*}}
            sed 's/^>\(.*\)$/>'"$sp"'_\\1/' $i >> {output.prot}
        done

        echo -n "" > {output.genome}
        for i in {input.genome}; do
            # Prepend strain-specific id to each scaffold
            sp=$(basename $i)
            sp=${{sp%%_*}}
            sed 's/^>\(.*\)$/>'"$sp"'_\\1/' $i >> {output.genome}
        done

        echo -n "" > {output.annot}
        for i in {input.annot}; do
            sp=$(basename $i)
            sp=${{sp%%_*}}
            echo "----------"
            echo $i
            echo $sp
            echo "----------"
            sed 's/^/'"$sp"'_/' $i |
                sed 's/ID=/ID='"$sp"'_/' >> {output.annot}
        done
        """

# 06 Get collinearity blocks between amoeba and viruses or bacteria
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
rule circos:
    input:
        ref = join(TMP, 'merged', 'genome.fa'),
        #sighunt = join(OUT, '{amoeba}_sighunt.bed.agg'),
        #interpro = join(OUT, '{amoeba}_interpro_bac_hist.txt'),
        #candidates = join(OUT, '{amoeba}_HGT_candidates.txt'),
        mcscx_flag = join(OUT, 'MCScanX', 'MCScanX.done')
    output:
        karyotype = join(OUT, 'plots', 'circos.svg')
    params:
        mcsx_prefix = join(OUT, 'MCScanX', 'MCScanX_in'),
        circos_dir = join(IN, 'misc', 'circos_conf')
    conda: '../envs/circos.yaml'
    shell:
        """
        bash scripts/04_gen_circos_files.sh {input.ref} {params.mcsx_prefix} {params.circos_dir}
        bundlelinks -strict -min_bundle_size 1000 -max_gap 1000 \
                    -links {params.circos_dir}/mcsx.txt > {params.circos_dir}/bundles.txt
        circos -conf {CIRCOS}/circos.conf -outputfile {output}
        """