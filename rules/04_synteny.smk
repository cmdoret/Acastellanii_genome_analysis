
# 06 Get collinearity blocks between amoeba and viruses or bacteria
rule mcscanx_amoeba:
    input:
        prot = join(ANNOT, "amoeba", "Neff_proteins.fa"),
        annot = join(ANNOT, "amoeba", "Neff.gff")
    output: touch(join(OUT, "MCScanX", "MCScanX.done"))
    threads: 12
    params:
        out_dir = join(OUT, "MCScanX")
    shell:
        """
        bash scripts/gen_mcscanx_input.sh -g {input.annot} \
                                          -o {params.out_dir} \
                                          -f {input.prot} \
                                          -c {threads}
        MCScanX -s 3 {params.out_dir}/MCScanX_in
        """

# 09 Generate circos plot and required inputs
rule circos:
    input:
        ref = join(GENOMES, 'amoeba', '{amoeba}.fa'),
        sighunt = join(OUT, '{amoeba}_sighunt.bed.agg'),
        interpro = join(OUT, '{amoeba}_interpro_bac_hist.txt'),
        candidates = join(OUT, '{amoeba}_HGT_candidates.txt'),
        mcscx_flag = join(OUT, 'MCScanX', 'MCScanX.done')
    output:
        karyotype = join(OUT, 'plots', 'circos_{amoeba}.svg')
    params:
        mcsx_prefix = join(OUT, 'MCScanX', 'MCScanX_in')
    shell:
        """
        bash scripts/gen_circos_files.sh {input.ref} {input.sighunt} \
                                         {input.interpro} {input.candidates} \
                                         {params.mcsx_prefix}
        circos -conf {CIRCOS}/circos.conf -outputfile {output}
        """