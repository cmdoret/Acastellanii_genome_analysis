# Analysing position of rDNA and relationship to 3D dna contact

# For now, files are generated externally by hicstuff (cool) and rnammer (rdna)
# TODO: Add file generation to analysis pipeline
rule plot_matrix_rDNA:
    input: 
        cool = lambda w: samples.cool[f'{w.amoeba}'],
        rdna = lambda w: samples.rdna[f'{w.amoeba}']
    output: join(OUT, 'plots', 'rdna_mat_{amoeba}.svg')
    params:
        res = 16000
    conda: '../envs/viz.yaml'
    shell:
        """
        python scripts/03_rdna_hic.py \
            {input.cool}::/resolutions/{params.res} \
            {input.rdna} \
            {output}
        """
