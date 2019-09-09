
# Build pangenome and protein tree from amoeba
rule:
    input: 
        of_dir = join(OUT, 'amoeba_proteomes'),
        ac_proteomes = samples.proteome
    output: touch(join(OUT, 'orthofinder.flag'))
    threads: NCPUS
    shell:
        """
        cp {input.ac_proteomes} {input.of_dir}
        orthofinder -f {input.of_dir} -S diamond -t {threads} 
        """

