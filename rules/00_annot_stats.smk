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
    input: join(IN, 'genomes', 'NEFF_v1.fa')
    output: join(TMP, 'renamed', 'NEFF_v1_genome.fa')
    shell: "cp {input} {output}"


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
        tbl = join(OUT, 'stats', 'NEFF_v1_oldannot_stats.tsv'),
        plt = join(OUT, 'plots', 'NEFF_v1_oldannot_stats.svg')
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
    run:
        from collections import OrderedDict
        def read_busco(tbl, name):
            df = pd.read_csv(tbl, names=['number', 'type'], sep='\t', header=None)
            df['strain'] = name 
            df['busco'] = df['type'].str.extract(r'.*\(([A-Z])\).*')
            df.loc[pd.isnull(df.busco), 'busco'] = 'T'
            df['perc'] = 100 * np.round(
	        df['number'] / df.number[df.busco == 'T'].values[0],
		decimals=2
	    )
            return df
        busco = (
            pd.concat([read_busco(tbl, name) for name, tbl in input.items()])
            .reset_index()
            .sort_values('strain', ascending=False)
       )
        mpl.use('Agg')
        r = range(len(input[:]))
        col_dict = OrderedDict([
            ('M', '#F04441'),
            ('F', '#F0E441'),
            ('D', '#3592C6'),
            ('S', '#58B4E8'),
        ])
	# TODO: Make sure strains y values are matched with their labels
        for i, t in enumerate(col_dict.keys()):
            if i:
                plt.barh(
                    r, busco.perc[busco.busco == t].values, color=col_dict[t], left=cum_height
                )
                cum_height += busco.perc[busco.busco == t].values
            else:
                cum_height = busco.perc[busco.busco == t].values
                plt.barh(r, busco.perc[busco.busco == t].values, color=col_dict[t])
        strains = np.unique(busco.strain.values)
        def busco_string(strain):
            str_df = busco.loc[busco.strain == strain, :]
            mis = str_df.number[str_df.busco == 'M'].values[0]
            fra = str_df.number[str_df.busco == 'F'].values[0]
            com = str_df.number[str_df.busco == 'C'].values[0]
            tot = str_df.number[str_df.busco == 'T'].values[0]
            return f'{strain} - C:{com}, F:{fra}, M:{mis}, n={tot}'
        
        plt.yticks(r, [busco_string(s) for s in strains])
        plt.xlabel('% BUSCOs')
        plt.ylabel('Strain')
        plt.savefig(str(output))

