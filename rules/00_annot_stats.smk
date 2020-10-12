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
    shell: "Rscript scripts/00_annot_stats.R {input.v1} {input.neff} {input.c3} {output.tbl} {output.plt}"


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
    run:
        from collections import OrderedDict
        def read_busco(tbl, name):
            """
            Read a busco summary table into a
            dataframe and create relevant columns
            """
            df = pd.read_csv(tbl, names=['number', 'type'], sep='\t', header=None)
            df['strain'] = name 
            # Single letter busco class as a column
            df['busco'] = df['type'].str.extract(r'.*\(([A-Z])\).*')
            df.loc[pd.isnull(df.busco), 'busco'] = 'T'
            # Transform number of buscos to percentages
            df['perc'] = 100 * np.round(
            df['number'] / df.number[df.busco == 'T'].values[0],
                decimals=2
            )
            return df
        # Concatenate busco tables from all strains and set relevant
        # variables as indices
        busco = (
            pd.concat([read_busco(tbl, name) for name, tbl in input.items()])
            .reset_index(drop=True)
            .sort_values('strain', ascending=False)
            .set_index(['busco', 'strain'])
        )
        mpl.use('Agg')
        # Keep track of the order in which strains will be plotted
        str_to_num = OrderedDict({'v1': 0, 'neff': 1, 'c3': 2})
        # Map each BUSCO class to a color
        col_dict = OrderedDict([
            ('M', '#F04441'),
            ('F', '#F0E441'),
            ('D', '#3592C6'),
            ('S', '#58B4E8'),
        ])
        # We plot one busco class per iteration
        for i, t in enumerate(col_dict.keys()):
            # Subset busco table for current strain
            busco_type = busco.loc[t]
            # Retrieve x index from strains (deterministic order)
            r = [str_to_num[s] for s in busco_type.index.values]
            # Stacked barplot, uses cumulative height of previous bars
            if i:
                plt.barh(
                    r, busco_type.perc.values, color=col_dict[t], left=cum_height
                )
                cum_height += busco_type.perc.values
            else:
                cum_height = busco_type.perc.values
                plt.barh(r, busco_type.perc.values, color=col_dict[t])

        def busco_string(strain):
            """Use the BUSCO summary table to construct the standard busco string"""
            mis = busco.loc[('M', strain), 'number']
            fra = busco.loc[('F', strain), 'number']
            com = busco.loc[('C', strain), 'number']
            sin = busco.loc[('S', strain), 'number']
            dup = busco.loc[('D', strain), 'number']
            tot = busco.loc[('T', strain), 'number']
            return f'{strain} - C:{com}[S:{sin}, D:{dup}], F:{fra}, M:{mis}, n={tot}'
        # Write the strain name and busco string besides each bar
        plt.yticks(r, [busco_string(s) for s in str_to_num.keys()])
        plt.xlabel('% BUSCOs')
        plt.ylabel('Strain')
        plt.savefig(str(output))


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