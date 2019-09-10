
# Build pangenome and protein tree from amoeba
rule orthofinder:
    input: 
        of_dir = join(OUT, 'amoeba_proteomes'),
        ac_proteomes = samples.proteome
    output: directory(join(OUT, 'orthofinder'))
    threads: NCPUS
    shell:
        """
        cp {input.ac_proteomes} {input.of_dir}
        orthofinder -f {input.of_dir} -o {output} -S diamond -t {threads} 
        """

# Get presence / absence info into GLOOME compatible format
rule format_presence_mat:
    input: join(OUT, 'orthofinder')
    output:
        matrix = join(TMP, 'GLOOME', 'gene_matrix.fa'),
        tree = join(TMP, 'GLOOME', 'tree.newick')
    shell:
        """
        cat {input}/Results_Sep10/Species_Tree/SpeciesTree_rooted.txt > {output.tree}
        python scripts/orthocount_to_fasta.py \
            {input}/Results_Sep10/Orthogroups/Orthogroups.GeneCount.tsv \
            {output.matrix}
        """


# Infer genes gain / loss events
rule GLOOME:
    input:
        matrix = join(TMP, 'GLOOME', 'gene_matrix.fa'),
        tree = join(TMP, 'GLOOME', 'tree.newick')
    params:
        conf_file = join(TMP, 'GLOOME', 'gloome.conf')
    output: directory(join(OUT, 'GLOOME'))
    shell: 
        """
        echo "_seqFile {input.matrix}" > {params.conf_file}
        # Tree seems to have incompatible format with CLI GLOOME but works in web version...
        #echo "_treeFile {input.tree}" >> {params.conf_file}
        echo "_outDir {output}" >> {params.conf_file}
        ./bin/gainLoss.VR01.266 {params.conf_file}
        """

# Retrieve A. castellanii-specific genes and strain specific genes
rule acastellanii_specific:
    input: join(OUT, 'orthofinder')
    output: 
        c3 = join(OUT, 'specific_genes', 'C3_specific.txt'),
        neff = join(OUT, 'specific_genes', 'Neff_specific.txt'),
        ac = join(OUT, 'specific_genes', 'Ac_specific.txt'),
        venn = join(OUT, 'specific_genes', 'acastellanii.svg')
    run:
        ortho = pd.read_csv(join(input[0], "Results_Sep10", "Orthogroups", "Orthogroups.tsv"), sep='\t')
        c3_abs = ortho.C3_proteins.isnull().values
        neff_abs = ortho.Neff_proteins.isnull().values
        ac_abs = c3_abs & neff_abs
        # Remove species of interest
        spec_df = ortho.drop(columns=['Orthogroup', 'C3_proteins', 'Neff_proteins'])
        # Binarize variables: Is there any gene for species c in each orthogroup
        spec_df = spec_df.apply(lambda c: c.isnull(), axis=0)
        # Summarize results: Does any species have at least one gene in each orthogroup
        out_pres = spec_df.apply(np.any, axis=1)
        out_abs = ~ out_pres
        c3_ortho = ortho.C3_proteins[neff_abs & (~ c3_abs) & out_abs]
        c3_ortho.to_csv(output['c3'], header=False)
        neff_ortho = ortho.Neff_proteins[c3_abs & (~ neff_abs) & out_abs]
        neff_ortho.to_csv(output['neff'], header=False)
        ac_ortho = ortho.Orthogroup[(~ ac_abs) & out_abs]
        ac_ortho.to_csv(output['ac'], header=False)

        # Make Venn diagram. sets are A, B, AB, C, AC, BC, ABC
        fig = plt.figure()
        venn3(
            subsets = (
                len(c3_ortho),
                len(neff_ortho),
                len(ac_ortho),
                len(ortho.Orthogroup[out_pres & ac_abs]),
                len(ortho.Orthogroup[(~ c3_abs) & out_pres & neff_abs]),
                len(ortho.Orthogroup[(~ neff_abs) & out_pres & c3_abs]),
                len(ortho.Orthogroup[(~ c3_abs) & (~ neff_abs) & out_pres])
            ), set_labels = ('A. castellanii C3', 'A. castellanii Neff', 'Amoeba'))
        fig.savefig(output['venn'])
