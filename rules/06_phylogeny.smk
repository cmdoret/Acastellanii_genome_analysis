
# Build pangenome and protein tree from amoeba and bacteria, also including HGT inferred
# in the first genome paper
rule orthofinder:
    input: 
        org_proteomes = expand(
            join(OUT, 'filtered_proteomes', '{organism}.fa'),
            organism=organisms.clean_name[organisms['type'] == 'amoeba']
        ),
        ac_proteomes = expand(
            join(TMP, 'renamed', '{strain}_proteome.fa'),
            strain=samples.strain
        )
        #hgt_neff_v1 = join(OUT, 'hgt', 'NEFF_v1_hgt_cds.fa')
    output: directory(join(OUT, 'hgt', 'orthofinder'))
    threads: NCPUS
    params:
      orthofinder_dir = join(TMP, 'orthofinder')
    conda: '../envs/orthofinder.yaml'
    shell:
        """
        ulimit -n 10000
        # Move all organisms proteomes to orthofinder workdir
        mkdir -p {params.orthofinder_dir}
        rm -f {params.orthofinder_dir}/*
        cp {input.org_proteomes} {params.orthofinder_dir}
        # Move Ac proteomes as well, but trim proteome from filename
        for strain in {input.ac_proteomes}; do
            fname=$(basename $strain)
            new_fname="${{fname/_proteome/}}"
            cp $strain "{params.orthofinder_dir}/$new_fname"
        done
        #cp "{{input.hgt_neff_v1}}" "{params.orthofinder_dir}/"
        orthofinder -f {params.orthofinder_dir} \
                    -o {output} \
                    -S diamond \
                    -t {threads} \
                    -n "amoeba"
        """

# Get presence / absence info into GLOOME compatible format
rule format_presence_mat:
    input: join(OUT, 'hgt', 'orthofinder')
    output:
        matrix = join(TMP, 'GLOOME', 'gene_matrix.fa'),
        tree = join(TMP, 'GLOOME', 'tree.newick')
    shell:
        """
        cat {input}/Results_amoeba/Species_Tree/SpeciesTree_rooted.txt > {output.tree}
        python scripts/06_orthocount_to_fasta.py \
            {input}/Results_amoeba/Orthogroups/Orthogroups.GeneCount.tsv \
            {output.matrix}
        """


# Infer genes gain / loss events
# Due to the lack of packaging options, this currently calls
# an embedded GLOOME binary which only works on linux
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

# Retrieve A. castellanii-specific + bact genes and strain specific genes
rule acastellanii_specific:
    input: join(OUT, 'hgt', 'orthofinder')
    output: 
        pres_mat = join(OUT, 'hgt', 'presence_matrix.tsv'),
        pres_compact = join(OUT, 'hgt', 'presence_compact.tsv')
    params:
        focus_st = samples.strain,
        amoeba = organisms.clean_name[organisms['type'] == 'amoeba']
        #bact = organisms.clean_name[organisms['type'] == 'bacteria']
    run:
        ortho_file = join(input[0], "Results_amoeba", "Orthogroups", "Orthogroups.tsv")
        ortho = pd.read_csv(ortho_file, sep='\t')
        # Get gene families absent in all of A. castellanii strains
        get_absent = lambda b: b.isnull().values
        st_abs = {strain: get_absent(ortho[strain]) for strain in params['focus_st']}
        ac_abs = np.bitwise_and.reduce([*st_abs.values()])
        
        # Convert count to presence/absence for all species
        pres_mat = ortho.drop(columns='Orthogroup').apply(lambda c: ~c.isnull(), axis=0)
        pres_mat.index = ortho.Orthogroup
        # Save presence matrix for further analysis (encore true/false as 0/1)
        pres_mat.astype(int).to_csv(output['pres_mat'], index=True, header=True, sep='\t')
        
        # Exclude strains of interest when counting
        amoeba_df = ortho.loc[:, list(params['amoeba'])]
        #bact_df = ortho.loc[:, list(params['bact'])]
        # Binarize variables: Is there any gene for species c in each orthogroup
        any_pres = lambda df: df.apply(lambda c: ~ c.isnull(), axis=0).apply(np.any, axis=1)
        #bact_pres = any_pres(bact_df)
        #bact_abs = ~bact_pres
        amoeba_pres = any_pres(amoeba_df)
        amoeba_abs = ~amoeba_pres
        # Record  genes found only in amoeba but not bacteria
        # Note: hgt events inferred in Clarke et al. also included for comparison)
        vdf = pd.DataFrame({
            'neff'  : ~st_abs['Neff'],
            #'NEFF_v1_hgt_cds'  : ~get_absent(ortho.NEFF_v1_hgt_cds),
            'c3'    : ~st_abs['C3'],
            'ac'    : ~ac_abs,
            #'bact'  : bact_pres,
            'amoeba': amoeba_pres,
        })
        vdf.index = ortho.Orthogroup
        
        vdf.astype(int).to_csv(output['pres_compact'], sep='\t', index=True, header=True)

rule plot_orthogroups_venn:
    input: 
        pres_compact = join(OUT, 'hgt', 'presence_compact.tsv')
    output:
        venn = join(OUT, 'plots', 'gene_families_venn.svg')
    run:
        mpl.use('Agg')
        vdf = pd.read_csv(input['pres_compact'], sep='\t').drop('Orthogroup', axis=1)
        vdf = vdf.astype(bool)
        #vdf = vdf.rename(columns={'NEFF_v1_hgt_cds': 'v1'})
        # Make Venn diagram. sets are A, B, AB, C, AC, BC, ABC
        fig, ax = plt.subplots(1, 2, figsize=(15, 12))
        venn3(
            subsets = (
                len(vdf.query('    c3 and not neff and not amoeba')), #A
                len(vdf.query('not c3 and     neff and not amoeba')), #B
                len(vdf.query('    c3 and     neff and not amoeba')), #AB
                len(vdf.query('not c3 and not neff and     amoeba')), #C
                len(vdf.query('    c3 and not neff and     amoeba')), #AC
                len(vdf.query('not c3 and     neff and     amoeba')), #BC
                len(vdf.query('    c3 and     neff and     amoeba')), #ABC
            ),
            set_labels = ('C3', 'Neff', 'Amoeba'),
            ax=ax[0]
        )
        venn2(
            subsets = (
                len(vdf.query('    c3 and not neff')),
                len(vdf.query('not c3 and     neff')),
                len(vdf.query('    c3 and     neff')),
            ),
        set_labels = ("C3", "Neff"),
        ax=ax[1]
        )
        fig.savefig(output['venn'])
