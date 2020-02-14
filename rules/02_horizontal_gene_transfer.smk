
# Build pangenome and protein tree from amoeba and bacteria
rule orthofinder:
    input: 
        org_proteomes = expand(
            join(OUT, 'proteomes', '{organism}.fa'),
            organism=organisms.clean_name
        ),
        ac_proteomes = expand(
            join(TMP, 'renamed', '{strain}_proteome.fa'),
            strain=samples.strain
        )
    output: directory(join(OUT, 'orthofinder'))
    threads: NCPUS
    params:
      orthofinder_dir = join(TMP, 'orthofinder')
    shell:
        """
        # Move all organisms proteomes to orthofinder workdir
        mkdir -p {params.orthofinder_dir}
        cp {input.org_proteomes} {params.orthofinder_dir}
        # Move Ac proteomes as well, but trim _proteome from filename
        for strain in {input.ac_proteomes}; do
            fname=$(basename $strain)
            new_fname="${{fname/_proteome/}}"
            cp $strain "{params.orthofinder_dir}/$new_fname"
        done
        orthofinder -f {params.orthofinder_dir} \
                    -o {output} \
                    -S diamond \
                    -t {threads} \
                    -n "amoeba"
        """

# Get presence / absence info into GLOOME compatible format
rule format_presence_mat:
    input: join(OUT, 'orthofinder')
    output:
        matrix = join(TMP, 'GLOOME', 'gene_matrix.fa'),
        tree = join(TMP, 'GLOOME', 'tree.newick')
    shell:
        """
        cat {input}/Results_amoeba/Species_Tree/SpeciesTree_rooted.txt > {output.tree}
        python scripts/02_orthocount_to_fasta.py \
            {input}/Results_amoeba/Orthogroups/Orthogroups.GeneCount.tsv \
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

# Retrieve A. castellanii-specific + bact genes and strain specific genes
rule acastellanii_specific:
    input: join(OUT, 'orthofinder')
    output: 
        venn = join(OUT, 'plots', 'gene_families_venn.svg'),
        pres_mat = join(OUT, 'specific_genes', 'presence_matrix.tsv'),
        pres_compact = join(OUT, 'specific_genes', 'presence_compact.tsv')
    params:
        focus_st = samples.strain,
        amoeba = organisms.clean_name[organisms['type'] == 'amoeba'],
        bact = organisms.clean_name[organisms['type'] == 'bacteria']
    run:
        mpl.use('TkAgg')
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
        bact_df = ortho.loc[:, list(params['bact'])]
        # Binarize variables: Is there any gene for species c in each orthogroup
        any_pres = lambda df: df.apply(lambda c: ~ c.isnull(), axis=0).apply(np.any, axis=1)
        bact_pres = any_pres(bact_df)
        bact_abs = ~bact_pres
        amoeba_pres = any_pres(amoeba_df)
        amoeba_abs = ~amoeba_pres
        # Record  genes found only in amoeba but not bacteria
        vdf = pd.DataFrame({
            'neff'  : ~st_abs['Neff'],
            'c3'    : ~st_abs['C3'],
            'ac'    : ~ac_abs,
            'bact'  : bact_pres,
            'amoeba': amoeba_pres,
        })
        vdf.index = ortho.Orthogroup
        vdf.astype(int).to_csv(output['pres_compact'], sep='\t', index=True, header=True)
        # Make Venn diagram. sets are A, B, AB, C, AC, BC, ABC
        fig, ax = plt.subplots(1, 2)
        venn3(
            subsets = (
                len(vdf.query('    ac and not bact and not amoeba')), #A
                len(vdf.query('not ac and     bact and not amoeba')), #B
                len(vdf.query('    ac and     bact and not amoeba')), #AB
                len(vdf.query('not ac and not bact and     amoeba')), #C
                len(vdf.query('    ac and not bact and     amoeba')), #AC
                len(vdf.query('not ac and     bact and     amoeba')), #BC
                len(vdf.query('    ac and     bact and     amoeba')), #ABC
            ),
            set_labels = ('A. castellanii', 'Bacteria', 'Amoeba'),
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

# Retrieve sequences from orthogroups belonging to specific groups
# for now, groups will be either hgt or background
rule get_spec_seq:
    input:
        pres = join(OUT, 'specific_genes', 'presence_compact.tsv'),
        orthofinder_dir = join(OUT, 'orthofinder')
    output:
        join(TMP, 'merged', '{group}_orthogroup_seq.fa')
    run:
        condition = {
            'hgt': 'ac and bact and not amoeba',
            'background': 'not bact and amoeba or ac'
        }
        # Retrieve HGT candidates (present in AC and bact,
        # but not other amoeba)
        pres = pd.read_csv(input['pres'], sep='\t')
        pres = pres.set_index(pres.Orthogroup)
        hgt_cand = pres.astype(bool).query(condition[f'{wildcards.group}'])
        out_fa = open(output[0], 'w')
        # Combine all sequences from candidate orthogroup into
        # a single fasta
        for og in hgt_cand.index:
            og_fa = join(
                input['orthofinder_dir'],
                'Results_amoeba',
                'Orthogroup_Sequences',
                f'{og}.fa'
            )
            og_handle = open(og_fa, 'r')
            for line in og_handle:
                out_fa.write(line)
            og_handle.close()
        out_fa.close()


# Filter Ac. genes from HGT candidiates or background sets
rule filter_ac_orthologues:
    input: join(TMP, 'merged', '{group}_orthogroup_seq.fa')
    output: join(TMP, 'merged', '{group}_orthogroup_seq_filtered.fa')
    params:
        id_filter = "^(C3|Neff)_FUN_"
    conda: '../envs/seqkit.yaml'
    shell: "seqkit grep -r -p '{params.id_filter}' {input} > {output}"


# Compute similarity profile between HGT candidates and bacteria 
# versus the rest of amoeba genes
rule bact_similarity:
    input:
        join(TMP, 'merged', '{group}_orthogroup_seq_filtered.fa')
    output:
        join(OUT, 'orthofinder_blast', '{group}_orthogroup_bact.blast')
    threads: NCPUS
    params:
        db = join(DB, 'blast', 'nr_v5', 'nr_v5')
        #db = join(DB, 'blast', 'nr_v5', 'nr.v5')
    singularity: 'docker://cmdoret/blast:2.9.0'
    shell:
        """
        blastp -num_threads {threads} \
               -max_target_seqs 5 \
               -taxids 2 \
               -db {params.db} \
               -query {input} \
               -outfmt 6 \
               -out {output}
        """

rule compute_similarity_profile:
    input:
        hgt_sim = join(OUT, 'orthofinder_blast', 'hgt_orthogroup_bact.blast'),
        bg_sim = join(OUT, 'orthofinder_blast', 'background_orthogroup_bact.blast')
    output: join(OUT, 'orthofinder_blast', 'similarity_profile_bact.svg')
    conda: "../envs/r.yaml"
    shell: "Rscript scripts/02_similarity_profile.R {input.hgt_sim} {input.bg_sim} {output}"


rule compute_exon_per_gene:
    input:
    output:
    shell:

# TODO: filter Ac_specifig genes based on similarity with bacteria
rule select_HGT_candidates:
    input: join(OUT, 'specific_genes', 'presence_compact.tsv')
    output: join(OUT, 'HGT_candidates.txt')
    shell: "cp {input} {output}"
