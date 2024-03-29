
# Remove isoforms from Acas proteins to keep only the longest version
# of each gene
rule select_primary_transcript:
    input: join(TMP, 'renamed', '{strain}_proteome.fa')
    output: join(TMP, 'renamed', '{strain}_proteome_filtered.fa')
    shell:
        """
        # Convert fasta to tabular, sort by seqname, keep longest seq
        # and convert back to fasta
        sed 's/>[^ ]* \(.*\)$/>\\1\t/' {input} \
            | tr -d '\\n' \
            | sed 's/>/\\n>/g' \
            | sed '/^$/d' \
            | awk -vFS='\t' -vOFS='\t' '{{print $2,length($2),$1}}' \
            | sort -k3,3 -k2,2n \
            | uniq -f2 \
            | awk '{{print $3,$1}}' \
            | sed 's/ /\\n/' \
            > {output}
        """


# Build pangenome and protein tree from amoeba and bacteria, also including HGT inferred
# in the first genome paper
rule orthofinder:
    input: 
        org_proteomes = expand(
            join(OUT, 'filtered_proteomes', '{organism}.fa'),
            organism=organisms.clean_name[organisms['type'] == 'amoeba']
        ),
        ac_proteomes = expand(
            join(TMP, 'renamed', '{strain}_proteome_filtered.fa'),
            strain=samples.strain
        )
        #hgt_neff_v1 = join(OUT, 'hgt', 'NEFF_v1_hgt_cds.fa')
    output: directory(join(OUT, 'hgt', 'orthofinder'))
    threads: NCPUS
    params:
      orthofinder_dir = join(TMP, 'orthofinder')
    conda: '../envs/orthofinder.yaml'
    priority: 10
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
            new_fname="${{fname/_proteome_filtered/}}"
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
        ortho_dir = join(input[0], "Results_amoeba", "Orthogroups")
        ortho = pd.read_csv(join(ortho_dir, 'Orthogroups.tsv'), sep='\t')
        unassigned = pd.read_csv(join(ortho_dir, 'Orthogroups_UnassignedGenes.tsv'), sep='\t')
        # Add unassigned genes as self-contained orthogroups
        ortho = pd.concat([ortho, unassigned], axis=0)
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
    conda: '../envs/viz.yaml'
    script: '../scripts/06_plot_orthogroups_venn.py'
