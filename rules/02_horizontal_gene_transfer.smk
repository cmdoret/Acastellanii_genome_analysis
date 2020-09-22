
# Retrieve HGT sequences from Neff v1 publication and translate them
rule get_v1_hgt_fa:
    output: join(OUT, 'hgt', 'NEFF_v1_hgt_cds.fa')
    params:
        hgt_ids = join(IN, 'misc', 'NEFF_v1_HGT.tsv'),
        cds = join(IN, 'cds', 'NEFF_v1.43.fa')
    conda: '../envs/seqkit.yaml'
    shell: "seqkit grep -r -f {params.hgt_ids} {params.cds} | seqkit translate > {output}"

# Same for the gff3 file...
rule get_v1_hgt_gff:
    output: join(OUT, 'hgt', 'NEFF_v1_hgt.tsv')
    params:
        gff = join(IN, 'annotations', 'NEFF_v1.43.gff'),
        hgt = join(IN, 'misc', 'NEFF_v1_HGT.tsv')
    run:
        hgt_ids = {i.strip() for i in open(params['hgt'])}
        outf = open(output[0], 'w')
        with open(params['gff'], 'r') as gff:
            hgt = False
            n_exons = 0
            for line in gff:
                fields = line.split('\t')
                if re.match('^#', line):
                    continue
                elif fields[2] == 'gene':
                    # ID chrom start end type strand n_exon gene_len mrna_len attributes
                    if hgt:
                        outf.write(f"{curr_hgt.strip()}\t{n_exons}\n")
                    gene_id = re.search('gene:([^;]+);', fields[8]).group(1)
                    if gene_id in hgt_ids:
                        hgt = True 
                        n_exons = 0
                        curr_hgt = line
                    else:
                        hgt = False
                elif fields[2] == 'exon':
                    n_exons += 1

# Build pangenome and protein tree from amoeba and bacteria, also including HGT inferred
# in the first genome paper
rule orthofinder:
    input: 
        org_proteomes = expand(
            join(OUT, 'filtered_proteomes', '{organism}.fa'),
            organism=organisms.clean_name
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
    input: join(OUT, 'hgt', 'orthofinder')
    output: 
        pres_mat = join(OUT, 'hgt', 'presence_matrix.tsv'),
        pres_compact = join(OUT, 'hgt', 'presence_compact.tsv')
    params:
        focus_st = samples.strain,
        amoeba = organisms.clean_name[organisms['type'] == 'amoeba'],
        bact = organisms.clean_name[organisms['type'] == 'bacteria']
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
        bact_df = ortho.loc[:, list(params['bact'])]
        # Binarize variables: Is there any gene for species c in each orthogroup
        any_pres = lambda df: df.apply(lambda c: ~ c.isnull(), axis=0).apply(np.any, axis=1)
        bact_pres = any_pres(bact_df)
        bact_abs = ~bact_pres
        amoeba_pres = any_pres(amoeba_df)
        amoeba_abs = ~amoeba_pres
        # Record  genes found only in amoeba but not bacteria
        # Note: hgt events inferred in Clarke et al. also included for comparison)
        vdf = pd.DataFrame({
            'neff'  : ~st_abs['Neff'],
            #'NEFF_v1_hgt_cds'  : ~get_absent(ortho.NEFF_v1_hgt_cds),
            'c3'    : ~st_abs['C3'],
            'ac'    : ~ac_abs,
            'bact'  : bact_pres,
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
# for now, groups will be hgt or background
rule get_spec_seq:
    input:
        pres = join(OUT, 'hgt', 'presence_compact.tsv'),
        orthofinder_dir = join(OUT, 'hgt', 'orthofinder')
    output:
        join(TMP, 'merged', '{group}_orthogroup_seq.fa')
    run:
        condition = {
            'hgt': 'ac and bact and not amoeba',
            'background': '(not bact) and amoeba or ac',
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


# Filter Ac. genes in HGT candidiates or out from background set
rule filter_ac_orthologues:
    input: join(TMP, 'merged', '{group}_orthogroup_seq.fa')
    output: join(TMP, 'merged', '{group}_orthogroup_seq_filtered.fa')
    params:
        id_filter = "^(C3|Neff)_FUN_",
        invert = lambda w: '-v' if w.group == 'background' else ''
    conda: '../envs/seqkit.yaml'
    shell: "seqkit grep {params.invert} -r -p '{params.id_filter}' {input} > {output}"


  

# Descriptive statistics of HGT candidates
rule compare_HGT_candidates:
    input:
        pres = join(OUT, 'hgt', 'presence_compact.tsv'),
        anno = expand(join(OUT, 'stats', '{strain}_annot_stats.tsv'), strain=["C3", "Neff"])
    params:
        og_to_gene = join(OUT, 'hgt', 'orthofinder', 'Results_amoeba', 'Orthogroups', 'Orthogroups.txt')
    output:
        plt = join(OUT, 'plots', 'hgt_stats.svg'),
        tbl = join(OUT, 'hgt', 'hgt_stats.tsv')
    script: "../scripts/02_validate_hgt.py"


rule get_new_hgt:
    input:
        tbl = join(OUT, 'hgt', 'hgt_stats.tsv'),
        fa = expand(join(TMP, 'renamed', '{strain}_cds.fa'), strain=samples.index)
    output: join(OUT, 'hgt', 'ac_hgt_cds.fa')
    params:
        hgt_names = temp(join(TMP, 'hgt', 'hgt_names.tsv')),
        all_cds = temp(join(TMP, 'hgt', 'all_cds.fa'))
    conda: '../envs/seqkit.yaml'
    shell:
        """
        mkdir -p $(dirname {params.all_cds})
        # Find all proteins which match HGT names in the table
        # HGT names are those where the last column (hgt) is equal to 1
        grep "1$" {input.tbl} | cut -f1 > {params.hgt_names}
        cat {input.fa} > {params.all_cds}
        seqkit grep -r -f {params.hgt_names} {params.all_cds} > {output}
        rm {params.hgt_names} {params.all_cds}
        """

# Blast v1 and v2 HGT against bacteria and amoeba backgrounds
rule blast_v1_v2:
    input:
        v2 = join(TMP, 'merged', 'hgt_orthogroup_seq_filtered.fa'),
        v1 = join(OUT, 'hgt', 'NEFF_v1_hgt_cds.fa'),
        amoeba = join(TMP, 'merged', 'background_orthogroup_seq_filtered.fa')
    output: join(OUT, 'hgt', 'hgt_similarity.blast')
    params:
        bact_db = join(DB, 'blast', 'nr_v5', 'nr_v5'),
        amoeba_db = join(TMP, 'amoeba_blast_db')
    threads: NCPUS
    conda: '../envs/blast.yaml'
    singularity: 'docker://cmdoret/blast:2.9.0'
    shell:
        """
        blastp -num_threads {threads} \
               -max_target_seqs 5 \
               -taxids 2 \
               -db {params.bact_db} \
               -query {input.v1} \
               -outfmt 6 \
            | awk -vOFS='\t' '{{print "bact","v1",$0}}' > {output} 

        blastp -num_threads {threads} \
               -max_target_seqs 5 \
               -taxids 2 \
               -db {params.bact_db} \
               -query {input.v2} \
               -outfmt 6 \
            | awk -vOFS='\t' '{{print "bact","v2",$0}}' >> {output} 
        
        makeblastdb -dbtype prot -in {input.amoeba} -title amoeba -out {params.amoeba_db}
        blastp -outfmt 6 -db {params.amoeba_db} -num_threads {threads} -query {input.v1} \
            | awk -vOFS='\t' '{{print "amoeba","v1",$0}}' >> {output} 
        blastp -outfmt 6 -db {params.amoeba_db} -num_threads {threads} -query {input.v2} \
            | awk -vOFS='\t' '{{print "amoeba","v2",$0}}' >> {output} 

        rm -f {params.amoeba_db}.p{{hr,in,sq}}
        """

rule blast_v1_hgt_vs_ac:
    input:
        v1_hgt = join(OUT, 'hgt', 'NEFF_v1_hgt_cds.fa')
    output:
        blast = join(OUT, 'blast', 'hgt_v1_vs_neff.blast')
    params:
        neff = samples.proteome["Neff"],
        neff_db = temp(join(TMP, 'blast', 'neff_db'))
    threads: 12
    shell:
        """
        makeblastdb -dbtype prot -in {params.neff} -title neff -out {params.neff_db}
        blastp -evalue 0.001 -outfmt 6 -db {params.neff_db} -num_threads {threads} -query {input.v1_hgt} > {output.blast}
        """


# GO enrichment analysis on HGT candidates
rule go_enrichment:
    input: join(OUT, 'hgt', 'hgt_stats.tsv')
    output: join(OUT, 'go_enrich', 'hgt_go_enrich.tsv')
    conda: '../envs/r.yaml'
    shell: 'Rscript scripts/go_enrich.R {input} {output}'

rule quant_xenolog_species:
    input:
        ortho_dir = join(OUT, 'hgt', 'orthofinder'),
        fastas = expand(
            join(OUT, 'filtered_proteome', '{organism}.fa'),
            organism=organisms.clean_name
        )
    output:
        plt = join(OUT, 'plots', 'xenolog_species_{strain}.svg'),
        tbl = join(OUT, 'hgt', 'xenolog_tbl_{strain}.tsv')
    run:
        # Generate a mapping from protein ID's to species name
        sp_dicts = []
        for fa in input['fastas']:
            clean_sp = basename(fa)[:-3]
            sp = organisms.loc[organisms.clean_name == clean_sp, 'name'].values[0]
            sp_dicts.append({prot.name: sp for prot in SeqIO.parse(fa, 'fasta')})
        # Flatten the list of dictionaries into a single dict
        prot2sp = {p: s for sp_dict in sp_dicts for p, s in sp_dict.items()}
        xeno_path = join(input['ortho_dir'], 'Results_amoeba', 'Putative_Xenologs', f'{wildcards.strain}.tsv')
        xeno = pd.read_csv(xeno_path, sep='\t')
        # Count the number of occurences of each donor in putative Xenolog for A. castellanii
        donors = defaultdict(list)
        for prot, og in zip(xeno.Other, xeno.Orthogroup):
            for donor_prot in prot.split(','):
                donor_prot = donor_prot.strip()
                if donor_prot.startswith(("Neff", "C3")):
                    continue
                # Append the orthogroup and protein names of the current protein
                donors[prot2sp[donor_prot]].append((og, donor_prot))
        labels, all_prots, all_ogs, genera = [], [], [], []
        for name, prots in donors.items():
            split = name.split(" ")
            if len(split) > 1:
                genus = split[0][0].upper()
                species = ' '.join(split[1:])
                labels.append(genus + ". " + species)
                genera.append(split[0])
            else:
                labels.append(name)
                genera.append(name)
            all_ogs.append([p[0].strip() for p in prots])
            all_prots.append([p[1].strip() for p in prots])
        donors_df = pd.DataFrame({'label': labels, 'orthogroups': all_ogs, 'prots': all_prots, 'genus': genera})
        mpl.use("Agg")
        fig, axes = plt.subplots(1, 2, figsize=(15, 12))
        # Count the number of different orthogroups in each genus and plot it
        (
            donors_df
            .groupby('genus')
            .orthogroups
            .apply(lambda g: len(np.unique([v for s in g for v in s])))
            .plot(kind='barh', ax=axes[0])
        )
        donors_df['n_prots'] = donors_df.prots.apply(lambda x: len(np.unique(x)))
        donors_df.plot(x='label', y='n_prots', kind='barh', ax=axes[1])
        axes[0].set_title('Number of distinct orthogroups with putative xenologs per genus.')
        axes[1].set_title('Number of putative xenologs per species.')
        plt.savefig(output['plt'])
        donors_df = donors_df.loc[:, ['genus', 'label', 'prots', 'orthogroups', 'n_prots']]
        donors_df.prots = donors_df.prots.apply(lambda x: ','.join(x))
        donors_df.orthogroups = donors_df.orthogroups.apply(lambda x: ','.join(x))
        donors_df.to_csv(output['tbl'], sep='\t', index=False)




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


# Visualise distribution of sequence similarity between HGT predicted in v1 and v2 
# vs bacterial genes
rule plot_similarity_hgt:
    input: join(OUT, 'hgt', 'hgt_similarity.blast')
    output: join(OUT, 'figures', 'hgt_similarity.svg')
    run:
        import seaborn as sns
        sim = pd.read_csv(input[0],
            sep='\t',
            header=None,
            usecols=[0, 1, 4],
            names=['comp', 'set', 'similarity']
        )
        sns.violinplot(data=sim, x='comp', y='similarity', split=True, hue='set')
        plt.savefig(output[0])


# Filter the best hit in the Neff proteome for every HGT
rule filter_blast_hits:
    input: join(OUT, 'blast', 'hgt_v1_vs_neff.blast')
    output: join(OUT, 'blast', 'hgt_v1_vs_neff_filtered.blast')
    run:
        df = pd.read_csv(input[0], sep='\t', header=None, usecols=[0, 1, 10], names=['query', 'target', 'evalue'])
        keep_idx = df.groupby('query')['evalue'].transform(min) == df['evalue']
        df[keep_idx].to_csv(output[0], sep='\t', index=False)


# Compute proportion of matching HGT between Clarke et al and our calls
rule plot_overlap_hgt_v1_v2:
    input:
        v1 = join(OUT, 'blast', 'hgt_v1_vs_neff_filtered.blast'),
        v2 = join(OUT, 'hgt', 'hgt_stats.tsv')
    output: join(OUT, 'plots', 'v1_v2_hgt_overlap.svg')
    run:
        v1 = pd.read_csv(input['v1'], sep='\t')
        v2 = pd.read_csv(input['v2'], sep='\t' )
        v2 = v2.loc[v2.ID.str.contains("Neff"), :]
        v2.ID = v2.ID.str.replace("Neff_", "")
        v1.target = v1.target.str.replace('-T[0-9]+$', '')
        merged = v1.merge(v2, left_on="target", right_on="ID", how='outer')
        mpl.use('Agg')
        fig, ax = plt.subplots(1, 2)
        venn2(
            subsets = (
                ((merged.hgt != 1) & ~merged['query'].isnull()).sum(), #v1 only
                (merged['query'].isnull() & (merged.hgt == 1)).sum(), #v2 only
                ((merged.hgt == 1) & ~merged['query'].isnull()).sum(), #both v1 and v2
            ),
        set_labels = ("HGT_v1", "HGT_v2"),
        ax=ax[0]
        )
        ax[0].set_title("Overlap between our HGT calls and Clarke et al.")
        vdf = merged.loc[~merged['query'].isnull(), ['ac', 'bact', 'amoeba']].astype(bool)
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
            ax=ax[1]
        )
        ax[1].set_title("Group distribution of HGT calls from Clarke et al.")
        plt.savefig(output[0])
