#!/usr/bin/env python3
# Given a list of HGTs and a table of annotation statistics, compute and visualize
# exploratory stats on the HGT and non-HGT sets.
# cmdoret, 20200215
import re
import seaborn as sns
import pandas as pd
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt

# Read input orthogroup presence matrix
pres = pd.read_csv(snakemake.input["pres"], sep="\t")
# Read all input annotation stats and concatenate into a single df
anno = pd.concat(
    [pd.read_csv(m, sep="\t") for m in snakemake.input["anno"]], axis=0
)

og2g = {}
g2og = {}
with open(snakemake.params["og_to_gene"], "r") as ortho_out:
    for line in ortho_out:
        og, genes = line.split(":")
        genes = genes.split(" ")
        # Trim newline character and funannotate-specific suffix to match
        # orthofinder output with funnanotate annotation ids. Also exclude potential
        # empty strings with if filter
        genes = [re.sub("-T1", "", g.rstrip("\n")) for g in genes if g]
        og2g[og] = genes
        for g in genes:
            g2og[g] = og

# map each gene to it's orthogroup
anno["Orthogroup"] = anno.ID.apply(lambda g: g2og.get(g, None))
# Define HGT based on presence/absence
pres["hgt"] = (pres.bact) & (pres.ac) & (~pres.amoeba)
anno = anno.merge(pres, on="Orthogroup", how="inner")


st.ttest_ind(anno.n_exon[anno.hgt == 1], anno.n_exon[anno.hgt == 0]).pvalue
anno[["hgt", "n_exon", "gene_len", "mrna_len"]].groupby("hgt").agg("mean")
anno_tmp = anno.copy()
anno_tmp["n_exon"] = st.zscore(anno_tmp["n_exon"], nan_policy="omit")
anno_tmp["gene_len"] = st.zscore(anno_tmp["gene_len"], nan_policy="omit")
long_anno = anno_tmp[["ID", "hgt", "gene_len", "n_exon"]].melt(["hgt", "ID"])
sns.violinplot(
    data=long_anno,
    x="variable",
    y="value",
    inner="quartiles",
    split=True,
    hue="hgt",
)
plt.ylabel("z-score")
plt.title(
    "Comparison of A. castellanii HGT candidates versus background genes"
)
plt.savefig(snakemake.output[0])
