import pandas as pd
from matplotlib_venn import venn2, venn3
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.use("Agg")
vdf = pd.read_csv(snakemake.input["pres_compact"], sep="\t").drop(
    "Orthogroup", axis=1
)
vdf = vdf.astype(bool)
# vdf = vdf.rename(columns={'NEFF_v1_hgt_cds': 'v1'})
# Make Venn diagram. sets are A, B, AB, C, AC, BC, ABC
fig, ax = plt.subplots(1, 2, figsize=(15, 12))
venn3(
    subsets=(
        len(vdf.query("    c3 and not neff and not amoeba")),  # A
        len(vdf.query("not c3 and     neff and not amoeba")),  # B
        len(vdf.query("    c3 and     neff and not amoeba")),  # AB
        len(vdf.query("not c3 and not neff and     amoeba")),  # C
        len(vdf.query("    c3 and not neff and     amoeba")),  # AC
        len(vdf.query("not c3 and     neff and     amoeba")),  # BC
        len(vdf.query("    c3 and     neff and     amoeba")),  # ABC
    ),
    set_labels=("C3", "Neff", "Amoeba"),
    ax=ax[0],
)
venn2(
    subsets=(
        len(vdf.query("    c3 and not neff")),
        len(vdf.query("not c3 and     neff")),
        len(vdf.query("    c3 and     neff")),
    ),
    set_labels=("C3", "Neff"),
    ax=ax[1],
)
fig.savefig(snakemake.output["venn"])

