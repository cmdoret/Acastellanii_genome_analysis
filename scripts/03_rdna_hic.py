# Visualize Hi-C matrices and overlay rDNA positions on top
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches as mpatches
from scipy.ndimage import gaussian_filter
import cooler
from scipy.sparse import coo_matrix
import hicstuff.hicstuff as hcs

## LOAD CLI ARGS
cool_path = sys.argv[1]
rdna_path = sys.argv[2]
out_path = sys.argv[3]

## LOAD INPUT FILES
c = cooler.Cooler(cool_path)
gff = pd.read_csv(rdna_path, sep="\t", comment="#", usecols=range(9))
gff.columns = [
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute",
]

## EXTRACT COOL DATA
mat = c.pixels()[:]
bins = c.bins()
chroms = c.chroms()[:]

## TRANSFORM

# Transform rDNA coords to bins
gff["minbin"] = gff.apply(
    lambda r: np.min(bins.fetch(f"{r.seqname}:{r.start}-{r.end}").index.values), axis=1
)

gff["maxbin"] = gff.apply(
    lambda r: np.max(bins.fetch(f"{r.seqname}:{r.start}-{r.end}").index.values), axis=1
)
gff = gff.drop_duplicates(subset=["minbin", "maxbin"])

## VISUALIZE
mat = c.matrix(balance=False)[:]
np.fill_diagonal(mat, 0)

# Define color saturation threshold
# sat_thr = np.percentile(mat[mat > 0], 85)

# Plot heatmap with rDNA overlay
fig = plt.figure()
coldict = {"18s_rRNA": ("b", "-."), "28s_rRNA": ("g", ":")}  # "8s_rRNA": ("r", "--")}
legend_entries = [
    mpatches.Patch(color=v[0], label=k, ls=v[1]) for k, v in coldict.items()
]
plt.legend(handles=legend_entries, loc="lower right")
# mat[np.isnan(mat)] = 0
plt.imshow(np.log(mat), cmap="Reds")

for r in gff.iterrows():
    if r[1]["attribute"] != "8s_rRNA":
        plt.axvline(
            x=r[1]["maxbin"],
            alpha=0.4,
            linestyle=coldict[r[1]["attribute"]][1],
            lw=2,
            c=coldict[r[1]["attribute"]][0],
        )
        plt.axhline(
            y=r[1]["maxbin"],
            alpha=0.4,
            linestyle=coldict[r[1]["attribute"]][1],
            lw=2,
            c=coldict[r[1]["attribute"]][0],
        )

# Add chromosome grid

# for r in chroms.iterrows():
#    plt.axvline(x=c.extent(r[1]["name"])[1], alpha=0.4, c="black", lw=0.5)
#    plt.axhline(y=c.extent(r[1]["name"])[1], alpha=0.4, c="black", lw=0.5)

fig.savefig(out_path, dpi=900)
