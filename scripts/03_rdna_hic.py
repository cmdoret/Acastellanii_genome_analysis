# Visualize Hi-C matrices and overlay rDNA positions on top
import sys
import os.path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches as mpatches
from scipy.ndimage import gaussian_filter
import cooler
from scipy.sparse import coo_matrix
import hicstuff.hicstuff as hcs
import matplotlib as mpl
import matplotlib.gridspec as gridspec

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

plt.style.use("seaborn-white")
gs = gridspec.GridSpec(4, 1, height_ratios=[5, 1, 1, 1])
mpl.rcParams["figure.figsize"] = (5, 10)
fig = plt.figure()
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])
ax3 = fig.add_subplot(gs[2])
ax4 = fig.add_subplot(gs[3])
fig.align_xlabels()
ax2.margins(x=0)
ax3.margins(x=0)
ax4.margins(x=0)

coldict = {
    "18s_rRNA": ("b", "-.", "18S", ax3),
    "28s_rRNA": ("g", ":", "28S", ax2),
    "8s_rRNA": ("r", "--", "5S", ax4),
}
legend_entries = [
    mpatches.Patch(color=v[0], label=v[2], ls=v[1]) for k, v in coldict.items()
]
plt.legend(handles=legend_entries, loc="upper right")
# mat[np.isnan(mat)] = 0
ax1.imshow(np.log(mat), cmap="Reds")

# Add vertical lines to line plots at rdna positions
for r in gff.iterrows():
    coldict[r[1]["attribute"]][3].axvline(
        x=r[1]["maxbin"],
        alpha=0.4,
        # linestyle=coldict[r[1]["attribute"]][1],
        lw=0.3,
        c=coldict[r[1]["attribute"]][0],
    )

# Add chromosome grid
for r in chroms.iterrows():
    ax2.axvline(x=c.extent(r[1]["name"])[1], alpha=0.4, c="grey", lw=0.5)
    ax3.axvline(x=c.extent(r[1]["name"])[1], alpha=0.4, c="grey", lw=0.5)
    ax4.axvline(x=c.extent(r[1]["name"])[1], alpha=0.4, c="grey", lw=0.5)


# Subset bins containing rRNA
rdna_bins = {
    sub: np.unique(gff.maxbin[gff.attribute == sub]) for sub in np.unique(gff.attribute)
}
# Compute contact sum per bin for each subunit
rdna_sig = {}
rdna_sig["28S"] = mat[rdna_bins["28s_rRNA"], :].mean(axis=0)
rdna_sig["18S"] = mat[rdna_bins["18s_rRNA"], :].mean(axis=0)
rdna_sig["5S"] = mat[rdna_bins["8s_rRNA"], :].mean(axis=0)

ax2.plot(range(mat.shape[0]), np.log10(rdna_sig["28S"]), label="28S", lw=0.5, c="g")
ax3.plot(range(mat.shape[0]), np.log10(rdna_sig["18S"]), label="18S", lw=0.5, c="b")
ax4.plot(range(mat.shape[0]), np.log10(rdna_sig["5S"]), label="5S", lw=0.5, c="r")

fig.savefig(out_path, dpi=900)
