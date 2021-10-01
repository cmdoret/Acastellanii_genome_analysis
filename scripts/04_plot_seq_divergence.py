import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

plt.figure(figsize=(12, 12))
div = pd.read_csv(
    snakemake.input[0], sep="\t", names=["chrom", "start", "end", "div"]
)
div = div.sort_values(["chrom", "start"])
div["align_len"] = abs(div.end - div.start)
div["log_len"] = np.log10(div.align_len)
jp = sns.jointplot(
    data=div,
    x="log_len",
    y="div",
    kind="kde",
    shade=True,
    cmap="plasma",
    levels=200,
)
jp.set_axis_labels(
    "Log10 length of aligned block [bp]",
    "Per base divergence in aligned blocks",
)
avg_div = np.average(div["div"], weights=div["align_len"])
avg_log_len = np.mean(div.log_len)
jp.ax_joint.axvline(
    x=avg_log_len, linestyle="--", label=f"",
)
jp.ax_joint.axhline(y=avg_div, linestyle="--")
jp.ax_marg_x.axvline(x=avg_log_len)
jp.ax_marg_y.axhline(y=avg_div)
jp.fig.suptitle(
    "Neff - C3 gap excluded sequence divergence\n"
    f"mean divergence length: {100*avg_div:.2f}%"
)
plt.savefig(snakemake.output[0])
