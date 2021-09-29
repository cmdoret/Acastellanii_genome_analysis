# Generate a nice BUSCO visualization from the BUSCO table
from collections import OrderedDict
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt


def read_busco(tbl: str, name: str) -> pd.DataFrame:
    """
    Read a busco summary table into a
    dataframe and create relevant columns
    """
    df = pd.read_csv(tbl, names=["number", "type"], sep="\t", header=None)
    df["strain"] = name
    # Single letter busco class as a column
    df["busco"] = df["type"].str.extract(r".*\(([A-Z])\).*")
    df.loc[pd.isnull(df.busco), "busco"] = "T"
    # Transform number of buscos to percentages
    df["perc"] = 100 * np.round(
        df["number"] / df.number[df.busco == "T"].values[0], decimals=2
    )
    return df


# Concatenate busco tables from all strains and set relevant
# variables as indices
busco = (
    pd.concat([read_busco(tbl, name) for name, tbl in snakemake.input.items()])
    .reset_index(drop=True)
    .sort_values("strain", ascending=False)
    .set_index(["busco", "strain"])
)
mpl.use("Agg")
# Keep track of the order in which strains will be plotted
str_to_num = OrderedDict({"v1": 0, "neff": 1, "c3": 2})
# Map each BUSCO class to a color
col_dict = OrderedDict(
    [("M", "#F04441"), ("F", "#F0E441"), ("D", "#3592C6"), ("S", "#58B4E8"),]
)
# We plot one busco class per iteration
for i, t in enumerate(col_dict.keys()):
    # Subset busco table for current strain
    busco_type = busco.loc[t]
    # Retrieve x index from strains (deterministic order)
    r = [str_to_num[s] for s in busco_type.index.values]
    # Stacked barplot, uses cumulative height of previous bars
    if i:
        plt.barh(r, busco_type.perc.values, color=col_dict[t], left=cum_height)
        cum_height += busco_type.perc.values
    else:
        cum_height = busco_type.perc.values
        plt.barh(r, busco_type.perc.values, color=col_dict[t])


def busco_string(strain):
    """Use the BUSCO summary table to construct the standard busco string"""
    mis = busco.loc[("M", strain), "number"]
    fra = busco.loc[("F", strain), "number"]
    com = busco.loc[("C", strain), "number"]
    sin = busco.loc[("S", strain), "number"]
    dup = busco.loc[("D", strain), "number"]
    tot = busco.loc[("T", strain), "number"]
    return f"{strain} - C:{com}[S:{sin}, D:{dup}], F:{fra}, M:{mis}, n={tot}"


# Write the strain name and busco string besides each bar
plt.yticks(r, [busco_string(s) for s in str_to_num.keys()])
plt.xlabel("% BUSCOs")
plt.ylabel("Strain")
plt.savefig(str(snakemake.output))
