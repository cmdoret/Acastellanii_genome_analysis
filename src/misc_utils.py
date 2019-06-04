# Helper or aesthetics functions
# cmdoret, 20190502
import sys
import numpy as np
import pandas as pd


def progbar(curr, tot, status="", size=40):
    """Displays a progress bar in the prompt"""
    filled = int(round(size * curr / float(tot)))

    percents = round(100.0 * curr / float(tot), 1)
    bar = "=" * filled + "-" * (size - filled)

    sys.stdout.write("  [%s] %s%s - %s\r" % (bar, percents, "%", status))
    sys.stdout.flush()


def bin_tsv_file(
    bed_in, bed_out, bin_size=1000, chrom=0, start=1, end=2, val=3, fun=np.mean
):
    """
    Bins a BED-like file based on genomic position, applying a user
    defined on a value column.

    Parameters
    ----------
    bed_in : str
        Path to the input BED-like file.
    bed_out : str
        Path to the output binned file
    bin_size : int
        Bin size to use, in base pairs.
    chrom : int
        0-based index of the chromosome column.
    start : int
        0-based index of the start column.
    end : int
        0-based index of the end column.
    val : int
        0-based index of the value column on which the function will be applied.
    fun : function
        A valid function to apply on the value column. Could be a lambda function
        or a numpy function.
    """
    bed = pd.read_csv(bed_in, sep="\t")
    bed["window"] = bed.iloc[:, start] // bin_size
    chrom_col, val_col = bed.columns[[chrom, val]]
    agg_val = (
        bed.groupby([chrom_col, "window"], sort=False)[val_col].agg(fun).reset_index()
    )
    agg_pos = (
        bed.groupby([chrom_col, "window"], sort=False)
        .agg({bed.columns[start]: min, bed.columns[end]: max})
        .reset_index()
    )
    new_val_col = val_col + "_" + str(bin_size) + "bp_bins"
    agg_pos[new_val_col] = agg_val[val_col].apply(lambda x: np.round(x, 4))

    agg_pos.to_csv(bed_out, sep="\t", index=False)

