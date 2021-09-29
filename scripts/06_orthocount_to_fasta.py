import sys
import pandas as pd

gene_count = sys.argv[1]
out_fasta = sys.argv[2]

gc = pd.read_csv(gene_count, sep="\t")


def get_presence(freqs):
    """
    Given an array of frequences, return gene presence code

    Parameters
    ----------
    freqs : numpy array of ints
        1D Array of gene frequences.

    Returns
    -------
    str :
        Gene presence code, in the form of a string of 1 (present) and 0 (absent)

    Examples
    --------
    >>> get_presence(np.array([32, 1, 0, 0, 2]))
    '11001'
    """
    freqs[freqs > 1] = 1
    freqs = freqs.astype(str)
    return "".join(freqs)


# Store species names and associated ortholog presence code in a dict
sp_presence = {
    gc.columns[i]: get_presence(gc.iloc[:, i]) for i in range(1, gc.shape[1] - 1)
}

with open(out_fasta, "w") as fa:
    for sp, presence in sp_presence.items():
        fa.write(">" + str(sp) + "\n")
        fa.write(presence + "\n")

