"""
Computes the gap-excluded sequence divergence from a PAF file
cmdoret, 202009
"""

import sys
import re


def get_cigar_gaps(cigar):
    # Match all indels and capture indel size
    gap_regex = re.compile(r"([0-9]+)[ID]")
    total_gap_size = 0
    for gap_size in re.findall(gap_regex, cigar):
        total_gap_size += int(gap_size)
    return total_gap_size


with open(sys.argv[1]) as paf, open(sys.argv[2], "w") as outf:
    for line in paf:
        align = line.split("\t")
        if align[16] != "tp:A:P":  # Only consider primary alignments
            continue
        cigar = align[-1]
        n_gaps = get_cigar_gaps(cigar)
        n_ambiguous = int(align[15].lstrip("nn:i:"))
        n_bad = int(align[12].lstrip("NM:i:"))
        n_matches = int(align[9])
        n_mismatches = n_bad - (n_ambiguous + n_gaps)
        div = n_mismatches / (n_matches + n_mismatches)
        outf.write(f"{align[5]}\t{align[7]}\t{align[8]}\t{div}\n")
