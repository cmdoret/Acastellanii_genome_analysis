# Retrieve informations about laterally transferred genes in Neff
import re

hgt_ids = {i.strip() for i in open(snakemake.params["hgt"])}
outf = open(snakemake.output[0], "w")
with open(snakemake.params["gff"], "r") as gff:
    hgt = False
    n_exons = 0
    for line in gff:
        fields = line.split("\t")
        if re.match("^#", line):
            continue
        elif fields[2] == "gene":
            # ID chrom start end type strand n_exon gene_len mrna_len attributes
            if hgt:
                outf.write(f"{curr_hgt.strip()}\t{n_exons}\n")
            gene_id = re.search("gene:([^;]+);", fields[8]).group(1)
            if gene_id in hgt_ids:
                hgt = True
                n_exons = 0
                curr_hgt = line
            else:
                hgt = False
        elif fields[2] == "exon":
            n_exons += 1
