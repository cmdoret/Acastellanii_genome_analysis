# Computing general statistics on input genome annotations.
# Usage: annots_stats input.gff output.pdf
# cmdoret, 20190509

library(tidyverse)
args <- commandArgs(trailingOnly=TRUE)

# input
gff <- args[1]

# output
out_plot <- args[2]

# load gff
gff <- read_tsv('data/input/annotations/amoeba/Neff.gff', comment="#", col_names=F)
colnames(gff) <- c("chrom", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
# compute stats

# visualise
# Length distribution of different types
pdf(out_plot)
ggplot(data=gff, aes(x=type, y=log(end - start))) + 
    geom_violin()
dev.off()
# N exon / gene
