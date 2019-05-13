# Computing general statistics on input genome annotations.
# Usage: annots_stats input.gff output.pdf
# cmdoret, 20190509

library(tidyverse)
args <- commandArgs(trailingOnly=TRUE)

# input
in_gff <- args[1]

# output
out_plot <- args[2]

# load gff
gff <- read_tsv(in_gff, comment="#", col_names=F)
colnames(gff) <- c("chrom", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
gff <- gff %>% filter(type %in% c("mRNA", "CDS", "gene", "exon", "tRNA"))
# compute stats

clean_name <- strsplit(basename(in_gff), "\\.")[[1]]
# visualise
# Length distribution of different types
n_obs <- function(x){return(data.frame(y=0, label=paste0("n=", length(x))))}
pdf(out_plot)
ggplot(data=gff, aes(x=type, y=log10(end - start))) + 
    geom_violin(scale='width') + 
    stat_summary(fun.data = n_obs, geom = "text") +
    theme_bw() +
    xlab("") + 
    ylab("log10 length") + 
    ggtitle(clean_name)

dev.off()
# N exon / gene
