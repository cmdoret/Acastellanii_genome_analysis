# Visualise blast similarity against bacteria from target and background gene families.

library(readr)
library(ggplot2)
library(dplyr)
library(ggpubr)

args <- commandArgs(trailingOnly=T)

target <- read_tsv(args[1], col_names=F)
background <- read_tsv(args[2], col_names=F)
out_figure <- args[3]

target$gene_group <- paste0('A. castellanii-specific (', nrow(target), ")")
background$gene_group <- paste0('background (', nrow(background), ")")
all_hits <- bind_rows(target, background)

# pident is column 3
svg(out_figure)
gghistogram(all_hits, x = "X3",
    y = '..density..',
    add = "median", rug = F,
    color = "gene_group", fill = "gene_group",
    palette = c("#00AFBB", "#E7B800")) + 
    xlab("BLAST identity") + 
    ylab("Density") +
    ggtitle("Similarity of A. castellanii-specific genes to bacteria")
dev.off()