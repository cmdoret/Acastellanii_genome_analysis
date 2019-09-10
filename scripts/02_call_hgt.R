# Combine results from sighunt and interpro hits to get HGT candidate genes
# cmdoret, 20190602

library(tidyverse)
library(GenomicRanges)

# Parse CL args
args <- commandArgs(trailingOnly=T)
sighunt <- read_tsv(args[1])
interpro <- read_tsv(args[2]) %>%
    dplyr::rename(chrom=Contig)

# Overlap candidates genes (85% of hits in prokaryotes) with sighunt 
# windows in top 5% highest DIAS
dias_cutoff <- unname(quantile(sighunt$dias, 0.90))
sighunt_ranges <- makeGRangesFromDataFrame(sighunt %>% filter(dias >= dias_cutoff))
interpro_ranges <- makeGRangesFromDataFrame(interpro %>% filter(prop_bac >= 0.85),
                                            keep.extra.columns=T)
overlaps <- subsetByOverlaps(interpro_ranges, sighunt_ranges)

# Extract positional data from Granges into a dataframe
candidates <- data.frame(chrom=seqnames(overlaps),
  start=start(overlaps)-1,
  stop=end(overlaps))

# Add metadata (incl gene ID), sort genes and order columns
candidates <- candidates %>%
    cbind(mcols(overlaps)) %>%
    as_tibble %>%
    distinct(GeneID, chrom, start, stop) %>%
    arrange(chrom, start) %>%
    select(chrom, start, stop, GeneID)

write_tsv(candidates, args[3])