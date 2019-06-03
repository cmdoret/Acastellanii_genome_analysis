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
dias_cutoff <- unname(quantile(sighunt$dias, 0.95))
sighunt_ranges <- makeGRangesFromDataFrame(sighunt %>% filter(dias >= dias_cutoff))
interpro_ranges <- makeGRangesFromDataFrame(interpro %>% filter(prop_bac >= 0.85))
overlaps <- GenomicRanges::findOverlaps(interpro_ranges, sighunt_ranges)

candidates <- interpro[queryHits(overlaps), ] %>% 
    distinct(GeneID, chrom, Start, Stop, GeneID) %>% 
    arrange(chrom, Start) %>%
    select(chrom, Start, Stop, GeneID)

write_tsv(candidates, args[3])