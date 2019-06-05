# Testing for enrichment of GO Terms among HGT candidates
# 20190604, cmdoret

library(tidyverse)
library(topGO)

# Parse CL args
args <- commandArgs(trailingOnly=T)

# Extract gene IDs and GO terms.
annot_tbl <- read_tsv(args[1], col_types='ccciiccccccccccccccc')
# Write geneID - GO terms mapping to file
write_tsv(annot_tbl, tmp_dir)
hgt_candidates <- read_tsv(args[2]) %>%
  pull(GeneID)
out_file <- args[4]
tmp_mapfile <- paste(dirname(out_file), 'id2go.tsv', sep='/')

annot_tbl %>%
  dplyr::select(GeneID, `GO Terms`) %>%
  mutate(`GO Terms` = str_replace_all(`GO Terms`, ";", ", ")) %>%
  write_tsv(tmp_mapfile)


geneID2GO <- readMappings(tmp_mapfile)

GOdata <- new("topGOdata",
              description = "HGT analyis in A. castellanii",
              ontology = "BP",
              allGenes = annot_tbl %>% pull(GeneID), 
              geneSel = hgt_candidates,
              nodeSize = 10,
              gene2GO = geneID2GO)

showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = 'all')
