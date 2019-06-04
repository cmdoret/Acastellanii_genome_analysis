# Testing for enrichment of GO Terms among HGT candidates
# 20190604, cmdoret

library(tidyverse)
library(topGO)

# Parse CL args
args <- commandArgs(trailingOnly=T)

annot_tbl <- rea
candidates <-
out <-

# Input format is:
# 133103  GO:0015031, GO:0005794, GO:0016020, GO:0017119, GO:0000139
# 121005  GO:0005576
# 155158  GO:0005488
# 160828  GO:0005488
# 105778  GO:0016021, GO:0016020
# 166881  GO:0003674, GO:0016021, GO:0016020, GO:0008150


geneID2GO <- readMappings(map_file)

sampleGOdata <- new("topGOdata",
                    description = "HGT analyis in A. castellanii",
                    ontology = "BP",
                    allGenes = geneList, 
                    geneSel = topDiffGenes,
                    nodeSize = 10,
                    annot = annFUN.db,
                    affyLib = affyLib)