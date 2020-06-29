# Testing for enrichment of GO Terms among HGT candidates
# 20190604, cmdoret

library(tidyverse)
library(topGO)
library(viridis)

# Parse CL args
args <- commandArgs(trailingOnly=T)

# Extract gene IDs and GO terms.
annot_tbl <- read_tsv(args[1], col_types='cciicciiicciiiiiii') %>%
  dplyr::select(ID, attributes, hgt) %>%
  mutate(GO=str_replace(attributes, '.*Ontology_term=([^;]*);', "\\1"))
hgt_candidates <- annot_tbl %>%
  filter(hgt == 1) %>%
  pull(ID)

out_file <- args[2]
tmp_mapfile <- paste(dirname(out_file), 'id2go.tsv', sep='/')
out_fig <- paste0(tools::file_path_sans_ext(out_file), ".svg")

# Write geneID - GO terms mapping to file
annot_tbl %>%
  dplyr::select(ID, GO) %>%
  mutate(GO = str_replace_all(GO, ";", ", ")) %>%
  write_tsv(tmp_mapfile)


geneID2GO <- readMappings(tmp_mapfile)
geneNames <- names(geneID2GO)
geneList <- factor(as.integer(geneNames %in% hgt_candidates))
names(geneList) <- geneNames
GOdata <- new(
  "topGOdata",
  description = "HGT analyis in A. castellanii",
  ontology = "BP",
  allGenes = geneList, 
  nodeSize = 10,
  gene2GO = geneID2GO, 
  annot = annFUN.gene2GO
)


fisher.stat <- new(
  "classicCount",
  testStatistic = GOFisherTest, 
  name = "Fisher test"
)
resultFisher <- getSigGroups(GOdata, fisher.stat)
weight.stat <- new(
  "weightCount",
  testStatistic = GOFisherTest, 
  name = "Fisher test with weight algorithm",
  sigRatio = 'ratio'
)
resultWeight <- getSigGroups(GOdata, weight.stat)

#showSigOfNodes(GOdata, score(resultWeight), firstSigNodes = 9, useInfo = 'all')

# Get 'significant' terms
pvals <- score(resultWeight)
sig <- pvals[pvals<0.01]

termStat(GOdata, names(sig))

# Report top enriched GO terms
allRes <- GenTable(
  GOdata,
  classic = resultFisher,
  weight = resultWeight, 
  orderBy = "weight",
  ranksOf = "classic",
  topNodes = 200
)

write_tsv(allRes, out_file)
allRes <- allRes %>%
	mutate(GO.ID=paste(GO.ID, Term, sep=', '))

tidy_res = as.tibble(allRes) %>%
  mutate(
    weight = as.numeric(weight),
    classic = as.numeric(classic),
    GO.ID = factor(
      GO.ID,
      ordered=T,
      levels=unique(GO.ID[ordered(weight)])
    )
  )

# Visualise top enriched GO terms
svg(out_fig)
ggplot(
  data=tidy_res %>% top_n(30, -weight), 
  aes(x=GO.ID, y=-log10(weight))
) + 
  geom_segment(
    aes(xend=GO.ID, yend=min(-log10(weight))),
    size=1.1
  ) +
  geom_point(aes(size=Annotated, color=Significant / Expected)) + 
  geom_hline(aes(yintercept=2), lty=2, col='red') +
  theme_minimal() + 
  xlab("") +
  ylab("-log10 pvalue") +
  ggtitle("GO enrichment test in HGT candidates\n(Fisher exact test, weight algorithm)") +
  scale_color_viridis() +
  coord_flip()
dev.off()
