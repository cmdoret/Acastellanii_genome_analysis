# Plot density of foreign origin genes along genome
# cmdoret, 20190521

library(tidyverse)
library(circlize)

args <- commandArgs(trailingOnly=TRUE)

### LOAD DATA ###
annot <- args[1]
out_fig <- args[2]

### FORMAT DATA ###
# Get approx contig size from annotations
contigs <- annot %>%
    group_by(Contig) %>%
    summarise(size=max(Stop))

good_contigs <- contigs$Contig[contigs$size > 5e5]
# remove genes with no known domains
filtered <- annot %>% 
    filter(n_all > 0) %>%
    group_by(Contig) %>%
    mutate(size=max(Stop)) %>%
    ungroup() %>%
    mutate(Contig=factor(Contig, 
                         ordered=T, 
                         levels=unique(filtered$Contig[
                             order(filtered$size, 
                                   decreasing=T)]))) %>%
    filter(size > 5e5) %>% 
    mutate(ori_bac=ifelse(prop_bac > 0.75, 1, 0),
           ori_vir=ifelse(prop_vir > 0.75, 1, 0))
    

# Running mean using stat_smooth with degree=0
ggplot(data=filtered, aes(x=(Start+Stop)/2)) +
    facet_wrap(~Contig) +
    geom_point(aes(y=100*prop_bac, size=n_all), alpha=0.3) +
    stat_smooth(aes(y=100*prop_bac, col="% bacterial domains"), 
                degree=0,method=loess, size=1.5, span=0.01, se=F) +
    stat_smooth(aes(y=100*prop_vir, col="% viral domains"), 
                degree=0,method=loess, size=1.5, span=0.01, se=F) +
    theme_light() +
    xlab("Genomic position [bp]") +
    ylab("Percentage of domains") +
    ggtitle("Origin of protein domains in A. castellanii Neff", 
            "On contigs larger than 500kb using a sliding mean of 10% genes")

ggsave(out_fig)