# Computing general statistics on input genome annotations.
# Usage: annots_stats input.gff output.pdf
# cmdoret, 20190509

library(readr)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(gridExtra)
args <- commandArgs(trailingOnly=TRUE)

### GET ARGS ###
v1_gff <- args[1]
neff_gff <- args[2]
c3_gff <- args[3]
out_tbl <- args[4]
out_plot <- args[5]

### LOAD DATA ###
load_gff <- function(in_gff, name){
    gff <- read_tsv(in_gff, comment="#", col_names=F)
    colnames(gff) <- c("chrom", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
    gff <- gff %>% filter(type %in% c("mRNA", "CDS", "gene", "exon"))

    ### Clean data
    gff <- gff %>% mutate(attributes=gsub("Parent=", "ID=", attributes)) %>%
        mutate(ID=str_extract(attributes, 'ID=[^;\\.]*')) %>%
        mutate(ID=gsub("=[a-zA-Z]{1,}:", "=", ID)) %>%
        mutate(ID=sapply(ID, function(x){str_split(x,regex('[=-]'))[[1]][2]})) %>%
        mutate(name=name)
    return(gff)
}
gff <- load_gff(v1_gff, "NEFF_v1") %>%
    bind_rows(load_gff(neff_gff, "Neff")) %>%
    bind_rows(load_gff(c3_gff, "C3"))

### COMPUTE STATS ###

gff_gene_len <- gff %>%
    filter(type == 'gene') %>%
    mutate(gene_len = end - start)

# N exon / gene (drop duplicated exons)
gff_exons <- gff %>% 
    filter(type != 'gene') %>%
    group_by(ID) %>%
    distinct(chrom, start, end, type, .keep_all=T) %>%
    mutate(n_exon=sum(type == 'exon')) %>%
    filter(type == 'mRNA') %>% 
    mutate(mrna_len = end - start)

# Combine informations and reorder columns
gff_stats <- gff_gene_len %>%
    select(ID, gene_len) %>%
    inner_join(gff_exons, by='ID') %>%
    select(ID, chrom, start, end, type, strand, n_exon, gene_len, mrna_len, attributes)


### VISUALISE ###

# Label generator to get number of samples
n_obs <- function(x){return(data.frame(y=0, label=paste0("n=", length(x))))}

feature_len <- ggplot(data=gff %>% filter(type %in% c("CDS", "gene")), 
                      aes(x=name, y=log10(end - start), fill=name)) +
    geom_violin() +
    stat_summary(fun.data = n_obs, geom = "text") +
    theme_minimal() +
    xlab("") + 
    ylab("log10 length") + 
    ggtitle("Gene length") +
    facet_wrap(.~type) +
    scale_y_continuous(limits = c(0, 6)) + 
    geom_boxplot(width=0.1)

exon_quantiles <- gff_exons %>% 
  group_by(name) %>% 
  summarise(x=list(enframe(quantile(n_exon, probs=c(0.25,0.5,0.75)), "quantiles", "exons"))) %>% 
  unnest(x) %>%
  spread(key=quantiles, value=exons)

# Assign a quantile label to each exon number based on its strain's exon number distribution
exons_tbl <- gff_exons %>%
    inner_join(exon_quantiles, by="name") %>%
    mutate(
        n_exon_quant = case_when(
            n_exon <= `25%` ~ "25% quantile",
            n_exon >= `75%` ~ "75% quantile",
            TRUE            ~ "Median"
        )
    ) %>%
    select(name, n_exon, n_exon_quant)

exons_tbl <- exons_tbl %>% 
    mutate(
        n_exon_quant = factor(
            n_exon_quant,
            ordered=T, 
            levels=c("25% quantile", "Median", "75% quantile")
        )
    )

exon_per_gene <- ggplot(data=exons_tbl, aes(x=n_exon, fill=n_exon_quant)) +
    scale_fill_manual(values=c("#6666ee", "#999999", "#ee6666")) +
    geom_histogram(binwidth=1) +
    theme_minimal() +
    stat_summary(aes(x = 0.1, y = n_exon, xintercept = stat(y), group = name), 
               fun.y = median, geom = "vline") +
    ylab("Genes") +
    xlab("Exons per genes") +
    facet_grid(name~.) +
    scale_x_continuous(limits = c(-0, 100))

### SAVE OUTPUT ###
svg(out_plot, height=10, width=12)
grid.arrange(feature_len, exon_per_gene)
dev.off()
write_tsv(gff_stats, out_tbl)

