# Computing general statistics on input genome annotations.
# Usage: annots_stats input.gff output.pdf
# cmdoret, 20190509

library(readr)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
args <- commandArgs(trailingOnly=TRUE)

### GET ARGS ###
in_gff <- args[1]
out_tbl <- args[2]
out_plot <- args[3]

### LOAD DATA ###
gff <- read_tsv(in_gff, comment="#", col_names=F)
colnames(gff) <- c("chrom", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
gff <- gff %>% filter(type %in% c("mRNA", "CDS", "gene", "exon"))

### Clean data
gff <- gff %>% mutate(attributes=gsub("Parent=", "ID=", attributes)) %>%
    mutate(ID=str_extract(attributes, 'ID=[^;]*')) %>%
    mutate(ID=gsub("=[a-zA-Z]{1,}:", "=", ID)) %>%
    mutate(ID=sapply(ID, function(x){str_split(x,regex('[=-]'))[[1]][2]}))

### COMPUTE STATS ###

gff_gene_len <- gff %>%
    filter(type == 'gene') %>%
    mutate(gene_len = end - start)

# N exon / gene
gff_exons <- gff %>% 
    filter(type != 'gene') %>%
    group_by(ID) %>%
    mutate(n_exon=sum(type == 'exon')) %>%
    filter(type == 'mRNA') %>% 
    mutate(mrna_len = end - start)

# Combine informations and reorder columns
gff_stats <- gff_gene_len %>%
    select(ID, gene_len) %>%
    inner_join(gff_exons, by='ID') %>%
    select(ID, chrom, start, end, type, strand, n_exon, gene_len, mrna_len, attributes)


### VISUALISE ###

# Get sample name from file name
clean_name <- strsplit(basename(in_gff), "\\.")[[1]][1]
# Label generator to get number of samples
n_obs <- function(x){return(data.frame(y=0, label=paste0("n=", length(x))))}

feature_len <- ggplot(data=gff %>%filter(type %in% c("CDS", "gene")), 
                      aes(x=type, y=log10(end - start))) +
    geom_violin(scale='width') +
    stat_summary(fun.data = n_obs, geom = "text") +
    theme_minimal() +
    xlab("") + 
    ylab("log10 length") + 
    ggtitle(clean_name) +
    scale_y_continuous(limits = c(0, 6))

q25 <- unname(quantile(gff_exons$n_exon, 0.25))
med <- median(gff_exons$n_exon)
q75 <- unname(quantile(gff_exons$n_exon, 0.75))
q_labs <- c(
    "q25"=paste0("25% quantile: ", q25), 
    "q75"=paste0("75% quantile: ", q75), 
    "med"=paste0("Median: ", med)
)
gff_exons <- gff_exons %>% 
  mutate(n_exon_quant = case_when(
    n_exon <= q25 ~ q_labs["q25"],
    n_exon >= q75 ~ q_labs["q75"],
    TRUE          ~ q_labs["med"]))
gff_exons <- gff_exons %>% 
    mutate(n_exon_quant = factor(n_exon_quant, ordered=T, 
                                 levels=unname(q_labs[c('q25', 'med', 'q75')])))

exon_per_gene <- ggplot(data=gff_exons, aes(x=n_exon, fill=n_exon_quant)) +
    scale_fill_manual(values=c("#6666ee", "#999999", "#ee6666")) +
    geom_histogram(binwidth=1) +
    geom_vline(xintercept=med) +
    theme_minimal() +
    ylab("Genes") +
    xlab("Exons per genes") +
    scale_x_continuous(limits = c(-0, 100))

### SAVE OUTPUT ###
svg(out_plot)
grid.arrange(feature_len, exon_per_gene)
dev.off()
write_tsv(gff_stats, out_tbl)

