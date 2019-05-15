# Computing general statistics on input genome annotations.
# Usage: annots_stats input.gff output.pdf
# cmdoret, 20190509

library(tidyverse)
library(gridExtra)
args <- commandArgs(trailingOnly=TRUE)

### GET ARGS ###
in_gff <- args[1]
out_plot <- args[2]

### LOAD DATA ###
gff <- read_tsv(in_gff, comment="#", col_names=F)
colnames(gff) <- c("chrom", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
gff <- gff %>% filter(type %in% c("mRNA", "CDS", "gene", "exon", "tRNA"))

### COMPUTE STATS ###

# N exon / gene
gff_exons <- gff %>% 
    filter(type != 'gene') %>%
    mutate(attributes=gsub("Parent=", "ID=", attributes)) %>%
    mutate(ID=str_extract(attributes, 'ID=[^;]*')) %>%
    mutate(ID=gsub("=[a-zA-Z]{1,}:", "=", ID)) %>%
    mutate(ID=sapply(ID, function(x){str_split(x,regex('[=-]'))[[1]][2]})) %>%
    group_by(ID) %>%
    mutate(n_exon=sum(type == 'exon')) %>%
    filter(type == 'mRNA')


### VISUALISE ###

# Get sample name from file name
clean_name <- strsplit(basename(in_gff), "\\.")[[1]][1]
# Label generator to get number of samples
n_obs <- function(x){return(data.frame(y=0, label=paste0("n=", length(x))))}

feature_len <- ggplot(data=gff, aes(x=type, y=log10(end - start))) + 
    geom_violin(scale='width') + 
    stat_summary(fun.data = n_obs, geom = "text") +
    theme_minimal() +
    xlab("") + 
    ylab("log10 length") + 
    ggtitle(clean_name)

q25 <- unname(quantile(gff_exons$n_exon, 0.25))
med <- median(gff_exons$n_exon)
q75 <- unname(quantile(gff_exons$n_exon, 0.75))
gff_exons <- gff_exons %>% mutate(n_exon_quant = case_when(
    n_exon <= q25 ~ paste0("25% quantile: ", q25),
    n_exon >= q75 ~ paste0("75% quantile: ", q75),
    TRUE          ~ paste0("Median: ", med)
))

exon_per_gene <- ggplot(data=gff_exons, aes(x=n_exon, fill=n_exon_quant)) +
    scale_fill_manual(values=c("#6666ee", "#999999", "#ee6666")) +
    geom_histogram(binwidth=1) +
    geom_vline(xintercept=med) + 
    theme_minimal() + 
    ylab("Genes") +
    xlab("Exons per genes")


### SAVE OUTPUT ###
pdf(out_plot)
grid.arrange(feature_len, exon_per_gene)
dev.off()


