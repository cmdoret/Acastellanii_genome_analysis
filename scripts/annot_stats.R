# Computing general statistics on input genome annotations.
# Usage: annots_stats input.gff output.pdf
# cmdoret, 20190509

library(tidyverse)
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
gff %>% 
    mutate(ID=str_split(str_extract(attributes, 'ID=[^;]*'), '=', n=2)[[1]][2]) %>%
    group_by(ID) %>%
    mutate(N_exons=tally(type))
# Annotated product (Y/N)

# N GO term

### VISUALISE ###

clean_name <- strsplit(basename(in_gff), "\\.")[[1]]
# Length distribution of different types
n_obs <- function(x){return(data.frame(y=0, label=paste0("n=", length(x))))}

feature_len <- ggplot(data=gff, aes(x=type, y=log10(end - start))) + 
    geom_violin(scale='width') + 
    stat_summary(fun.data = n_obs, geom = "text") +
    theme_bw() +
    xlab("") + 
    ylab("log10 length") + 
    ggtitle(clean_name)

exon_per_gene <- ggplot(data=gff, aes()) +
    geom_point()

### SAVE OUTPUT ###
pdf(out_plot)


dev.off()


