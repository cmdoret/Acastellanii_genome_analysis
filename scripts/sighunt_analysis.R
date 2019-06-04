# Running HGT detection based on kmer frequencies using sighunt.
# Returns windows with associated discrete interval accumulative score  (DIAS)
# representing deviations of 4-mer frequencies from local densities.
# cmdoret, 20190529

library(sighunt)
library(tidyverse)
### Parse args
args <- commandArgs(trailingOnly=T)

genome <- args[1]
out_candidates <- args[2]

# function that executes the sighunt analysis on one contig
get_candidate <- function(sequence, win=5000, step=1000){
  # Only including contigs larger than two windows
  if(nchar(sequence) > 2*win){
    signature <- get_signature(sequence, window = win, step = step)
    dias <- global_density(signature)
  }
  else{
    dias <- c()
  }
  #dias <- dias[dias > threshold]
  return(dias)
}

# Read input genome
myseq <- read_fasta(genome)
candidates <- unlist(lapply(myseq, get_candidate))

# Parse the candidate names to obtain chromosome and positions.
# This will break if chromosome names contain . or -
cand_tbl <- tibble(name=names(candidates), dias=candidates) %>%
    mutate(pos=purrr::map(name, function(x)unlist(str_split(x, "[\\.-]")))) %>% 
    mutate(chrom=map_chr(pos, function(x) x[1]),
           start=map_int(pos, function(x) as.integer(x[2])), 
           end  =map_int(pos, function(x) as.integer(x[3]))) %>%
    select(chrom, start, end, dias, -pos, -name)

# Aggregate windows in 50kb bins
agg_tbl <- cand_tbl %>%
  group_by(chrom, window=round(start / 50000)) %>%
  summarise(start=min(start), end=max(end), agg_dias = mean(dias)) %>%
  select(chrom, start, end, agg_dias)

write_tsv(cand_tbl, out_candidates, col_names=T)
write_tsv(agg_tbl, paste0(out_candidates, ".agg"), col_names=T)
