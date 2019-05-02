library(tidyverse)
library(GenomicRanges)
library(ggbio)
args <- commandArgs(trailingOnly=TRUE)

# input
blast_file <- args[1]
names_file <- args[2]
# output
figure <- args[3]

# Loading data
blast_fmt6_cols <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
virus_names <- read_tsv(names_file, col_names=FALSE)
blast_output <- read_tsv(blast_file, col_names=FALSE)

# Giving meaningful names to BLAST columns, settin start to lowest genomic
# pos and end to highest (regardless of strand) and subsetting cols
blast_simple <- blast_output %>% 
  setNames(blast_fmt6_cols) %>%
  mutate(start=pmin(sstart, send), end=pmax(sstart, send)) %>%
  select(virus_id=qseqid, chrom=sseqid, start, end, evalue)


# Count occurences of virus types in BLAST hits
blast_simple %>% 
  group_by(virus_id) %>% 
  count() %>% 
  left_join(virus_names, by=c("virus_id"="X2")) %>%
  ggplot(aes(x="", y=n, fill=X1))+
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0)+
  theme_minimal() + 
  ggtitle("A.castellanii Old virus hits")

# Intersecting all genomic ranges to get non-overlapping virus sequences
blast_interval <- makeGRangesFromDataFrame(blast_simple)
no_overlaps <- GenomicRanges::reduce(blast_interval)
ggbio(no_overlaps)
gr <- GRanges(seqnames = rep('chr1', 6), 
              IRanges(start = c(1, 500, 1000, 2500, 10000, 20000), 
                      end = c(499, 999, 2499, 9999, 19999, 30000) ), 
              strand = rep('*', 6),
              name = sample(c('A', 'B'), size = 6, replace =T ) )

data <- GRanges(seqnames = rep('chr1', 100), 
                IRanges(start = sample(runif(100, min = 0, max = max(end(gr) ) ) ), 
                        width = 50), 
                d= runif(100, min = 0, max = 30) )

a <- ggbio() + circle(no_overlaps,  geom = 'rect', aes(fill = seqnames), space.skip = 0.01) 
a + circle(data , geom  = 'bar', aes(y = d) )
