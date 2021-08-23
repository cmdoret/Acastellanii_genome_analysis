# This script generates coverage plot to inform us about ploidy
# It takes the bedgraph output from tinycov (github.com/cmdoret/tinycov)
# which contains mean coverage in sliding windows.
# cmdoret, 20210820

library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)

args <- commandArgs(trailingOnly = T)
bedgraph_in <- args[1]
plot_out <- args[2]

bg <- read_tsv(
  bedgraph_in,
	col_names = c('chrom', 'start', 'end', 'cov')
)

win_size <- 100000
win_step <- 10000

# Only keep chromosomes longer than win_size
# Trim windows in the top 0.01% coverage (e.g. rDNA)
# Compute coverage relative to genome average
df <- bg %>%
  group_by(chrom) %>%
  filter(any(max(end > win_size))) %>%
  ungroup() %>%
  filter(cov < quantile(cov, .99)) %>%
  mutate(cov_ratio = cov/median(cov)) %>%
  mutate(chrom=chrom)

# Order chromosomes by size
chrom_levels <- df %>%
  group_by(chrom) %>%
  filter(end == max(end)) %>%
  ungroup() %>%
  arrange(-end) %>%
  pull(chrom)

df$chrom = factor(df$chrom, ordered = T, levels = chrom_levels)
  
#ggplot(data=df, aes(x=start, y=log10(cov))) + geom_point() + facet_wrap(~chrom)
boxes <- ggplot(
    data=df,
    aes(x=chrom, y=cov_ratio, col=chrom)
  ) +
  geom_jitter(alpha=0.3) +
  geom_boxplot(col='black', fill=NA) +
  ylim(0.25, 2.0) +
  ylab("Coverage relative to genome median") + 
  xlab("Chromosome") +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90))

chrom_avg <- df %>%
  group_by(chrom) %>%
  summarise(avg = median(cov_ratio))


rug <- ggplot(data=chrom_avg, aes(x=avg)) +
  geom_rug() +
  xlim(0.25, 2.0) +
  coord_flip() +
  theme_minimal()

svg(plot_out)
boxes + rug + plot_layout(widths = c(5, 1))
dev.off()