# analysing distribution of gene familises among species of amoeba and A. castellanii
# cmdoret, 20191014

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggupset)

pres_mat <- read_tsv('data/out/specific_genes/presence_matrix.tsv')
pres_mat$outgroup <-  apply(X=pres_mat[, -c(1:3)], FUN=any, MARGIN=1)
pres_mat <- pres_mat[, c("Orthogroup", "C3_proteome", "Neff_proteome", "outgroup")]

pres_tidy <- pres_mat %>%
  as_tibble() %>%
  gather(Organism, Member, -Orthogroup) %>%
  filter(Member) %>%
  select(- Member)

pres_tidy %>%
  group_by(Orthogroup) %>%
  summarize(Organisms = list(Organism))

pres_tidy %>%
  group_by(Orthogroup) %>%
  summarize(Organisms = list(Organism)) %>%
  ggplot(aes(x = Organisms)) +
  geom_bar() +
  scale_x_upset(order_by='freq', n_intersections = 24)
