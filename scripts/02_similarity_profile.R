# Visualise blast similarity against bacteria from target and background gene families.

library(tidyverse)

args <- commandArgs(trailingOnly=T)

target read_tsv(args[1])
background <- read_tsv(args[2])
out_figure <- args[3]

