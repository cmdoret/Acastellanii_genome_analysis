# Plot density of foreign origin genes along genome
# cmdoret, 20190521

library(tidyverse)
library(circlize)

args <- commandArgs(trailingOnly=TRUE)

### LOAD DATA ###
annot <- args[1]
out_fig <- args[2]

### FORMAT DATA ###
