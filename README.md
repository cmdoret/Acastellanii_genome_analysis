# _Acanthamoeba castellanii_ genome assembly and analysis

This repository contains the Assembly and comparative of two strains of Acanthamoeba castellanii.
The assembly process was done by combining Oxford Nanopore long reads with illumina shotgun libraries and Hi-C data. It is documented [in a separate file](doc/assembly_pipeline.md).

The comparative analysis can be run using the Snakefile.

Requirements:
* python3
    + snakemake
    + Biopython
    + BCBio
    + numpy
    + pandas
* R 3.5.x
    + tidyverse
    + GenomicRanges
    + ggbio
    + gridExtra
