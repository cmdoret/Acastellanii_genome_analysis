# _Acanthamoeba castellanii_ genome assembly and analysis

This repository contains the scripts and documentation for the assembly and comparative genomics of two strains of Acanthamoeba castellanii (Neff and C3). The assembly process was done by combining Oxford Nanopore long reads with illumina shotgun libraries and Hi-C data. The assembly and annotation pipelines are documented in separate files:

Currently, the genome of the C3 strain is still at the polishing stage it should should soon be finished to allow comparative analyses.

## [Assembly pipeline](doc/assembly_pipeline.md)
## [Annotation pipeline](doc/annotation_pipeline.md)

The comparative analysis can be run using snakemake (data not available yet).

Requirements:
* python3:
    + snakemake
    + Biopython
    + BCBio
    + numpy
    + pandas
    + pyomadb
    + requests
* R 3.6.x:
    + tidyverse
    + GenomicRanges
    + ggbio
    + gridExtra
    + [sighunt](https://github.com/KamilSJaron/sighunt)
    + topGO
* Softwares:
    + MCScanX
    + circos
