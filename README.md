# _Acanthamoeba castellanii_ genome analysis

> TODO: Define conda environments for the remaining rules to minimize dependencies

## Description 
This repository contains scripts and documentation related to the analysis and comparison of the _Acanthamoeba castellanii_ genome from strains C3 and Neff.

## Installation

The pipeline is written using snakemake and manages dependencies using conda. Most of the pipeline steps are run inside self-contained conda environments, which are automatically built upon execution. There are two dependencies which cannot be installed using conda and need to be installed separately.

**Dependencies:**

* python3.7+
    + snakemake
    + pandas
    + numpy
* conda
* dnaglider
* MCScanX

Additionally, the input data must be provided (annotations, genomes, )
## Usage

 The analyses are separated into distinct workflows in the `rules` directory.
 The whole analysis pipeline can be run using snakemake as follows:

 ```snakemake --use-conda -j4```

## Structure

The master script `Snakefile` will call each workflow one after the other. Each workflow contains rules with input and output files, which execute code or external scripts. Each rule is executed in its own conda environment and will download its dependencies on the first execution. The overall workflow can be represented as a graph:

![pipeline graph](doc/dag.svg)

The `envs` directory contains conda environment build specifications for the different rules.

General parameters for the pipeline are stored in the `config.yaml` file and can be modified. The strains to analyze as well as the path to their sequence files are defined in `samples.tsv`. All external scripts executed by rules are stored in the `scripts` folder. Custom python utility libraries imported in the pipeline are stored in `src`.

The `doc` directory contains jupyter notebook with general analyses of the pipeline results.


