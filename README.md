# CRDG interview (Alexander Blakes)
A public repository containing the analysis pipeline for the CRDG postdoc pre-interview task.

This analysis quantifies differences in variant density and allele frequency between snRNA genes and pseudogenes, and explores the basis for heterozygote excess in certain snRNAs.

## Analysis pipeline
This analysis is run as a Snakemake pipeline, driven by `crdg_interview/Snakefile`. The pipeline has a simple linear logic, as defined in the Snakefile. The order in which each analysis step is run is defined in the default target rule, `all`, at the bottom of the Snakefile. Running this pipeline will reproduce all the statistics and plots used in the analysis.

## Directory structure
The "worker" .smk files used to run each step of the analysis can be found in `crdg_interview/workflow/`.

The scripts which are executed by Snakemake can be found under `crdg_interview/experiments/`.

Some utility functions and scripts are stored under `crdg_interview/src/`.

The `crdg_interview` directory is used as the working directory for this analysis.

The layout of the `data/` directories used in this analysis is recorded through placeholder `.gitignore` files.

## Virtual environments
Mamba v2.0.5 was used to install and manage virtual environments. The virtual environments used in this analysis are documented in `.yml` files under `crdg_interview/env/`.

The Snakemake pipeline depends on a local install of these environments, so be sure to install them with `conda install -f <yaml_file_name>` before running the pipeline. Alternatively, the `.smk` scripts can readily be updated to accept the `.yml` definition files.

## gnomAD v3.1.2 files
The gnomAD v3.1.2 download step is not included in the Snakemake pipeline. But the scripts with which it was accomplished are available in `crdg_interview/experiments/download_gnomad_v3_genomes/`.

The absolute path to the directory containing the gnomAD data is hard-coded in `experiments/extract_variants/extract_gnomad_variants_in_snrnas.sh`, and should be changed before running the pipeline.

