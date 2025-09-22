# CRDG interview (Alexander Blakes)
Pre-interview task for the CRDG postdoc.

## Analysis pipeline
The analysis is run as a Snakemake pipeline, driven by `crdg_interview/Snakefile`. The pipeline has a simple linear logic, as defined in the Snakefile. 

The "worker" .smk files used to run each step of the analysis can be found in `crdg_interview/workflow/`.

The scripts which are used to run the analysis can be found under `crdg_interview/experiments/`.

Some utility functions and scripts are stored under `crdg_interview/src/`.

The `crdg_interview` directory is used as the working directory for this analysis.

## Virtual environments
Package management is with Mamba v2.0.5. The virtual environments used in this analysis are documented in `crdg_interview/env/`.

The Snakemake pipeline depends on a local install of these environments, so be sure to install them with `conda install -f <yaml_file_name>` before running the pipeline.

## gnomAD v3.1.2 download
The gnomAD v3.1.2 download step is not included in the Snakemake pipeline. But the scripts with which it was accomplished are saved in `crdg_interview/experiments/download_gnomad_v3_genomes/`.

