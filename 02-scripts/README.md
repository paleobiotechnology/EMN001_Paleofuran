## Overview of the folder `02-scripts`

This folder contains all the scripts and Jupyter notebooks to prepare and analyse the data.

1. `PREP`: preparing of the data

  - `PREP_download_from_ENA.Snakefile`: download the sequencing data from ENA with the help of the
    Python program [ffq](https://github.com/pachterlab/ffq)
  - `PREP_preprocessing_sequencing_data.Snakefile`: prepare input file for pre-processing the
    sequecing data with the Nextflow pipeline nf-core/eager and extract and merge all non-human
    reads per individual

X. `ENVS`: conda environments for use with Snakemake

  - `ENVS_bioawk.yaml`: [bioawk](https://github.com/lh3/bioawk)
  - `ENVS_samtools.yaml`: [samtools](https://github.com/samtools/samtools)
