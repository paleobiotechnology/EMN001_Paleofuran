## Overview of the folder `05-results`

This folder contains the results of the analyses that are used to produce the figures and tables
presented in the manuscript.

1. `PREP`: results related to the data preparation

  - `PREP_Nextflow_EAGER_noReads_per_sample.tsv`: the overview of the number of reads per individual
    and library type (all data and non-UDG data only) that were available for *de novo* assembly
  - `PREP_Nextflow_EAGER_noReads_Weyrich2017_Neanderthals.tsv`: the overview of the number of reads
    available for the four Neanderthal calculus samples published by Weyrich *et al.* (2017)
  - `PREP_Nextflow_EAGER_noReads_ElMiron_sediments.tsv`: the overview of the number of the number of
    reads available for the six El Miron sediment samples generated for this study

2. `ASMB`: results related to the *de novo* assembly and non-reference binning

  - `ASMB_assemblystats_calN50_metaQUAST.tsv`: the overview of the performance of the *de novo*
    assembly summarised using calN50, metaQUAST, Prokka and pyDamage

3. `QUAL`: results related to the data quality

  - `QUAL_fragmentlength_distribution.tsv`: the overview of the distribution of DNA molecule lengths
    across the sequencing data generated on 2x 75 bp paired-end Illumina sequencing runs
  - `QUAL_freeBayes_nSubsts.tsv`: the overview of the frequency of each allele difference between the
    allele exported by MEGAHIT from the assembly graph and the genotype call inferred by freeBayes
    was observed
  - `QUAL_freeBayes_distributionMAF.tsv`: the distribution of the allele frequencies for the allele
    that freeBayes inferred to be the genotype by 5% frequency bins
  - `QUAL_freeBayes_nContigs.tsv`: the overview of the number of contigs per sample for which a
    difference between the MEGAHIT and freeBayes sequence could be observed
