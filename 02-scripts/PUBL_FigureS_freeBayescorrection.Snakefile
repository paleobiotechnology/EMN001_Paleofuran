################################################################################
# Project: Natural products from the Palaeolithic
# Part: Publication
# Step: Prepare the Supplementary Figure investigating the contig correction
#       using freeBayes
#
# Dependent on:
#   - QUAL_freeBayes_nomismatches_MEGAHIT.Snakefile
#
# Alex Huebner, 17/06/22
################################################################################

import os

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

rule all:
    output:
        pdf = "06-figures_tables/FigureS_ASMB01.pdf",
        png = "06-figures_tables/FigureS_ASMB01.png"
    message: "Plot the results investigating the contig correction using freeBayes"
    params:
        maf = "05-results/QUAL_freeBayes_distributionMAF.tsv",
        subst = "05-results/QUAL_freeBayes_nSubsts.tsv"
    threads: 1
    script:
        "rscripts/PUBL_FigureS_freeBayescorrection.R"
