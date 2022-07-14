################################################################################
# Project: Natural products from the Palaeolithic
# Part: Publication
# Step: Prepare the Supplementary Figure visualising the taxonomic assignment
#       of the MAGs based on the GTDB classification
#
# Dependent on:
#   - PUBL_Dataset_S3.Snakefile
#
# Alex Huebner, 14/07/22
################################################################################

rule all:
    output:
        pdf = "06-figures_tables/FigureS_ASMB06.pdf",
        png = "06-figures_tables/FigureS_ASMB06.png"
    message: "Plot the summary of taxonomic classification based on the GTDB"
    params:
        dataset_s3 = "06-figures_tables/Dataset_S3.xlsx",
    threads: 1
    script:
        "rscripts/PUBL_FigureS_MAG_taxonomicclassification.R"
