################################################################################
# Project: Natural products from the Palaeolithic
# Part: Publication
# Step: Prepare the Supplementary Figure visualisng the quality of the metaWRAP
#       refined MAGs based on the estimates of checkM, BUSCO, and GUNC
#
# Dependent on:
#   - PUBL_Dataset_S3.Snakefile
#
# Alex Huebner, 05/07/22
################################################################################

rule all:
    output:
        pdf = "06-figures_tables/FigureS_ASMB04.pdf",
        png = "06-figures_tables/FigureS_ASMB04.png"
    message: "Plot the summary of MAG quality estimated by checkM, BUSCO, and GUNC"
    params:
        dataset_s3 = "06-figures_tables/Dataset_S3.xlsx",
    threads: 1
    script:
        "rscripts/PUBL_FigureS_MAGquality_prefilter.R"
