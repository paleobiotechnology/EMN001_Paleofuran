################################################################################
# Project: Natural products from the Palaeolithic
# Part: Publication
# Step: Prepare the Supplementary Figure visualisng the quality of the metaWRAP
#       refined MAGs based on the estimates of checkM and GUNC after the
#       automatic refinement
#
# Dependent on:
#   - PUBL_Dataset_S3.Snakefile
#
# Alex Huebner, 06/07/22
################################################################################

rule all:
    output:
        pdf = "06-figures_tables/FigureS_ASMB05.pdf",
        png = "06-figures_tables/FigureS_ASMB05.png"
    message: "Plot the summary of MAG quality estimated by checkM and GUNC"
    params:
        dataset_s3 = "06-figures_tables/Dataset_S3.xlsx",
    threads: 1
    script:
        "rscripts/PUBL_FigureS_MAGquality_postfilter.R"
