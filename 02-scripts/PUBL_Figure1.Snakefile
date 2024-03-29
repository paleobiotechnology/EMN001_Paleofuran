################################################################################
# Project: Natural products from the Palaeolithic
# Part: Publication
# Step: Prepare Figure 1 summarising the sample set, the de novo assembly, and
#       MAG binning
#
# Dependent on:
#   - ASMB_denovo_assembly_binning.Snakefile
#   - PUBL_Dataset_S1.Snakefile
#   - PUBL_Dataset_S2.Snakefile
#
# Alex Huebner, 22/06/22
################################################################################

rule all:
    output:
        pdf = "06-figures_tables/Figure1.pdf",
        png = "06-figures_tables/Figure1.png"
    message: "Prepare the draft of Figure 1"
    params: 
        dataset_s1 = "06-figures_tables/Dataset_S1.xlsx",
        dataset_s2 = "06-figures_tables/Dataset_S2.xlsx",
        dataset_s3 = "06-figures_tables/Dataset_S3.xlsx"
    script:
        "rscripts/PUBL_Figure1.R"
