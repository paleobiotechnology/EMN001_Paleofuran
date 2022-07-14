################################################################################
# Project: Natural products from the Palaeolithic
# Part: Publication
# Step: Prepare the Supplementary Figure visualising the ancient DNA damage at
#       the 5' end of the non-UDG treated reads against the MAGs assigned to
#       the genera Chlorobium and Flexilinea
#
# Dependent on:
#   - PUBL_Dataset_S3.Snakefile
#
# Alex Huebner, 14/07/22
################################################################################

rule all:
    output:
        pdf = "06-figures_tables/FigureS_ASMB07.pdf",
        png = "06-figures_tables/FigureS_ASMB07.png"
    message: "Plot the ancient DNA at the 5' end of reads"
    params:
        dataset_s3 = "06-figures_tables/Dataset_S3.xlsx",
    threads: 1
    script:
        "rscripts/PUBL_FigureS_aDNA_Chlorobium_Flexilinea.R"
