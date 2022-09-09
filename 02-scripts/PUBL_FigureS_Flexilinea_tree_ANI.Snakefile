################################################################################
# Project: Natural products from the Palaeolithic
# Part: Publication
# Step: Prepare the Supplementary Figure visualising tree of the genus
#       Flexilinea and the pairwise ANI analysis
#
# Disclaimer: The phylogenetic tree produced by RAxML was visualised using iTOL
#             and later combined with the pairwise ANI plot.
#
# Dependent on:
#   - PUBL_Dataset_S8.Snakefile
#
# Alex Huebner, 09/09/22
################################################################################

rule all:
    output:
        pdf = "06-figures_tables/FigureS_PHYL10_pANI.pdf"
    message: "Plot the pairwise ANI analysis "
    params:
        dataset_s8 = "06-figures_tables/Dataset_S8.xlsx",
    threads: 1
    script:
        "rscripts/PUBL_FigureS_Flexilinea_tree_ANI.R"
