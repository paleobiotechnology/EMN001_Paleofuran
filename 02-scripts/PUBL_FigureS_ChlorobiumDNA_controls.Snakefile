################################################################################
# Project: Natural products from the Palaeolithic
# Part: Publication
# Step: Prepare the Supplementary Figure visualisng the relative abundance of
#       Chlorobium DNA in the dental calculus samples, environmental samples
#       from El Miron, and the laboratory negative controls
#
# Dependent on:
#   - QUAL_refalignment_Chlorobiaceae.Snakefile
#   - PUBL_Dataset_S7.Snakefile
#
# Alex Huebner, 30/08/22
################################################################################

rule all:
    output:
        pdf = "06-figures_tables/FigureS_QUAL01.pdf",
        png = "06-figures_tables/FigureS_QUAL01.png"
    message: "Plot the relative abundance of Chlorobium DNA in different sample types"
    params:
        dentalcalculus = "05-results/QUAL_dentalcalculus_Chlorobiaceae_refalignment.tsv",
        dataset_s7 = "06-figures_tables/Dataset_S7.xlsx",
    threads: 1
    script:
        "rscripts/PUBL_FigureS_ChlorobiumDNA_controls.R"
