################################################################################
# Project: Natural products from the Palaeolithic
# Part: Publication
# Step: Prepare the Supplementary Figure visualisng the pyDamage results on
#       contig level.
#
# Dependent on:
#   - PUBL_Dataset_S2.Snakefile
#
# Alex Huebner, 19/06/22
################################################################################

rule all:
    output:
        pdf = "06-figures_tables/FigureS_ASMB03.pdf",
        png = "06-figures_tables/FigureS_ASMB03.png"
    message: "Plot the summary of the pyDamage results on contig level"
    params:
        dataset_s1 = "06-figures_tables/Dataset_S1.xlsx",
        dataset_s2 = "06-figures_tables/Dataset_S2.xlsx"
    threads: 1
    script:
        "rscripts/PUBL_FigureS_pyDamage_contigs.R"
