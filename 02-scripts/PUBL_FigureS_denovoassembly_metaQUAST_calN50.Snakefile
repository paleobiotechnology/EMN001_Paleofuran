################################################################################
# Project: Natural products from the Palaeolithic
# Part: Publication
# Step: Prepare the Supplementary Figure evaluating the de novo assembly
#       performance using metaQUAST and calN50
#
# Dependent on:
#   - ASMB_denovo_assembly_binning.Snakefile
#   - PUBL_Dataset_S2.Snakefile
#
# Alex Huebner, 19/06/22
################################################################################

rule all:
    output:
        pdf = "06-figures_tables/FigureS_ASMB02.pdf",
        png = "06-figures_tables/FigureS_ASMB02.png"
    message: "Plot the results investigating de novo assembly performance using freeBayes"
    params:
        assembly = "06-figures_tables/Dataset_S2.xlsx"
    threads: 1
    script:
        "rscripts/PUBL_FigureS_denovoassembly_metaQUAST_calN50.R"
