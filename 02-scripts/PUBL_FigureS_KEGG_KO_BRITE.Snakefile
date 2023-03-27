################################################################################
# Project: Natural products from the Palaeolithic
# Part: Publication
# Step: Prepare heatmaps showing the frequency of specific KEGG orthologs
#       clustered at the KEGG BRITE 3rd-level hierarchies and more specifically
#       of the photosynthesis pathway
#
# Dependent on:
#   - FUNC_pangenome_Climicola.Snakefile
#   - PUBL_Dataset_S8.Snakefile
#
# Alex Huebner, 27/03/23
################################################################################

rule all:
    output:
        pdf = "06-figures_tables/FigureS_FUNC001.pdf",
        png = "06-figures_tables/FigureS_FUNC001.png"
    message: "Plot the frequency of KEGG orthologs specific for either C. limicola or the ancient Chlorobium MAGs and a heatmap for the presence/absence of photosynthesis genes"
    params:
        dataset_s8 = "06-figures_tables/Dataset_S8.xlsx",
        ko_photo = "05-results/FUNC_ko00194_kegg_pathway.txt"
    threads: 1
    script:
        "rscripts/PUBL_FigureS_KEGG_KO_BRITE.R"
