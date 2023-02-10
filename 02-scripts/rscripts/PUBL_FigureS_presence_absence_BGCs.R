###############
# Title: Presence and absence of BGC in modern ref genomes and ancient MAGs
###############
# Author: Anan Ibrahim = darcy220606
###############
# Date: July 2022
###############
# Libraries required ####
library(tidyverse)
library(pheatmap)
library('readr') # if not installed already
library(dplyr)
library(ggplot2)
######################################################
# Set-WD ####
setwd("/media/aibrahim/INTENSO/Deep-Evo/final/antismash_mags")
######################################################
# Grab the presence absence of BGCs across all drep genomes ####
table <- read_table("antismash_drepgenomes_table2.txt", col_names = TRUE)

# Order of the tree (drep chlorobiales 99)
samples <- c(
  "GCA_013824245.1",  "GCA_013824085.1",
  "GCA_014380515.1",  "GCA_002763895.1",
  "GCA_002794105.1",  "GCF_000020525.1",
  "GCA_013334725.1",
  "GCA_013334605.1",
  "GCF_002240205.1",
  "GCA_016776555.1",
  "GCA_010912745.1",
  "GCA_004028745.1",
  "GCA_015482705.2",
  "GCF_003182595.1",
  "GCF_001687065.1",
  "GCF_003344365.2",
  "GCF_002113825.1",
  "GCA_010912735.1",
  "GCF_000020625.1",
  "GCF_003344375.1",
  "GCA_903929135.1",
  "GCA_013334925.1",
  "GCA_013335825.1",
  "GCA_013335115.1",
  "GCA_013334865.1",
  "GCA_013334965.1",
  "GCA_013334525.1",
  "GCA_903878795.1",
  "GCA_013334975.1",
  "GCA_016931495.1",
  "GCF_000020505.1",
  "GCA_011372505.1",
  "GCA_011372555.1",
  "GCA_001602925.1",
  "GCF_000006985.1",
  "GCF_001747405.1",
  "GCF_004118635.1",
  "GCF_006265165.1",
  "GCA_013334645.1",
  "GCA_903901395.1",
  "GCA_013334855.1",
  "GCA_013334655.1",
  "GCA_013336025.1",
  "GCA_019163275.1",
  "GCF_000015125.1",
  "GCA_019162945.1",
  "GCA_016927285.1",
  "GCA_013334435.1",
  "GCF_000020465.1",
  "GCF_001509575.1",
  "PLV001_002",
  "PES001_018",
  "GOY005_001",
  "GOY006_RA",
  "TAF017_RA",
  "PLV001_001",
  "EMN001_021",
  "RIG001_014",
  "GCA_013203905.1",
  "GCA_003517435.1",
  "GCF_008634005.1",
  "GCF_004332115.1",
  "GCA_001622165.1",
  "GCF_000012485.1",
  "GCA_005862285.1",
  "GCA_013334225.1",
  "GCA_903928225.1",
  "GCA_013334565.1",
  "GCA_013335525.1",
  "GCA_013335405.1",
  "GCA_903914905.1",
  "GCA_013334485.1",
  "GCF_011601315.1",
  "GCF_001975465.1",
  "GCF_000168715.1",
  "GCA_903839405.1",
  "GCA_005843805.1",
  "GCA_903893775.1",
  "GCA_903886455.1",
  "GCA_012972835.1",
  "GCA_903994805.1",
  "GCA_012972825.1",
  "GCA_903836935.1",
  "GCA_903940125.1",
  "GCA_903904495.1",
  "GCA_903876135.1",
  "GCA_903895195.1",
  "GCA_903994515.1",
  "GCA_903908555.1",
  "GCA_903851385.1",
  "GCA_903936995.1",
  "GCA_903829495.1",
  "GCA_903994485.1",
  "GCA_903920375.1",
  "GCA_903903455.1",
  "GCA_005843825.1",
  "GCA_903830695.1",
  "GCA_903890365.1",
  "GCA_903934085.1",
  "GCA_013334755.1",
  "GCF_000020645.1",
  "GCA_903830225.1",
  "GCA_903865705.1",
  "GCA_020035305.1",
  "GCA_903942355.1",
  "GCA_903957875.1",
  "GCA_903836805.1",
  "GCA_903822595.1",
  "GCA_903876775.1",
  "GCA_903939765.1",
  "GCA_903838025.1",
  "GCA_903875625.1",
  "GCA_903904565.1",
  "GCA_903916475.1",
  "GCA_903851535.1",
  "GCA_903840935.1",
  "GCA_903835515.1",
  "GCA_903904115.1",
  "GCA_903913415.1",
  "GCA_903994415.1",
  "GCA_903874585.1",
  "GCA_903859055.1",
  "GCA_903994295.1",
  "GCA_005843815.1",
  "GCA_903866255.1",
  "GCA_903994365.1",
  "GCA_019163115.1",
  "GCA_903994465.1",
  "GCA_903884705.1",
  "GCA_005862225.1",
  "GCA_903838095.1",
  "GCA_903932555.1",
  "GCA_903843305.1",
  "GCA_903842325.1",
  "GCA_903940785.1",
  "GCA_903994455.1",
  "GCA_903926375.1",
  "GCA_903870225.1",
  "GCA_903994655.1",
  "GCA_903836465.1",
  "GCA_903836115.1",
  "GCA_903898135.1",
  "GCA_903847795.1")

pdf("BGCS_in_drep_mags.pdf")
ggplot(table, aes(x = BGC_core_genes, y = factor(Sample_name, levels = samples))) +
  geom_point(size=1, shape=20, color="black", stroke=0.2) +
  scale_y_discrete(limits=rev)+
  theme(
    aspect.ratio = 2,
    axis.text.x = element_text(angle = 75, vjust = 1, hjust = 1),
    axis.text.y = element_text(size=3.5),
    axis.line = element_line(size = 0.1, colour = "black", linetype=1)) 
# +  theme(panel.grid.minor = element_line(colour="blue", size=0.5)) 
dev.off()

jpeg("BGCS_in_drep_mags.jpeg")
ggplot(table, aes(x = BGC_core_genes, y = factor(Sample_name, levels = samples))) +
  geom_point(size=1, shape=20, color="black", stroke=0.2) +
  scale_y_discrete(limits=rev)+
  theme(
    aspect.ratio = 2,
    axis.text.x = element_text(angle = 75, vjust = 1, hjust = 1),
    axis.text.y = element_text(size=3.5),
    axis.line = element_line(size = 0.1, colour = "black", linetype=1)) 
# +  theme(panel.grid.minor = element_line(colour="blue", size=0.5)) 
dev.off()
######################################################

