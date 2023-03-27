###############
# Title: Heatmap illustrating teh fastANI values of ancient and modern genomes
###############
# Author: Anan Ibrahim = darcy220606
###############
# Date: July 2022
###############
# Libraries required

library('splitstackshape')
library('pheatmap')
library('dplyr') # if not installed already
library('readr') # if not installed already
library('RColorBrewer')
######################################################
# Set-WD ####
setwd("/media/aibrahim/INTENSO/Deep-Evo/final/Tree-chlorobiales-genomes")
######################################################
# Pruned tree version2 - Fraglen 1000 ####
# pims <- read.table("ANI_OUT_MATRIX_fg1000.txt", header=TRUE, sep="\t",row.names = 1)
pims <- read.table("ANI_output_fg1000.txt.matrix", header=TRUE, sep="\t",row.names = 1)
pims <- read.table("ANI_OUT_MATRIX_fg1000_filter.txt", header=TRUE, sep="\t",row.names = 1)

#Version#2
samples <- c("GCF_000020505.1.fna","GCA_011372555.1.fna","GCA_011372505.1.fna",
             "GCF_000006985.1.fna","GCF_006265165.1.fna","GCF_001747405.1.fna","GCF_000020525.1.fna",
             "GCA_002794105.1.fna","GCF_003182595.1.fna","GCA_014384625.1.fna","GCA_016776555.1.fna",
             "GCF_006265245.1.fna","GCF_000020625.1.fna","GCA_011054385.1.fna","GCF_013357945.1.fna",
             "GCF_015233875.1.fna","GCF_015234275.1.fna","GCF_000168715.1.fna","GCF_011601315.1.fna",
             "GCA_019163155.1.fna","GCF_000020645.1.fna","GCF_000012485.1.fna","GCA_001622165.1.fna",
             "GCF_003970105.1.fna","GCF_008634005.1.fna","GCF_009807145.1.fna","GCF_008586845.1.fna",
             "GCF_003968655.1.fna","GCF_000015125.1.fna","GCA_019163275.1.fna","GCF_001509575.1.fna",
             "GCA_013335335.1.fna","GCA_013335765.1.fna","GCF_000020465.1.fna","PLV001_002.fna",
             "TAF017_RA.fna","PLV001_001.fna","EMN001_021.fna","RIG001_014.fna",
             "PES001_018.fna","GOY005_001.fna","GOY006_RA.fna")

#reoder columns according to vector
pims1 <- pims[, samples] #reorder columns according to order of the samples in the pylogentic tree
#reoder rows according to vector
pims1 <- pims1[match(samples, rownames(pims1)), ]  

indx <- sapply(pims1, is.factor)
pims1[indx] <- lapply(pims1[indx], function(x) as.numeric(as.character(x)))

col <- c("#F497B6","#E8749A","#D12E64", "#87153B","#5E0C27", "#97E4F4","#74C9E8", "#2D92D1","#174E86","#0C285E")

pdf("ANI_output_fg1000.pdf")
pheatmap(pims1, display_numbers = F,
         main = "ANI_output_fg1000.txt.matrix",color=col,
         show_rownames = TRUE,show_colnames = TRUE,fontsize  = 5,
         cluster_cols = F,cluster_rows=F)
dev.off()

jpeg("ANI_output_fg1000.jpg")
pheatmap(pims, display_numbers = F,
         main = "ANI_output_fg1000.txt.matrix",color=col,
         show_rownames = TRUE,show_colnames = TRUE,fontsize  = 5,
         cluster_cols = F,cluster_rows=F)
dev.off()

######################################################
## Pruned tree version2 - Fraglen 3000 ####

# pims <- read.table("ANI_OUT_MATRIX_fg1000.txt", header=TRUE, sep="\t",row.names = 1)
pims <- read.table("ANI_output_fg3000.txt.matrix", header=TRUE, sep="\t",row.names = 1)
pims <- read.table("ANI_OUT_MATRIX_fg3000_filter.txt", header=TRUE, sep="\t",row.names = 1)

#Version#2
samples <- c("GCF_000020505.1.fna","GCA_011372555.1.fna","GCA_011372505.1.fna",
             "GCF_000006985.1.fna","GCF_006265165.1.fna","GCF_001747405.1.fna","GCF_000020525.1.fna",
             "GCA_002794105.1.fna","GCF_003182595.1.fna","GCA_014384625.1.fna","GCA_016776555.1.fna",
             "GCF_006265245.1.fna","GCF_000020625.1.fna","GCA_011054385.1.fna","GCF_013357945.1.fna",
             "GCF_015233875.1.fna","GCF_015234275.1.fna","GCF_000168715.1.fna","GCF_011601315.1.fna",
             "GCA_019163155.1.fna","GCF_000020645.1.fna","GCF_000012485.1.fna","GCA_001622165.1.fna",
             "GCF_003970105.1.fna","GCF_008634005.1.fna","GCF_009807145.1.fna","GCF_008586845.1.fna",
             "GCF_003968655.1.fna","GCF_000015125.1.fna","GCA_019163275.1.fna","GCF_001509575.1.fna",
             "GCA_013335335.1.fna","GCA_013335765.1.fna","GCF_000020465.1.fna","PLV001_002.fna",
             "TAF017_RA.fna","PLV001_001.fna","EMN001_021.fna","RIG001_014.fna",
             "PES001_018.fna","GOY005_001.fna","GOY006_RA.fna")

#reoder columns according to vector
pims1 <- pims[, samples] #reorder columns according to order of the samples in the pylogentic tree
#reoder rows according to vector
pims1 <- pims1[match(samples, rownames(pims1)), ]  

indx <- sapply(pims1, is.factor)
pims1[indx] <- lapply(pims1[indx], function(x) as.numeric(as.character(x)))

col <- c("#F497B6","#E8749A","#D12E64", "#87153B","#5E0C27", "#97E4F4","#74C9E8", "#2D92D1","#174E86","#0C285E")

pdf("ANI_output_fg3000.pdf")
pheatmap(pims1, display_numbers = F,
         main = "ANI_output_fg3000.txt.matrix",color=col,
         show_rownames = TRUE,show_colnames = TRUE,fontsize  = 5,
         cluster_cols = F,cluster_rows=F)
dev.off()

jpeg("ANI_output_fg3000.jpg")
pheatmap(pims, display_numbers = F,
         main = "ANI_output_fg3000.txt.matrix",color=col,
         show_rownames = TRUE,show_colnames = TRUE,fontsize  = 5,
         cluster_cols = F,cluster_rows=F)
dev.off()
######################################################