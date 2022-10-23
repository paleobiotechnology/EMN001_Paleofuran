###############
# Title: Bigscape analysis
###############
# Author: Anan Ibrahim = darcy220606
###############
# Date: July 2022
###############
# Libraries required
###############

# Change the directory accordingly 
# arg = commandArgs(trailingOnly = TRUE)

library('tidyverse')
library('pheatmap')
library('readr') # if not installed already
library('dplyr')
library('plyr')
library('data.table')
library('tidyr')
library('stringr')
library('RColorBrewer')
library('reshape2')
library('GGally')
library('network')
library('shadowtext')
library('ggplot2')
library(igraph)
library(ggraph)
library(graphlayouts)
library(ggforce)


######################################################
# Set-WD ####
setwd("/media/aibrahim/INTENSO/Deep-Evo/final/BIGSCAPE/BUTYR-1/network_files/2022-07-21_17-08-22_hybrids_auto_chlorobiales-drep-butyrolactone/Others")
######################################################
# HOW many genomes contains THE BGC? ####
# Filter for the Butyrolactone BGC
network <- read.table("Network_Annotations_Others.tsv", header=TRUE, sep="\t")

but <- network[grep("butyrolactone", network$Product.Prediction),] # grabs all rows that have the butyrolactone BGC
no.of.genomes.with.BGC <- nrow(but)
######################################################
# HOW many links exist between the genomes? ####
# Filter for the Butyrolactone BGC
network <- read.table("Others_c1.00.network", header=TRUE, sep="\t")
but <- network[grep("butyrolactone - butyrolactone", network$Combined.group),] # grabs all rows that have the butyrolactone BGC
no.of.links.between.all <- nrow(but)
write.csv(but, "Others_c1.00.network.butyrolactone.only.csv") # write the butyrolactone table only to generate the network on cytoscape
write.table(but, "Others_c1.00.network.butyrolactone.only.txt")

but <- filter(but,Raw.distance<0.7)
write.table(but, "Others_c1.00.network.butyrolactone.only.0.7.txt")

anc <- but[Reduce(`|`, lapply(but[-1], grepl, pattern="c00")),] #retain only rows that have contig name in any column
anc <- anc[ - grep("c00011_NZ_JAAO...region001", anc$Clustername.1),] # one modern group with c00 remove 
anc <- anc[ - grep("c00011_NZ_JAAO...region001", anc$Clustername.2),] # one modern group with c00 remove 
no.of.links.from.anc <- nrow(anc)
write.table(anc, "Others_c1.00.network.butyrolactone.ancient.only.txt")

# plv = anc[Reduce(`|`, lapply(anc[-1], grepl, pattern="c00055_PLV001")),] #retain only rows that have contig name in any column

######################################################
# How many families and clans in the butyrolactone? ####

# Grab only the butyrolactone links
network <- read.table("Others_c1.00.network", header=TRUE, sep="\t")
but <- network[grep("butyrolactone - butyrolactone", network$Combined.group),] # grabs all rows that have the butyrolactone BGC
names <- data.frame(unique(rbind(data.frame(V1 = unique(but$Clustername.1)),
                                 data.frame(V1 = unique(but$Clustername.2))
)))

cf <- read.table("Others_clans_1.00_1.00.tsv", header=FALSE, sep="\t",stringsAsFactors=FALSE)
cf <-right_join(cf, names, by = c("V1")) #merges by names that only have butyrolactone links
# write.table(cf, "Others_clans_1.00_1.00_butyrolactone.txt")
no.of.clans <- length(unique(cf$V2))
no.of.families <- length(unique(cf$V3))

count = data.frame(table(cf$V3))
colnames(count) = c("family number", "no. of genomes in family")
count = t(count)

output <- data.frame(rbind(no.of.genomes.with.BGC, 
                           no.of.links.between.all,
                           no.of.links.from.anc,
                           no.of.clans,
                           no.of.families,
                           count))

write.table(output, "bigscape_stats_butyrolactone_network.txt")

######################################################
# Visualize in a  network! ####
# Raw distance <= 0.3 ####
network <- read.table("Others_c1.00.network", header=TRUE, sep="\t")
but <- network[grep("butyrolactone - butyrolactone", network$Combined.group),] # grabs all rows that have the butyrolactone BGC
# #remove all raw-distances greater than 0.7
but <- filter(but,Raw.distance<0.3)

# retain all unique BGC names 
names <- data.frame(unique(rbind(data.frame(V1 = unique(but$Clustername.1)),
                                 data.frame(V1 = unique(but$Clustername.2)) )))

# read in the clan and family tsv
cf <- read.table("Others_clans_1.00_1.00.tsv", header=FALSE, sep="\t",stringsAsFactors=FALSE)

# keep only the entries with butyrolactone
cf <-right_join(cf, names, by = c("V1")) #merges by names that only have butyrolactone links

f = data.frame(from = but$Clustername.1,
               to = but$Clustername.2, 
               but[,3:4])

g <- graph_from_data_frame(f, directed=FALSE, vertices=cf)

plot(g, edge.arrow.size = 0.2)
plot(g, vertex.size = 10, vertex.color = "red", vertex.shape = "circle", 
     edge.width = 1, edge.color = "black", edge.lty = 3, 
     edge.label = f$Raw.distance, edge.curved = FALSE)

autograph(g)

V(g)$clu <- as.character(membership(cluster_louvain(g)))
V(g)$size <- degree(g)

# Assign colouurs to thegroups
sort(unique(V(g)$V3))
# more directed colouring
got_palette2 <- c("119" = "#1A5878", 
                  "123" = "#C44237", 
                  "127" = "#AD8951",
                  "242" = "#A99083", 
                  "251" = "#50594B", 
                  "304" = "#F8766D", 
                  "84" = "#9ACD32")

ggraph(g, layout = "stress", bbox = 10) +
  geom_edge_link0(edge_width = 0.5, edge_colour = "grey66",
                  edge_alpha= 0.8) +
  geom_node_point(aes(fill = V3, size = V2), shape = 21) +
  geom_node_text(aes(filter = size >= 1, label = name), 
                 family="arial", size = 5) +
  scale_fill_manual(values = got_palette2) +
  theme_graph() +
  theme(legend.position = "bottom")
