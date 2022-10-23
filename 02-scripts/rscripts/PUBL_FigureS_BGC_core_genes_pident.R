###############
# Title: Heatmap for core genes of butyrolactone BGC
###############
# Author: Anan Ibrahim = @darcy220606
###############
# Date: May,2022
###############
# Libraries required
###############
library('pheatmap')
library('dplyr') # if not installed already
library('readr') # if not installed already
library('RColorBrewer')
######################################################
# Set-WD ####
setwd("/media/aibrahim/INTENSO/Deep-Evo/final/afsa_details")
######################################################
# Heatmap of ANIs ####

# Uncomment if you use the script in the bash script
#arg = commandArgs(trailingOnly = TRUE)

files <- list.files(path="/media/aibrahim/INTENSO/Deep-Evo/final/afsa_details/pims",
                    pattern="*.txt", full.names=TRUE, recursive=FALSE)
# fileNames <- Sys.glob("*.csv")
lapply(files, function(x) {
  pims <- read.table(x, sep = "",nrows= 9 , 
                     na.strings ="", stringsAsFactors= F ) # import dataset seperated by white space/ # load file
  # retrieve the filename basename and remove the extension
  name = basename(x) 
  name <- gsub(".pim.txt","",name)
  
  pims <- pims %>% mutate_all(~gsub("_AFSA_aa", "", .)) # remove specific strings
  pims <- pims %>% mutate_all(~gsub("_AFSA_nt", "", .)) # remove specific strings
  pims <- pims %>% mutate_all(~gsub("_OXI_aa", "", .)) # remove specific strings
  pims <- pims %>% mutate_all(~gsub("_OXI_nt", "", .)) # remove specific strings
  pims <- pims %>% mutate_all(~gsub("_SKI_aa", "", .)) # remove specific strings
  pims <- pims %>% mutate_all(~gsub("_SKI_nt", "", .)) # remove specific strings
  pims <- pims %>% mutate_all(~gsub("_TERP_aa", "", .)) # remove specific strings
  pims <- pims %>% mutate_all(~gsub("_TERP_nt", "", .)) # remove specific strings
  
  pims[nrow(pims) + 1,] <- c("xx",pims$V1) # add a row to the dataframe and add the row names
  
  samples<- c("EMN001",
              "RIG001",
              "PLV001_001",
              "PES001",
              "GOY005",
              "TAF017",
              "GOY006",
              "PLV001_002")
  
  colnames(pims) <- c(pims[9,]) #change coluumnames according to a row 
  
  pims <- pims[match(samples, pims$xx),] #reorder according to age of the samples 
  
  rownames(pims) <- pims$xx # change the rownames 
  
  pims = pims[,-c(1)] #remove the first column 
  
  pims1 <- as.data.frame(sapply(pims, as.numeric)) # change to numeric for heatmap
  rownames(pims1) <- rownames(pims)
  
  pims1 = pims1 %>% select(samples) #reorder column names according to the ordered samples 
  
  pdf(file = paste(name,".pdf",sep = "",collapse = NULL))
  par(mfrow=c(4,2))
  col <- c("#F497B6","#E8749A","#D12E64", "#87153B","#5E0C27", "#97E4F4","#74C9E8", "#2D92D1","#174E86","#0C285E")
  pheatmap(pims1, display_numbers = T,color=col,
           main = name,
           show_rownames = TRUE,show_colnames = TRUE,fontsize  = 10,
           # width = 10,height = 10,
           # cellwidth=10, cellheight=10,
           cluster_cols = F,cluster_rows=F)
  dev.off()
})
######################################################
