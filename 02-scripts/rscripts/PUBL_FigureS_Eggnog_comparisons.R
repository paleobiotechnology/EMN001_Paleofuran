###############
# Title: Heatmap illustrating teh fastANI values of ancient and modern genomes
###############
# Author: Anan Ibrahim = darcy220606
###############
# Date: July 2022
###############
# Libraries required
library('tidyverse')
library('pheatmap')
library('readr') # if not installed already
library('dplyr')
library('plyr')
library('ggplot2')
library('data.table')
library('tidyr')
library('stringr')
library('RColorBrewer')
library('reshape2')
library('readxl')

######################################################
# Set-WD ####
setwd("F:/review_folder/photosynthesis")


######################################################
# Load the eggnog and prokka files in the global env. ####

#eggNOG tables:
files <- list.files(path="F:/review_folder/photosynthesis/eggnog",
                    pattern="test.emapper*", full.names=TRUE, recursive=FALSE)
lapply(files, 
       #import dataset seperated by white space/ # load file
       function(x) {
         # table <- read.table(x, sep = "\t",header=FALSE ,
         #                     na.strings ="", stringsAsFactors= F ) 
         table <- read_excel(x)
         # remove the first two rows and assign the first as header and rename the column to locus tag
         table <- table[-c(1),]
         names(table) <- table[1,]
         table <- table[-1,]
         colnames(table)[1] <- "locus_tag"

         #grab the names of the samples
         name = basename(x)
         name <- gsub("test.emapper.annotations_","",name)
         name
         #write dataframes to workspace in r
         assign(name, table, envir=.GlobalEnv)
       } )

#prokka tables:
files <- list.files(path="F:/review_folder/photosynthesis/prokka",
                    pattern="*_prokka", full.names=TRUE, recursive=FALSE)
lapply(files, 
       #import dataset seperated by white space/ # load file
       function(x) {
         table <- read.delim(x, sep = "\t",header=TRUE ,
                             na.strings ="", stringsAsFactors= F )
         name = basename(x)
         assign(name, table, envir=.GlobalEnv)
       } )

######################################################
# Merge tables from EggNOG analysis with Prokka annotations ####

EMN001_021 <- merge(x = EMN001_021_prokka.tsv, y =EMN001_021.xlsx, by="locus_tag", all = TRUE)
EMN001_021 <- EMN001_021[-c(1,2,3),]
PES001_018 <- merge(x = PES001_018_prokka.tsv, y =PES001_018.xlsx, by="locus_tag", all = TRUE)
PES001_018 <- PES001_018[-c(1,2,3),]
TAF017_RA <- merge(x = TAF017_RA_prokka.tsv, y =TAF017_RA.xlsx, by="locus_tag", all = TRUE)
TAF017_RA <- TAF017_RA[-c(1,2,3),]
GOY005_001 <- merge(x = GOY005_001_prokka.tsv, y =GOY005_001.xlsx, by="locus_tag", all = TRUE)
GOY005_001 <- GOY005_001[-c(1,2,3),]
GOY006_RA <- merge(x = GOY006_RA_prokka.tsv, y =GOY006_RA.xlsx, by="locus_tag", all = TRUE)
GOY006_RA <- GOY006_RA[-c(1,2,3),]
RIG001_014 <- merge(x = RIG001_014_prokka.tsv, y =RIG001_014.xlsx, by="locus_tag", all = TRUE)
RIG001_014 <- RIG001_014[-c(1,2,3),]
PLV001_001 <- merge(x = PLV001_001_prokka.tsv, y =PLV001_001.xlsx, by="locus_tag", all = TRUE)
PLV001_001 <- PLV001_001[-c(1,2,3),]
PLV001_002 <- merge(x = PLV001_002_prokka.tsv, y =PLV001_002.xlsx, by="locus_tag", all = TRUE)
PLV001_002 <- PLV001_002[-c(1,2,3),]

GCA_013335335.1 <- merge(x = GCA_013335335.1_prokka.tsv, y =GCA_013335335.1.xlsx, by="locus_tag", all = TRUE)
GCA_013335335.1 <- GCA_013335335.1[-c(1,2,3),]
GCA_013335765.1 <- merge(x = GCA_013335765.1_prokka.tsv, y =GCA_013335765.1.xlsx, by="locus_tag", all = TRUE)
GCA_013335765.1 <- GCA_013335765.1[-c(1,2,3),]
GCA_019163275.1 <- merge(x = GCA_019163275.1_prokka.tsv, y =GCA_019163275.1.xlsx, by="locus_tag", all = TRUE)
GCA_019163275.1 <- GCA_019163275.1[-c(1,2,3),]
GCF_000015125.1 <- merge(x = GCF_000015125.1_prokka.tsv, y =GCF_000015125.1.xlsx, by="locus_tag", all = TRUE)
GCF_000015125.1 <- GCF_000015125.1[-c(1,2,3),]
GCF_000020465.1 <- merge(x = GCF_000020465.1_prokka.tsv, y =GCF_000020465.1.xlsx, by="locus_tag", all = TRUE)
GCF_000020465.1 <- GCF_000020465.1[-c(1,2,3),]
GCF_001509575.1 <- merge(x = GCF_001509575.1_prokka.tsv, y =GCF_001509575.1.xlsx, by="locus_tag", all = TRUE)
GCF_001509575.1 <- GCF_001509575.1[-c(1,2,3),]

# write to output folder

write.csv(EMN001_021,"EMN001_021.csv")
write.csv(PES001_018,"PES001_018.csv")
write.csv(TAF017_RA,"TAF017_RA.csv")
write.csv(GOY005_001,"GOY005_001.csv")
write.csv(GOY006_RA,"GOY006_RA.csv")
write.csv(RIG001_014,"RIG001_014.csv")
write.csv(PLV001_001,"PLV001_001.csv")
write.csv(PLV001_002,"PLV001_002.csv")
write.csv(GCA_013335335.1,"GCA_013335335.1.csv")
write.csv(GCA_013335335.1,"GCA_013335335.1.csv")
write.csv(GCA_013335765.1,"GCA_013335765.1.csv")
write.csv(GCA_013335765.1,"GCA_013335765.1.csv")
write.csv(GCA_019163275.1,"GCA_019163275.1.csv")
write.csv(GCA_019163275.1,"GCA_019163275.1.csv")
write.csv(GCF_000015125.1,"GCF_000015125.1.csv")
write.csv(GCF_000015125.1,"GCF_000015125.1.csv")
write.csv(GCF_000020465.1,"GCF_000020465.1.csv")
write.csv(GCF_000020465.1,"GCF_000020465.1.csv")
write.csv(GCF_001509575.1,"GCF_001509575.1.csv")
write.csv(GCF_001509575.1,"GCF_001509575.1.csv")


######################################################
# Heatmap plot the photosynthetic genes from KO00194 in MAGs - FigureS ####
#load the KO genes associated with Ko00194 (photosynthetic genes)
KO <- read.table("F:/review_folder/photosynthesis/ko00194_kegg_pathway.txt", sep = ";", header=FALSE ,
                 na.strings ="", stringsAsFactors= F )
ko_var <- c(KO$V1)

files <- list.files(path="F:/review_folder/photosynthesis/eggnog",
                    pattern="test.emapper*", full.names=TRUE, recursive=FALSE)

# Convert to long format 
lapply(files, 
       function(x) {
         # table <- read.table(x, sep = "\t",header=FALSE ,
         #                     na.strings ="", stringsAsFactors= F ) 
         table <- read_excel(x)
         # remove the first two rows and assign the first as header and rename the column to locus tag
         table <- table[-c(1),]
         names(table) <- table[1,]
         table <- table[-1,]
         # retrieve the filename basename and remove the extension
         name = basename(x)
         name <- gsub("test.emapper.annotations_","",name)
         name <- gsub(".xlsx","",name)
         
         new_table <- table[ grepl(paste(ko_var, collapse="|"), table$KEGG_ko),]
      
         name_table = data.frame(name = new_table$KEGG_ko)
         name_table = data.frame(name = table(name_table))
         name_table = data.frame(name_table, name=c(rep(name,nrow(name_table))))
         colnames(name_table) = c("ko", "freq", "name")
         
         name_table_gather <- name_table %>%
           gather(-ko, -name, key = ko, value = freq) %>% #join by (-), 
           group_by(ko) %>%
           # drop_na(otunum) %>%
           ungroup()
         
         assign(name, name_table, envir=.GlobalEnv) #write dataframes to workspace in r
       } )

# Grab the rows that contain the Ko (photosynthetic genes)

m1 = merge(x = EMN001_021, y =PES001_018, all = TRUE)
m2 = merge(x = RIG001_014, y =TAF017_RA, all = TRUE)
m3 = merge(x = PLV001_002, y =PLV001_001, all = TRUE)
m4 = merge(x = GOY005_001, y =GOY006_RA, all = TRUE)

m5 = merge(x = m1, y = m2, all = TRUE)
m6 = merge(x = m3, y = m4, all = TRUE)

m13 = merge(x = m5, y = m6, all = TRUE)

m8 = merge(x = GCA_013335335.1, y =GCA_013335765.1, all = TRUE)
m9 = merge(x = GCA_019163275.1, y =GCF_000015125.1, all = TRUE)
m10 = merge(x = GCF_000020465.1, y =GCF_001509575.1, all = TRUE)

m11 = merge(x = m8, y =m9, all = TRUE)
m12 = merge(x = m11, y =m10, all = TRUE)

m7 = merge(x = m12, y =m13, all = TRUE)

# remove the suffix from query
#m7$query = stringr::str_extract(m7$query, "[^_]*_[^_]*")
# remove the prefix from ko
m7$ko <- gsub('ko:', '', m7$ko)

m7$ko <- factor(m7$ko, levels=unique(m7$ko)) 

# merge the ko identifications
m7 =merge (m7, KO, by.x = 'ko', by.y ='V1')

name_order <- c('RIG001_014', 
                'EMN001_021', 
                'PLV001_001', 
                'TAF017_RA',
                'GOY005_001',
                'GOY006_RA',
                'PES001_018',
                'PLV001_002',
                'GCF_000020465.1',
                'GCA_013335765.1',
                'GCA_013335335.1',
                'GCF_001509575.1',
                'GCA_019163275.1',
                'GCF_000015125.1') 

# plot heatmap
# heatmap:::
pdf("photosynthesis_kegg_genes_heatmap_ABBREV.pdf",width=20,height=20)
ggplot(m7, aes(y=factor(name, level = name_order), x=factor(ko))) +
  geom_tile(aes(fill=freq), colour="Black", lwd =0.1, linetype=1) +
  # geom_text(aes(label = fre), color = "black", size = 1) +
  scale_fill_gradient(low = "#F497B6", high = "#0C285E") +
  guides(fill = guide_colourbar(title = "Gene proportion %",ticks = TRUE)) +
  coord_fixed() +
  xlab("KO genes") + ylab("Genomes") +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))
# theme_bw()
dev.off()

pdf("photosynthesis_kegg_genes_heatmap_IDS.pdf",width=20,height=20)
ggplot(m7, aes(y=factor(name, level = name_order), x=factor(V2))) +
  geom_tile(aes(fill=freq), colour="Black", lwd =0.1, linetype=1) +
  # geom_text(aes(label = fre), color = "black", size = 1) +
  scale_fill_gradient(low = "#F497B6", high = "#0C285E") +
  guides(fill = guide_colourbar(title = "Gene proportion %",ticks = TRUE)) +
  coord_fixed() +
  xlab("KO genes") + ylab("Genomes") +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))
# theme_bw()
dev.off()

######################################################
# Calculate the photosynthesis genes on the butyrolactone BGC contig ####
# Tabulate the Contig pyrothrin pathways : contigs that contain the BGC in a table
# AFSA gene number : Manually insert according to table
# Note:: Cannot be automated as old .gbk files do not follow the gbk format as they were not run with --compliant 
EMN021_021_AFSA = "00950"
PES001_018_AFSA = "00504"
GOY006_RA_AFSA = "00924"
GOY005_001_12_AFSA = "00900" 
GOY005_001_11_AFSA =  "01203"
TAF017_RA_AFSA = "00934"
RIG001_014_AFSA = "00328"
PLV001_002_AFSA = "01140"
PLV001_001_AFSA =  "00554"
GCF_000020465.1_AFSA = "00914"
GCF_000015125.1_AFSA = "01900"
GCF_001509575.1_AFSA = "00290"
GCA_013335765.1_AFSA = "00659"
# GCA_019163275.1_AFSA = "" #cant find it at all
GCA_013335335.1_AFSA = "02058"

# Gene number list (on BGC contig): : Manually insert according to gff file
# Note:: Cannot be automated as old .gbk files do not follow the gbk format as they were not run with --compliant 
#EMN
list = as.list.data.frame(00946:01038) # list the numbers in sequential order 
list=  formatC(list, width = 5, format = "d", flag = "0") #prevent r from removing the zeros prefixes : fixed lengths width 5
EMN021_021_list = paste(list, collapse= "|") 
#PES
list = as.list.data.frame(00495:00592) # list the numbers in sequential order 
list =  formatC(list, width = 5, format = "d", flag = "0") #prevent r from removing the zeros prefixes : fixed lengths width 5
PES001_018_list = paste(list, collapse= "|") 
#GOY006
list = as.list.data.frame(00920:01006) # list the numbers in sequential order 
list =  formatC(list, width = 5, format = "d", flag = "0") #prevent r from removing the zeros prefixes : fixed lengths width 5
GOY006_RA_list = paste(list, collapse= "|") 
#GOY005_11
list = as.list.data.frame(01199:01204) # list the numbers in sequential order 
list =  formatC(list, width = 5, format = "d", flag = "0") #prevent r from removing the zeros prefixes : fixed lengths width 5
GOY005_001_11_list = paste(list, collapse= "|") 
#GOY005_12
list = as.list.data.frame(00900:00910) # list the numbers in sequential order 
list =  formatC(list, width = 5, format = "d", flag = "0") #prevent r from removing the zeros prefixes : fixed lengths width 5
GOY005_001_12_list = paste(list, collapse= "|") 
#TAF017_RA
list = as.list.data.frame(00930:01016) # list the numbers in sequential order 
list =  formatC(list, width = 5, format = "d", flag = "0") #prevent r from removing the zeros prefixes : fixed lengths width 5
TAF017_RA_list = paste(list, collapse= "|") 
#RIG001
list = as.list.data.frame(00324:00430) # list the numbers in sequential order 
list =  formatC(list, width = 5, format = "d", flag = "0") #prevent r from removing the zeros prefixes : fixed lengths width 5
RIG001_014_list = paste(list, collapse= "|") 
#PLV001_002
list = as.list.data.frame(01135:01145) # list the numbers in sequential order 
list =  formatC(list, width = 5, format = "d", flag = "0") #prevent r from removing the zeros prefixes : fixed lengths width 5
PLV001_002_list = paste(list, collapse= "|") 
#PLV001_001
list = as.list.data.frame(00550:00558) # list the numbers in sequential order 
list =  formatC(list, width = 5, format = "d", flag = "0") #prevent r from removing the zeros prefixes : fixed lengths width 5
PLV001_001_list = paste(list, collapse= "|") 
#GCF_000020465.1
list = as.list.data.frame(00001:02550) # list the numbers in sequential order 
list =  formatC(list, width = 5, format = "d", flag = "0") #prevent r from removing the zeros prefixes : fixed lengths width 5
GCF_000020465.1_list = paste(list, collapse= "|") 
#GCF_000015125.1
list = as.list.data.frame(00001:02550) # list the numbers in sequential order 
list =  formatC(list, width = 5, format = "d", flag = "0") #prevent r from removing the zeros prefixes : fixed lengths width 5
GCF_000015125.1_list = paste(list, collapse= "|") 
#GCA_013335765.1
list = as.list.data.frame(00498:00682) # list the numbers in sequential order 
list =  formatC(list, width = 5, format = "d", flag = "0") #prevent r from removing the zeros prefixes : fixed lengths width 5
GCA_013335765.1_list = paste(list, collapse= "|") 
#GCF_001509575.1
list = as.list.data.frame(00281:00290) # list the numbers in sequential order 
list =  formatC(list, width = 5, format = "d", flag = "0") #prevent r from removing the zeros prefixes : fixed lengths width 5
GCF_001509575.1_list = paste(list, collapse= "|") 
#GCA_013335335.1
list = as.list.data.frame(02054:02059) # list the numbers in sequential order 
list =  formatC(list, width = 5, format = "d", flag = "0") #prevent r from removing the zeros prefixes : fixed lengths width 5
GCA_013335335.1_list = paste(list, collapse= "|") 

# Run table filtering
files_list <- list.files(path="/media/aibrahim/INTENSO/Deep-Evo/final/EGG_NOG-ancient-mags",
                         pattern="test.emapper.annotations_*", full.names=TRUE, recursive=FALSE)

for (file in files_list){
  if(grepl("test.emapper.annotations_EMN021", file)){
    lapply(file,
           function(x) {
             name = basename(x)
             name <- gsub("test.emapper.annotations_","",name)
             table <- read.table(x, header=FALSE, sep="\t")
             table$V1 <- gsub("^EMN001_021_", "", table$V1) # remove part of a string form a column in df
             # Step1 : Adjust the table according to the contig annotation list
             table1 = table[grepl(EMN021_021_list, table$V1),]
             # Step2 : grep only the KEGG pathways (map...) from the contig genes
             table1 <- data.frame(gene_id=table1$V1, pathways=table1$V13)
             table1$pathways <- strsplit(as.character(table1$pathways), split = ",",fixed=FALSE) #split a column based on commas
             table1 <- tidyr::unnest(table1, cols = pathways) #unlist a list of strings in a column made of a list
             table1 <- table1[- grep("ko", table1$pathways),] # remove the ko pathways
             genes <- table1[grep("map00860", table1$pathways),] # grabs all rows that have the pyroph genes
             positions <- data.frame(unique(genes$gene_id))
             # positions <- paste(positions$unique.genes.gene_id.,sep = ",")
             positions <- data.frame(val=paste0(positions$unique.genes.gene_id.,collapse = ', '),stringsAsFactors = F)
             # Step3 : Count the number of time perphorine metabolism occured (i.e. how many genes) on the contig
             sum1 <- length(unique(table1$gene_id))
             sum2 <- sum(str_count(table1, "map00860"))
             # STep4 : Get the proportion of the Porphyrine genes in the whole bin
             table2 <- data.frame(gene_id=table$V1, pathways=table$V13)
             table2$pathways <- strsplit(as.character(table2$pathways), split = ",",fixed=FALSE) #split a column based on commas
             table2 <- tidyr::unnest(table2, cols = pathways) #unlist a list of strings in a column made of a list
             table2 <- table2[- grep("ko", table2$pathways),]
             sum4 <- sum(str_count(table2, "map00860"))
             name1 <- data.frame(sample_name = name,
                                 prop_of_pyroph_in_tot_pyroph_in_mag = (sum2/sum4)*100,
                                 prop_of_pyroph_in_contig = (sum2/sum1)*100,
                                 sum_of_genes_in_bin = length(unique(table2$gene_id)),
                                 sum_of_genes_in_contig = length(unique(table1$gene_id)),
                                 AFSA_gene_number = EMN021_021_AFSA,
                                 gene_id_of_pyroph = c(positions)
             )
             print(name1)
             assign(name,name1, envir=.GlobalEnv)
           })
  } }
for (file in files_list){
  if(grepl("test.emapper.annotations_PES001", file)){
    lapply(file,
           function(x) {
             name = basename(x)
             name <- gsub("test.emapper.annotations_","",name)
             table <- read.delim(x, header=FALSE, sep="\t")
             
             table$V1 <- gsub("^PES001_018_", "", table$V1) # remove part of a string form a column in df
             # Step1 : Adjust the table according to the contig annotation list
             table1 = table[grepl(PES001_018_list, table$V1),]
             # Step2 : grep only the KEGG pathways (map...) from the contig genes
             table1 <- data.frame(gene_id=table1$V1, pathways=table1$V13)
             table1$pathways <- strsplit(as.character(table1$pathways), split = ",",fixed=FALSE) #split a column based on commas
             table1 <- tidyr::unnest(table1, cols = pathways) #unlist a list of strings in a column made of a list
             table1 <- table1[- grep("ko", table1$pathways),] # remove the ko pathways
             genes <- table1[grep("map00860", table1$pathways),] # grabs all rows that have the pyroph genes
             positions <- data.frame(unique(genes$gene_id))
             # positions <- paste(positions$unique.genes.gene_id.,sep = ",")
             positions <- data.frame(val=paste0(positions$unique.genes.gene_id.,collapse = ', '),stringsAsFactors = F)
             # Step3 : Count the number of time perphorine metabolism occured (i.e. how many genes) on the contig
             sum1 <- length(unique(table1$gene_id))
             sum2 <- sum(str_count(table1, "map00860"))
             # STep4 : Get the proportion of the Porphyrine genes in the whole bin
             table2 <- data.frame(gene_id=table$V1, pathways=table$V13)
             table2$pathways <- strsplit(as.character(table2$pathways), split = ",",fixed=FALSE) #split a column based on commas
             table2 <- tidyr::unnest(table2, cols = pathways) #unlist a list of strings in a column made of a list
             table2 <- table2[- grep("ko", table2$pathways),]
             sum4 <- sum(str_count(table2, "map00860"))
             name1 <- data.frame(sample_name = name,
                                 prop_of_pyroph_in_tot_pyroph_in_mag = (sum2/sum4)*100,
                                 prop_of_pyroph_in_contig = (sum2/sum1)*100,
                                 sum_of_genes_in_bin = length(unique(table2$gene_id)),
                                 sum_of_genes_in_contig = length(unique(table1$gene_id)),
                                 AFSA_gene_number = PES001_018_AFSA,
                                 gene_id_of_pyroph = c(positions)
             )
             print(name1)
             assign(name,name1, envir=.GlobalEnv)
           })
  } }
for (file in files_list){
  if(grepl("test.emapper.annotations_GOY006", file)){
    lapply(file,
           function(x) {
             name = basename(x)
             name <- gsub("test.emapper.annotations_","",name)
             table <- read.delim(x, header=FALSE, sep="\t")
             table$V1 <- gsub("^GOY006_RA_", "", table$V1) # remove part of a string form a column in df
             # Step1 : Adjust the table according to the contig annotation list
             table1 = table[grepl(GOY006_RA_list, table$V1),]
             # Step2 : grep only the KEGG pathways (map...) from the contig genes
             table1 <- data.frame(gene_id=table1$V1, pathways=table1$V13)
             table1$pathways <- strsplit(as.character(table1$pathways), split = ",",fixed=FALSE) #split a column based on commas
             table1 <- tidyr::unnest(table1, cols = pathways) #unlist a list of strings in a column made of a list
             table1 <- table1[- grep("ko", table1$pathways),] # remove the ko pathways
             genes <- table1[grep("map00860", table1$pathways),] # grabs all rows that have the pyroph genes
             positions <- data.frame(unique(genes$gene_id))
             # positions <- paste(positions$unique.genes.gene_id.,sep = ",")
             positions <- data.frame(val=paste0(positions$unique.genes.gene_id.,collapse = ', '),stringsAsFactors = F)
             # Step3 : Count the number of time perphorine metabolism occured (i.e. how many genes) on the contig
             sum1 <- length(unique(table1$gene_id))
             sum2 <- sum(str_count(table1, "map00860"))
             # STep4 : Get the proportion of the Porphyrine genes in the whole bin
             table2 <- data.frame(gene_id=table$V1, pathways=table$V13)
             table2$pathways <- strsplit(as.character(table2$pathways), split = ",",fixed=FALSE) #split a column based on commas
             table2 <- tidyr::unnest(table2, cols = pathways) #unlist a list of strings in a column made of a list
             table2 <- table2[- grep("ko", table2$pathways),]
             sum4 <- sum(str_count(table2, "map00860"))
             name1 <- data.frame(sample_name = name,
                                 prop_of_pyroph_in_tot_pyroph_in_mag = (sum2/sum4)*100,
                                 prop_of_pyroph_in_contig = (sum2/sum1)*100,
                                 sum_of_genes_in_bin = length(unique(table2$gene_id)),
                                 sum_of_genes_in_contig = length(unique(table1$gene_id)),
                                 AFSA_gene_number = GOY006_RA_AFSA,
                                 gene_id_of_pyroph = c(positions)
             )
             print(name1)
             assign(name,name1, envir=.GlobalEnv)
           })
  } }
for (file in files_list){
  if(grepl("test.emapper.annotations_TAF017", file)){
    lapply(file,
           function(x) {
             name = basename(x)
             name <- gsub("test.emapper.annotations_","",name)
             table <- read.delim(x, header=FALSE, sep="\t")
             table$V1 <- gsub("^TAF017_RA_", "", table$V1) # remove part of a string form a column in df
             # Step1 : Adjust the table according to the contig annotation list
             table1 = table[grepl(TAF017_RA_list, table$V1),]
             # Step2 : grep only the KEGG pathways (map...) from the contig genes
             table1 <- data.frame(gene_id=table1$V1, pathways=table1$V13)
             table1$pathways <- strsplit(as.character(table1$pathways), split = ",",fixed=FALSE) #split a column based on commas
             table1 <- tidyr::unnest(table1, cols = pathways) #unlist a list of strings in a column made of a list
             table1 <- table1[- grep("ko", table1$pathways),] # remove the ko pathways
             genes <- table1[grep("map00860", table1$pathways),] # grabs all rows that have the pyroph genes
             positions <- data.frame(unique(genes$gene_id))
             # positions <- paste(positions$unique.genes.gene_id.,sep = ",")
             positions <- data.frame(val=paste0(positions$unique.genes.gene_id.,collapse = ', '),stringsAsFactors = F)
             # Step3 : Count the number of time perphorine metabolism occured (i.e. how many genes) on the contig
             sum1 <- length(unique(table1$gene_id))
             sum2 <- sum(str_count(table1, "map00860"))
             # STep4 : Get the proportion of the Porphyrine genes in the whole bin
             table2 <- data.frame(gene_id=table$V1, pathways=table$V13)
             table2$pathways <- strsplit(as.character(table2$pathways), split = ",",fixed=FALSE) #split a column based on commas
             table2 <- tidyr::unnest(table2, cols = pathways) #unlist a list of strings in a column made of a list
             table2 <- table2[- grep("ko", table2$pathways),]
             sum4 <- sum(str_count(table2, "map00860"))
             name1 <- data.frame(sample_name = name,
                                 prop_of_pyroph_in_tot_pyroph_in_mag = (sum2/sum4)*100,
                                 prop_of_pyroph_in_contig = (sum2/sum1)*100,
                                 sum_of_genes_in_bin = length(unique(table2$gene_id)),
                                 sum_of_genes_in_contig = length(unique(table1$gene_id)),
                                 AFSA_gene_number = TAF017_RA_AFSA,
                                 gene_id_of_pyroph = c(positions)
             )
             print(name1)
             assign(name,name1, envir=.GlobalEnv)
           })
  } }
for (file in files_list){
  if(grepl("test.emapper.annotations_RIG001", file)){
    lapply(file,
           function(x) {
             name = basename(x)
             name <- gsub("test.emapper.annotations_","",name)
             table <- read.delim(x, header=FALSE, sep="\t")
             table$V1 <- gsub("^RIG001_014_", "", table$V1) # remove part of a string form a column in df
             # Step1 : Adjust the table according to the contig annotation list
             table1 = table[grepl(RIG001_014_list, table$V1),]
             # Step2 : grep only the KEGG pathways (map...) from the contig genes
             table1 <- data.frame(gene_id=table1$V1, pathways=table1$V13)
             table1$pathways <- strsplit(as.character(table1$pathways), split = ",",fixed=FALSE) #split a column based on commas
             table1 <- tidyr::unnest(table1, cols = pathways) #unlist a list of strings in a column made of a list
             table1 <- table1[- grep("ko", table1$pathways),] # remove the ko pathways
             genes <- table1[grep("map00860", table1$pathways),] # grabs all rows that have the pyroph genes
             positions <- data.frame(unique(genes$gene_id))
             # positions <- paste(positions$unique.genes.gene_id.,sep = ",")
             positions <- data.frame(val=paste0(positions$unique.genes.gene_id.,collapse = ', '),stringsAsFactors = F)
             # Step3 : Count the number of time perphorine metabolism occured (i.e. how many genes) on the contig
             sum1 <- length(unique(table1$gene_id))
             sum2 <- sum(str_count(table1, "map00860"))
             # STep4 : Get the proportion of the Porphyrine genes in the whole bin
             table2 <- data.frame(gene_id=table$V1, pathways=table$V13)
             table2$pathways <- strsplit(as.character(table2$pathways), split = ",",fixed=FALSE) #split a column based on commas
             table2 <- tidyr::unnest(table2, cols = pathways) #unlist a list of strings in a column made of a list
             table2 <- table2[- grep("ko", table2$pathways),]
             sum4 <- sum(str_count(table2, "map00860"))
             name1 <- data.frame(sample_name = name,
                                 prop_of_pyroph_in_tot_pyroph_in_mag = (sum2/sum4)*100,
                                 prop_of_pyroph_in_contig = (sum2/sum1)*100,
                                 sum_of_genes_in_bin = length(unique(table2$gene_id)),
                                 sum_of_genes_in_contig = length(unique(table1$gene_id)),
                                 AFSA_gene_number = RIG001_014_AFSA,
                                 gene_id_of_pyroph = c(positions)
             )
             print(name1)
             assign(name,name1, envir=.GlobalEnv)
           })
  } }
for (file in files_list){
  if(grepl("test.emapper.annotations_PLV001_002", file)){
    lapply(file,
           function(x) {
             name = basename(x)
             name <- gsub("test.emapper.annotations_","",name)
             table <- read.delim(x, header=FALSE, sep="\t")
             table$V1 <- gsub("^PLV001_002_", "", table$V1) # remove part of a string form a column in df
             # Step1 : Adjust the table according to the contig annotation list
             table1 = table[grepl(PLV001_002_list, table$V1),]
             # Step2 : grep only the KEGG pathways (map...) from the contig genes
             table1 <- data.frame(gene_id=table1$V1, pathways=table1$V13)
             table1$pathways <- strsplit(as.character(table1$pathways), split = ",",fixed=FALSE) #split a column based on commas
             table1 <- tidyr::unnest(table1, cols = pathways) #unlist a list of strings in a column made of a list
             table1 <- table1[- grep("ko", table1$pathways),] # remove the ko pathways
             genes <- table1[grep("map00860", table1$pathways),] # grabs all rows that have the pyroph genes
             positions <- data.frame(unique(genes$gene_id))
             # positions <- paste(positions$unique.genes.gene_id.,sep = ",")
             positions <- data.frame(val=paste0(positions$unique.genes.gene_id.,collapse = ', '),stringsAsFactors = F)
             # Step3 : Count the number of time perphorine metabolism occured (i.e. how many genes) on the contig
             sum1 <- length(unique(table1$gene_id))
             sum2 <- sum(str_count(table1, "map00860"))
             # STep4 : Get the proportion of the Porphyrine genes in the whole bin
             table2 <- data.frame(gene_id=table$V1, pathways=table$V13)
             table2$pathways <- strsplit(as.character(table2$pathways), split = ",",fixed=FALSE) #split a column based on commas
             table2 <- tidyr::unnest(table2, cols = pathways) #unlist a list of strings in a column made of a list
             table2 <- table2[- grep("ko", table2$pathways),]
             sum4 <- sum(str_count(table2, "map00860"))
             name1 <- data.frame(sample_name = name,
                                 prop_of_pyroph_in_tot_pyroph_in_mag = (sum2/sum4)*100,
                                 prop_of_pyroph_in_contig = (sum2/sum1)*100,
                                 sum_of_genes_in_bin = length(unique(table2$gene_id)),
                                 sum_of_genes_in_contig = length(unique(table1$gene_id)),
                                 AFSA_gene_number = PLV001_002_AFSA,
                                 gene_id_of_pyroph = c(positions)
             )
             print(name1)
             assign(name,name1, envir=.GlobalEnv)
           })
  } }
for (file in files_list){
  if(grepl("test.emapper.annotations_PLV001_001", file)){
    lapply(file,
           function(x) {
             name = basename(x)
             name <- gsub("test.emapper.annotations_","",name)
             table <- read.delim(x, header=FALSE, sep="\t")
             table$V1 <- gsub("^PLV001_001_", "", table$V1) # remove part of a string form a column in df
             # Step1 : Adjust the table according to the contig annotation list
             table1 = table[grepl(PLV001_001_list, table$V1),]
             # Step2 : grep only the KEGG pathways (map...) from the contig genes
             table1 <- data.frame(gene_id=table1$V1, pathways=table1$V13)
             table1$pathways <- strsplit(as.character(table1$pathways), split = ",",fixed=FALSE) #split a column based on commas
             table1 <- tidyr::unnest(table1, cols = pathways) #unlist a list of strings in a column made of a list
             table1 <- table1[- grep("ko", table1$pathways),] # remove the ko pathways
             genes <- table1[grep("map00860", table1$pathways),] # grabs all rows that have the pyroph genes
             positions <- data.frame(unique(genes$gene_id))
             # positions <- paste(positions$unique.genes.gene_id.,sep = ",")
             positions <- data.frame(val=paste0(positions$unique.genes.gene_id.,collapse = ', '),stringsAsFactors = F)
             # Step3 : Count the number of time perphorine metabolism occured (i.e. how many genes) on the contig
             sum1 <- length(unique(table1$gene_id))
             sum2 <- sum(str_count(table1, "map00860"))
             # STep4 : Get the proportion of the Porphyrine genes in the whole bin
             table2 <- data.frame(gene_id=table$V1, pathways=table$V13)
             table2$pathways <- strsplit(as.character(table2$pathways), split = ",",fixed=FALSE) #split a column based on commas
             table2 <- tidyr::unnest(table2, cols = pathways) #unlist a list of strings in a column made of a list
             table2 <- table2[- grep("ko", table2$pathways),]
             sum4 <- sum(str_count(table2, "map00860"))
             name1 <- data.frame(sample_name = name,
                                 prop_of_pyroph_in_tot_pyroph_in_mag = (sum2/sum4)*100,
                                 prop_of_pyroph_in_contig = (sum2/sum1)*100,
                                 sum_of_genes_in_bin = length(unique(table2$gene_id)),
                                 sum_of_genes_in_contig = length(unique(table1$gene_id)),
                                 AFSA_gene_number = PLV001_001_AFSA,
                                 gene_id_of_pyroph = c(positions)
             )
             print(name1)
             assign(name,name1, envir=.GlobalEnv)
           })
  } }
for (file in files_list){
  if(grepl("test.emapper.annotations_GOY005", file)){
    lapply(file,
           function(x) {
             name = basename(x)
             name <- gsub("test.emapper.annotations_","12_",name)
             table <- read.delim(x, header=FALSE, sep="\t")
             table$V1 <- gsub("^GOY005_001_", "", table$V1) # remove part of a string form a column in df
             # Step1 : Adjust the table according to the contig annotation list
             table1 = table[grepl(GOY005_001_12_list, table$V1),]
             # Step2 : grep only the KEGG pathways (map...) from the contig genes
             table1 <- data.frame(gene_id=table1$V1, pathways=table1$V13)
             table1$pathways <- strsplit(as.character(table1$pathways), split = ",",fixed=FALSE) #split a column based on commas
             table1 <- tidyr::unnest(table1, cols = pathways) #unlist a list of strings in a column made of a list
             table1 <- table1[- grep("ko", table1$pathways),] # remove the ko pathways
             genes <- table1[grep("map00860", table1$pathways),] # grabs all rows that have the pyroph genes
             positions <- data.frame(unique(genes$gene_id))
             # positions <- paste(positions$unique.genes.gene_id.,sep = ",")
             positions <- data.frame(val=paste0(positions$unique.genes.gene_id.,collapse = ', '),stringsAsFactors = F)
             # Step3 : Count the number of time perphorine metabolism occured (i.e. how many genes) on the contig
             sum1 <- length(unique(table1$gene_id))
             sum2 <- sum(str_count(table1, "map00860"))
             # STep4 : Get the proportion of the Porphyrine genes in the whole bin
             table2 <- data.frame(gene_id=table$V1, pathways=table$V13)
             table2$pathways <- strsplit(as.character(table2$pathways), split = ",",fixed=FALSE) #split a column based on commas
             table2 <- tidyr::unnest(table2, cols = pathways) #unlist a list of strings in a column made of a list
             table2 <- table2[- grep("ko", table2$pathways),]
             sum4 <- sum(str_count(table2, "map00860"))
             name1 <- data.frame(sample_name = name,
                                 prop_of_pyroph_in_tot_pyroph_in_mag = (sum2/sum4)*100,
                                 prop_of_pyroph_in_contig = (sum2/sum1)*100,
                                 sum_of_genes_in_bin = length(unique(table2$gene_id)),
                                 sum_of_genes_in_contig = length(unique(table1$gene_id)),
                                 AFSA_gene_number = GOY005_001_12_AFSA,
                                 gene_id_of_pyroph = c(positions)
             )
             print(name1)
             assign(name,name1, envir=.GlobalEnv)
           })
  } }
for (file in files_list){
  if(grepl("test.emapper.annotations_GOY005", file)){
    lapply(file,
           function(x) {
             name = basename(x)
             name <- gsub("test.emapper.annotations_","11_",name)
             table <- read.delim(x, header=FALSE, sep="\t")
             table$V1 <- gsub("^GOY005_001_", "", table$V1) # remove part of a string form a column in df
             # Step1 : Adjust the table according to the contig annotation list
             table1 = table[grepl(GOY005_001_11_list, table$V1),]
             # Step2 : grep only the KEGG pathways (map...) from the contig genes
             table1 <- data.frame(gene_id=table1$V1, pathways=table1$V13)
             table1$pathways <- strsplit(as.character(table1$pathways), split = ",",fixed=FALSE) #split a column based on commas
             table1 <- tidyr::unnest(table1, cols = pathways) #unlist a list of strings in a column made of a list
             table1 <- table1[- grep("ko", table1$pathways),] # remove the ko pathways
             genes <- table1[grep("map00860", table1$pathways),] # grabs all rows that have the pyroph genes
             positions <- data.frame(unique(genes$gene_id))
             # positions <- paste(positions$unique.genes.gene_id.,sep = ",")
             positions <- data.frame(val=paste0(positions$unique.genes.gene_id.,collapse = ', '),stringsAsFactors = F)
             # Step3 : Count the number of time perphorine metabolism occured (i.e. how many genes) on the contig
             sum1 <- length(unique(table1$gene_id))
             sum2 <- sum(str_count(table1, "map00860"))
             # STep4 : Get the proportion of the Porphyrine genes in the whole bin
             table2 <- data.frame(gene_id=table$V1, pathways=table$V13)
             table2$pathways <- strsplit(as.character(table2$pathways), split = ",",fixed=FALSE) #split a column based on commas
             table2 <- tidyr::unnest(table2, cols = pathways) #unlist a list of strings in a column made of a list
             table2 <- table2[- grep("ko", table2$pathways),]
             sum4 <- sum(str_count(table2, "map00860"))
             name1 <- data.frame(sample_name = name,
                                 prop_of_pyroph_in_tot_pyroph_in_mag = (sum2/sum4)*100,
                                 prop_of_pyroph_in_contig = (sum2/sum1)*100,
                                 sum_of_genes_in_bin = length(unique(table2$gene_id)),
                                 sum_of_genes_in_contig = length(unique(table1$gene_id)),
                                 AFSA_gene_number = GOY005_001_11_AFSA,
                                 gene_id_of_pyroph = c(positions)
             )
             print(name1)
             assign(name,name1, envir=.GlobalEnv)
           })
  } }

for (file in files_list){
  if(grepl("test.emapper.annotations_GCF_000020465.1", file)){
    lapply(file,
           function(x) {
             name = basename(x)
             name <- gsub("test.emapper.annotations_","",name)
             table <- read.delim(x, header=FALSE, sep="\t")
             table$V1 <- gsub("^GCF_000020465.1_", "", table$V1) # remove part of a string form a column in df
             # Step1 : Adjust the table according to the contig annotation list
             table1 = table[grepl(GCF_000020465.1_list, table$V1),]
             # Step2 : grep only the KEGG pathways (map...) from the contig genes
             table1 <- data.frame(gene_id=table1$V1, pathways=table1$V13)
             table1$pathways <- strsplit(as.character(table1$pathways), split = ",",fixed=FALSE) #split a column based on commas
             table1 <- tidyr::unnest(table1, cols = pathways) #unlist a list of strings in a column made of a list
             table1 <- table1[- grep("ko", table1$pathways),] # remove the ko pathways
             genes <- table1[grep("map00860", table1$pathways),] # grabs all rows that have the pyroph genes
             positions <- data.frame(unique(genes$gene_id))
             # positions <- paste(positions$unique.genes.gene_id.,sep = ",")
             positions <- data.frame(val=paste0(positions$unique.genes.gene_id.,collapse = ', '),stringsAsFactors = F)
             # Step3 : Count the number of time perphorine metabolism occured (i.e. how many genes) on the contig
             sum1 <- length(unique(table1$gene_id))
             sum2 <- sum(str_count(table1, "map00860"))
             # STep4 : Get the proportion of the Porphyrine genes in the whole bin
             table2 <- data.frame(gene_id=table$V1, pathways=table$V13)
             table2$pathways <- strsplit(as.character(table2$pathways), split = ",",fixed=FALSE) #split a column based on commas
             table2 <- tidyr::unnest(table2, cols = pathways) #unlist a list of strings in a column made of a list
             table2 <- table2[- grep("ko", table2$pathways),]
             sum4 <- sum(str_count(table2, "map00860"))
             name1 <- data.frame(sample_name = name,
                                 prop_of_pyroph_in_tot_pyroph_in_mag = (sum2/sum4)*100,
                                 prop_of_pyroph_in_contig = (sum2/sum1)*100,
                                 sum_of_genes_in_bin = length(unique(table2$gene_id)),
                                 sum_of_genes_in_contig = length(unique(table1$gene_id)),
                                 AFSA_gene_number = GCF_000020465.1_AFSA, 
                                 gene_id_of_pyroph = c(positions)
             )
             print(name1)
             assign(name,name1, envir=.GlobalEnv)
           })
  } }
for (file in files_list){
  if(grepl("test.emapper.annotations_GCF_000015125.1", file)){
    lapply(file,
           function(x) {
             name = basename(x)
             name <- gsub("test.emapper.annotations_","",name)
             table <- read.delim(x, header=FALSE, sep="\t")
             table$V1 <- gsub("^GCF_000015125.1_", "", table$V1) # remove part of a string form a column in df
             # Step1 : Adjust the table according to the contig annotation list
             table1 = table[grepl(GCF_000015125.1_list, table$V1),]
             # Step2 : grep only the KEGG pathways (map...) from the contig genes
             table1 <- data.frame(gene_id=table1$V1, pathways=table1$V13)
             table1$pathways <- strsplit(as.character(table1$pathways), split = ",",fixed=FALSE) #split a column based on commas
             table1 <- tidyr::unnest(table1, cols = pathways) #unlist a list of strings in a column made of a list
             table1 <- table1[- grep("ko", table1$pathways),] # remove the ko pathways
             genes <- table1[grep("map00860", table1$pathways),] # grabs all rows that have the pyroph genes
             positions <- data.frame(unique(genes$gene_id))
             # positions <- paste(positions$unique.genes.gene_id.,sep = ",")
             positions <- data.frame(val=paste0(positions$unique.genes.gene_id.,collapse = ', '),stringsAsFactors = F)
             # Step3 : Count the number of time perphorine metabolism occured (i.e. how many genes) on the contig
             sum1 <- length(unique(table1$gene_id))
             sum2 <- sum(str_count(table1, "map00860"))
             # STep4 : Get the proportion of the Porphyrine genes in the whole bin
             table2 <- data.frame(gene_id=table$V1, pathways=table$V13)
             table2$pathways <- strsplit(as.character(table2$pathways), split = ",",fixed=FALSE) #split a column based on commas
             table2 <- tidyr::unnest(table2, cols = pathways) #unlist a list of strings in a column made of a list
             table2 <- table2[- grep("ko", table2$pathways),]
             sum4 <- sum(str_count(table2, "map00860"))
             name1 <- data.frame(sample_name = name,
                                 prop_of_pyroph_in_tot_pyroph_in_mag = (sum2/sum4)*100,
                                 prop_of_pyroph_in_contig = (sum2/sum1)*100,
                                 sum_of_genes_in_bin = length(unique(table2$gene_id)),
                                 sum_of_genes_in_contig = length(unique(table1$gene_id)),
                                 AFSA_gene_number = GCF_000015125.1_AFSA,
                                 gene_id_of_pyroph = c(positions)
             )
             print(name1)
             assign(name,name1, envir=.GlobalEnv)
           })
  } }
for (file in files_list){
  if(grepl("test.emapper.annotations_GCF_001509575.1", file)){
    lapply(file,
           function(x) {
             name = basename(x)
             name <- gsub("test.emapper.annotations_","",name)
             table <- read.delim(x, header=FALSE, sep="\t")
             table$V1 <- gsub("^GCF_001509575.1_", "", table$V1) # remove part of a string form a column in df
             # Step1 : Adjust the table according to the contig annotation list
             table1 = table[grepl(GCF_001509575.1_list, table$V1),]
             # Step2 : grep only the KEGG pathways (map...) from the contig genes
             table1 <- data.frame(gene_id=table1$V1, pathways=table1$V13)
             table1$pathways <- strsplit(as.character(table1$pathways), split = ",",fixed=FALSE) #split a column based on commas
             table1 <- tidyr::unnest(table1, cols = pathways) #unlist a list of strings in a column made of a list
             table1 <- table1[- grep("ko", table1$pathways),] # remove the ko pathways
             genes <- table1[grep("map00860", table1$pathways),] # grabs all rows that have the pyroph genes
             positions <- data.frame(unique(genes$gene_id))
             # positions <- paste(positions$unique.genes.gene_id.,sep = ",")
             positions <- data.frame(val=paste0(positions$unique.genes.gene_id.,collapse = ', '),stringsAsFactors = F)
             # Step3 : Count the number of time perphorine metabolism occured (i.e. how many genes) on the contig
             sum1 <- length(unique(table1$gene_id))
             sum2 <- sum(str_count(table1, "map00860"))
             # STep4 : Get the proportion of the Porphyrine genes in the whole bin
             table2 <- data.frame(gene_id=table$V1, pathways=table$V13)
             table2$pathways <- strsplit(as.character(table2$pathways), split = ",",fixed=FALSE) #split a column based on commas
             table2 <- tidyr::unnest(table2, cols = pathways) #unlist a list of strings in a column made of a list
             table2 <- table2[- grep("ko", table2$pathways),]
             sum4 <- sum(str_count(table2, "map00860"))
             name1 <- data.frame(sample_name = name,
                                 prop_of_pyroph_in_tot_pyroph_in_mag = (sum2/sum4)*100,
                                 prop_of_pyroph_in_contig = (sum2/sum1)*100,
                                 sum_of_genes_in_bin = length(unique(table2$gene_id)),
                                 sum_of_genes_in_contig = length(unique(table1$gene_id)),
                                 AFSA_gene_number = GCF_001509575.1_AFSA,
                                 gene_id_of_pyroph = c(positions)
             )
             print(name1)
             assign(name,name1, envir=.GlobalEnv)
           })
  } }
for (file in files_list){
  if(grepl("test.emapper.annotations_GCA_013335765.1", file)){
    lapply(file,
           function(x) {
             name = basename(x)
             name <- gsub("test.emapper.annotations_","",name)
             table <- read.delim(x, header=FALSE, sep="\t")
             table$V1 <- gsub("^GCA_013335765.1_", "", table$V1) # remove part of a string form a column in df
             # Step1 : Adjust the table according to the contig annotation list
             table1 = table[grepl(GCA_013335765.1_list, table$V1),]
             # Step2 : grep only the KEGG pathways (map...) from the contig genes
             table1 <- data.frame(gene_id=table1$V1, pathways=table1$V13)
             table1$pathways <- strsplit(as.character(table1$pathways), split = ",",fixed=FALSE) #split a column based on commas
             table1 <- tidyr::unnest(table1, cols = pathways) #unlist a list of strings in a column made of a list
             table1 <- table1[- grep("ko", table1$pathways),] # remove the ko pathways
             genes <- table1[grep("map00860", table1$pathways),] # grabs all rows that have the pyroph genes
             positions <- data.frame(unique(genes$gene_id))
             # positions <- paste(positions$unique.genes.gene_id.,sep = ",")
             positions <- data.frame(val=paste0(positions$unique.genes.gene_id.,collapse = ', '),stringsAsFactors = F)
             # Step3 : Count the number of time perphorine metabolism occured (i.e. how many genes) on the contig
             sum1 <- length(unique(table1$gene_id))
             sum2 <- sum(str_count(table1, "map00860"))
             # STep4 : Get the proportion of the Porphyrine genes in the whole bin
             table2 <- data.frame(gene_id=table$V1, pathways=table$V13)
             table2$pathways <- strsplit(as.character(table2$pathways), split = ",",fixed=FALSE) #split a column based on commas
             table2 <- tidyr::unnest(table2, cols = pathways) #unlist a list of strings in a column made of a list
             table2 <- table2[- grep("ko", table2$pathways),]
             sum4 <- sum(str_count(table2, "map00860"))
             name1 <- data.frame(sample_name = name,
                                 prop_of_pyroph_in_tot_pyroph_in_mag = (sum2/sum4)*100,
                                 prop_of_pyroph_in_contig = (sum2/sum1)*100,
                                 sum_of_genes_in_bin = length(unique(table2$gene_id)),
                                 sum_of_genes_in_contig = length(unique(table1$gene_id)),
                                 AFSA_gene_number = GCA_013335765.1_AFSA,
                                 gene_id_of_pyroph = c(positions)
             )
             print(name1)
             assign(name,name1, envir=.GlobalEnv)
           })
  } }
for (file in files_list){
  if(grepl("test.emapper.annotations_GCA_013335335.1", file)){
    lapply(file,
           function(x) {
             name = basename(x)
             name <- gsub("test.emapper.annotations_","",name)
             table <- read.delim(x, header=FALSE, sep="\t")
             table$V1 <- gsub("^GCA_013335335.1_", "", table$V1) # remove part of a string form a column in df
             # Step1 : Adjust the table according to the contig annotation list
             table1 = table[grepl(GCA_013335335.1_list, table$V1),]
             # Step2 : grep only the KEGG pathways (map...) from the contig genes
             table1 <- data.frame(gene_id=table1$V1, pathways=table1$V13)
             table1$pathways <- strsplit(as.character(table1$pathways), split = ",",fixed=FALSE) #split a column based on commas
             table1 <- tidyr::unnest(table1, cols = pathways) #unlist a list of strings in a column made of a list
             table1 <- table1[- grep("ko", table1$pathways),] # remove the ko pathways
             genes <- table1[grep("map00860", table1$pathways),] # grabs all rows that have the pyroph genes
             positions <- data.frame(unique(genes$gene_id))
             # positions <- paste(positions$unique.genes.gene_id.,sep = ",")
             positions <- data.frame(val=paste0(positions$unique.genes.gene_id.,collapse = ', '),stringsAsFactors = F)
             # Step3 : Count the number of time perphorine metabolism occured (i.e. how many genes) on the contig
             sum1 <- length(unique(table1$gene_id))
             sum2 <- sum(str_count(table1, "map00860"))
             # STep4 : Get the proportion of the Porphyrine genes in the whole bin
             table2 <- data.frame(gene_id=table$V1, pathways=table$V13)
             table2$pathways <- strsplit(as.character(table2$pathways), split = ",",fixed=FALSE) #split a column based on commas
             table2 <- tidyr::unnest(table2, cols = pathways) #unlist a list of strings in a column made of a list
             table2 <- table2[- grep("ko", table2$pathways),]
             sum4 <- sum(str_count(table2, "map00860"))
             name1 <- data.frame(sample_name = name,
                                 prop_of_pyroph_in_tot_pyroph_in_mag = (sum2/sum4)*100,
                                 prop_of_pyroph_in_contig = (sum2/sum1)*100,
                                 sum_of_genes_in_bin = length(unique(table2$gene_id)),
                                 sum_of_genes_in_contig = length(unique(table1$gene_id)),
                                 AFSA_gene_number = GCA_013335335.1_AFSA,
                                 gene_id_of_pyroph = c(positions)
             )
             print(name1)
             assign(name,name1, envir=.GlobalEnv)
           })
  } }

# Merge into a table
merge = rbind(EMN021,PES001,RIG001,TAF017,PLV001_002,PLV001_001,`11_GOY005`,`12_GOY005`,GOY006,
              GCF_000020465.1,GCF_000015125.1,GCF_001509575.1,GCA_013335765.1,GCA_013335335.1)
write.csv(merge, "porphyrine_genes_across_genomes.csv")

# Tabulate KEGG KO gene identification for the pyrophyrin genes next to the AFSA
# Database KO read in
KO <- read.table("KO.txt", sep="\t", header=FALSE, comment.char="#",
                 na.strings=".", stringsAsFactors=FALSE,
                 quote="", fill=FALSE)
colnames(KO) <- c("KO", "label")


for (file in files_list){
  if(grepl("test.emapper.annotations_EMN021", file)){
    lapply(file,
           function(x) {
             name = basename(x)
             name <- gsub("test.emapper.annotations_","",name)
             table <- read.table("test.emapper.annotations_EMN021", sep="\t", header=FALSE, comment.char="#",
                                 na.strings=".", stringsAsFactors=FALSE,
                                 quote="", fill=FALSE) # this is to remove the EOF error/ warning that we get with read.delim
             table$V1 <- gsub("^EMN001_021_", "", table$V1) # remove part of a string form a column in df
             # Step1 : Adjust the table according to the contig annotation list
             table1 = table[grepl(EMN021_021_list, table$V1),]
             # Step2 : Grab the table with only genes identified in the prophyrin pathway
             table1$V13 <- strsplit(as.character(table1$V13), split = ",",fixed=FALSE) #split a column based on commas
             table1 <- tidyr::unnest(table1, cols = V13) #unlist a list of strings in a column made of a list
             table1 <- table1[- grep("ko", table1$V13),] # remove the ko pathways
             table1 <- table1[grep("map00860", table1$V13),] # grabs all rows that have the pyroph genes
             table1 <- data.frame(gene_id=table1$V1, KO=table1$V12)
             # Step3 : Merge the tables together by common column
             merge <-left_join(table1, KO, by = c("KO")) 
             # test <- merge(table1, KO, by = c("KO"), all.x = TRUE)
             name1 = data.frame(sample_name = name,
                                merge)
             
             print(name1)
             assign(name,name1, envir=.GlobalEnv)
           })
  } }
for (file in files_list){
  if(grepl("test.emapper.annotations_PES001", file)){
    lapply(file,
           function(x) {
             name = basename(x)
             name <- gsub("test.emapper.annotations_","",name)
             table <- read.table("test.emapper.annotations_PES001", sep="\t", header=FALSE, comment.char="#",
                                 na.strings=".", stringsAsFactors=FALSE,
                                 quote="", fill=FALSE) # this is to remove the EOF error/ warning that we get with read.delim
             table$V1 <- gsub("^PES001_018_", "", table$V1) # remove part of a string form a column in df
             # Step1 : Adjust the table according to the contig annotation list
             table1 = table[grepl(PES001_018_list, table$V1),]
             # Step2 : Grab the table with only genes identified in the prophyrin pathway
             table1$V13 <- strsplit(as.character(table1$V13), split = ",",fixed=FALSE) #split a column based on commas
             table1 <- tidyr::unnest(table1, cols = V13) #unlist a list of strings in a column made of a list
             table1 <- table1[- grep("ko", table1$V13),] # remove the ko pathways
             table1 <- table1[grep("map00860", table1$V13),] # grabs all rows that have the pyroph genes
             table1 <- data.frame(gene_id=table1$V1, KO=table1$V12)
             # Step3 : Merge the tables together by common column
             merge <-left_join(table1, KO, by = c("KO")) 
             # test <- merge(table1, KO, by = c("KO"), all.x = TRUE)
             name1 = data.frame(sample_name = name,
                                merge)
             
             print(name1)
             assign(name,name1, envir=.GlobalEnv)
           })
  } }
for (file in files_list){
  if(grepl("test.emapper.annotations_GOY006", file)){
    lapply(file,
           function(x) {
             name = basename(x)
             name <- gsub("test.emapper.annotations_","",name)
             table <- read.table("test.emapper.annotations_GOY006", sep="\t", header=FALSE, comment.char="#",
                                 na.strings=".", stringsAsFactors=FALSE,
                                 quote="", fill=FALSE) # this is to remove the EOF error/ warning that we get with read.delim
             table$V1 <- gsub("^GOY006_RA_", "", table$V1) # remove part of a string form a column in df
             # Step1 : Adjust the table according to the contig annotation list
             table1 = table[grepl(GOY006_RA_list, table$V1),]
             # Step2 : Grab the table with only genes identified in the prophyrin pathway
             table1$V13 <- strsplit(as.character(table1$V13), split = ",",fixed=FALSE) #split a column based on commas
             table1 <- tidyr::unnest(table1, cols = V13) #unlist a list of strings in a column made of a list
             table1 <- table1[- grep("ko", table1$V13),] # remove the ko pathways
             table1 <- table1[grep("map00860", table1$V13),] # grabs all rows that have the pyroph genes
             table1 <- data.frame(gene_id=table1$V1, KO=table1$V12)
             # Step3 : Merge the tables together by common column
             merge <-left_join(table1, KO, by = c("KO")) 
             # test <- merge(table1, KO, by = c("KO"), all.x = TRUE)
             name1 = data.frame(sample_name = name,
                                merge)
             
             print(name1)
             assign(name,name1, envir=.GlobalEnv)
           })
  } }
for (file in files_list){
  if(grepl("test.emapper.annotations_TAF017", file)){
    lapply(file,
           function(x) {
             name = basename(x)
             name <- gsub("test.emapper.annotations_","",name)
             table <- read.table("test.emapper.annotations_TAF017", sep="\t", header=FALSE, comment.char="#",
                                 na.strings=".", stringsAsFactors=FALSE,
                                 quote="", fill=FALSE) # this is to remove the EOF error/ warning that we get with read.delim
             table$V1 <- gsub("^TAF017_RA_", "", table$V1) # remove part of a string form a column in df
             # Step1 : Adjust the table according to the contig annotation list
             table1 = table[grepl(TAF017_RA_list, table$V1),]
             # Step2 : Grab the table with only genes identified in the prophyrin pathway
             table1$V13 <- strsplit(as.character(table1$V13), split = ",",fixed=FALSE) #split a column based on commas
             table1 <- tidyr::unnest(table1, cols = V13) #unlist a list of strings in a column made of a list
             table1 <- table1[- grep("ko", table1$V13),] # remove the ko pathways
             table1 <- table1[grep("map00860", table1$V13),] # grabs all rows that have the pyroph genes
             table1 <- data.frame(gene_id=table1$V1, KO=table1$V12)
             # Step3 : Merge the tables together by common column
             merge <-left_join(table1, KO, by = c("KO")) 
             # test <- merge(table1, KO, by = c("KO"), all.x = TRUE)
             name1 = data.frame(sample_name = name,
                                merge)
             
             print(name1)
             assign(name,name1, envir=.GlobalEnv)
           })
  } }
for (file in files_list){
  if(grepl("test.emapper.annotations_RIG001", file)){
    lapply(file,
           function(x) {
             name = basename(x)
             name <- gsub("test.emapper.annotations_","",name)
             table <- read.table("test.emapper.annotations_RIG001", sep="\t", header=FALSE, comment.char="#",
                                 na.strings=".", stringsAsFactors=FALSE,
                                 quote="", fill=FALSE) # this is to remove the EOF error/ warning that we get with read.delim
             table$V1 <- gsub("^RIG001_014_", "", table$V1) # remove part of a string form a column in df
             # Step1 : Adjust the table according to the contig annotation list
             table1 = table[grepl(RIG001_014_list, table$V1),]
             # Step2 : Grab the table with only genes identified in the prophyrin pathway
             table1$V13 <- strsplit(as.character(table1$V13), split = ",",fixed=FALSE) #split a column based on commas
             table1 <- tidyr::unnest(table1, cols = V13) #unlist a list of strings in a column made of a list
             table1 <- table1[- grep("ko", table1$V13),] # remove the ko pathways
             table1 <- table1[grep("map00860", table1$V13),] # grabs all rows that have the pyroph genes
             table1 <- data.frame(gene_id=table1$V1, KO=table1$V12)
             # Step3 : Merge the tables together by common column
             merge <-left_join(table1, KO, by = c("KO")) 
             # test <- merge(table1, KO, by = c("KO"), all.x = TRUE)
             name1 = data.frame(sample_name = name,
                                merge)
             
             print(name1)
             assign(name,name1, envir=.GlobalEnv)
           })
  } }
for (file in files_list){
  if(grepl("test.emapper.annotations_PLV001_002", file)){
    lapply(file,
           function(x) {
             name = basename(x)
             name <- gsub("test.emapper.annotations_","",name)
             table <- read.table("test.emapper.annotations_PLV001_002", sep="\t", header=FALSE, comment.char="#",
                                 na.strings=".", stringsAsFactors=FALSE,
                                 quote="", fill=FALSE) # this is to remove the EOF error/ warning that we get with read.delim
             table$V1 <- gsub("^PLV001_002_", "", table$V1) # remove part of a string form a column in df
             # Step1 : Adjust the table according to the contig annotation list
             table1 = table[grepl(PLV001_002_list, table$V1),]
             # Step2 : Grab the table with only genes identified in the prophyrin pathway
             table1$V13 <- strsplit(as.character(table1$V13), split = ",",fixed=FALSE) #split a column based on commas
             table1 <- tidyr::unnest(table1, cols = V13) #unlist a list of strings in a column made of a list
             table1 <- table1[- grep("ko", table1$V13),] # remove the ko pathways
             table1 <- table1[grep("map00860", table1$V13),] # grabs all rows that have the pyroph genes
             table1 <- data.frame(gene_id=table1$V1, KO=table1$V12)
             # Step3 : Merge the tables together by common column
             merge <-left_join(table1, KO, by = c("KO")) 
             # test <- merge(table1, KO, by = c("KO"), all.x = TRUE)
             name1 = data.frame(sample_name = name,
                                merge)
             
             print(name1)
             assign(name,name1, envir=.GlobalEnv)
           })
  } }
for (file in files_list){
  if(grepl("test.emapper.annotations_PLV001_001", file)){
    lapply(file,
           function(x) {
             name = basename(x)
             name <- gsub("test.emapper.annotations_","",name)
             table <- read.table("test.emapper.annotations_PLV001_001", sep="\t", header=FALSE, comment.char="#",
                                 na.strings=".", stringsAsFactors=FALSE,
                                 quote="", fill=FALSE) # this is to remove the EOF error/ warning that we get with read.delim
             table$V1 <- gsub("^PLV001_001_", "", table$V1) # remove part of a string form a column in df
             # Step1 : Adjust the table according to the contig annotation list
             table1 = table[grepl(PLV001_001_list, table$V1),]
             # Step2 : Grab the table with only genes identified in the prophyrin pathway
             table1$V13 <- strsplit(as.character(table1$V13), split = ",",fixed=FALSE) #split a column based on commas
             table1 <- tidyr::unnest(table1, cols = V13) #unlist a list of strings in a column made of a list
             table1 <- table1[- grep("ko", table1$V13),] # remove the ko pathways
             table1 <- table1[grep("map00860", table1$V13),] # grabs all rows that have the pyroph genes
             table1 <- data.frame(gene_id=table1$V1, KO=table1$V12)
             # Step3 : Merge the tables together by common column
             merge <-left_join(table1, KO, by = c("KO")) 
             # test <- merge(table1, KO, by = c("KO"), all.x = TRUE)
             name1 = data.frame(sample_name = name,
                                merge)
             
             print(name1)
             assign(name,name1, envir=.GlobalEnv)
           })
  } }
for (file in files_list){
  if(grepl("test.emapper.annotations_GOY005", file)){
    lapply(file,
           function(x) {
             name = basename(x)
             name <- gsub("test.emapper.annotations_","11_",name)
             table <- read.table("test.emapper.annotations_GOY005", sep="\t", header=FALSE, comment.char="#",
                                 na.strings=".", stringsAsFactors=FALSE,
                                 quote="", fill=FALSE) # this is to remove the EOF error/ warning that we get with read.delim
             table$V1 <- gsub("^GOY005_", "", table$V1) # remove part of a string form a column in df
             # Step1 : Adjust the table according to the contig annotation list
             table1 = table[grepl(GOY005_001_11_list, table$V1),]
             # Step2 : Grab the table with only genes identified in the prophyrin pathway
             table1$V13 <- strsplit(as.character(table1$V13), split = ",",fixed=FALSE) #split a column based on commas
             table1 <- tidyr::unnest(table1, cols = V13) #unlist a list of strings in a column made of a list
             table1 <- table1[- grep("ko", table1$V13),] # remove the ko pathways
             table1 <- table1[grep("map00860", table1$V13),] # grabs all rows that have the pyroph genes
             table1 <- data.frame(gene_id=table1$V1, KO=table1$V12)
             # Step3 : Merge the tables together by common column
             merge <-left_join(table1, KO, by = c("KO")) 
             # test <- merge(table1, KO, by = c("KO"), all.x = TRUE)
             name1 = data.frame(sample_name = name,
                                merge)
             
             print(name1)
             assign(name,name1, envir=.GlobalEnv)
           })
  } }
for (file in files_list){
  if(grepl("test.emapper.annotations_GOY005", file)){
    lapply(file,
           function(x) {
             name = basename(x)
             name <- gsub("test.emapper.annotations_","12_",name)
             table <- read.table("test.emapper.annotations_GOY005", sep="\t", header=FALSE, comment.char="#",
                                 na.strings=".", stringsAsFactors=FALSE,
                                 quote="", fill=FALSE) # this is to remove the EOF error/ warning that we get with read.delim
             table$V1 <- gsub("^GOY005_", "", table$V1) # remove part of a string form a column in df
             # Step1 : Adjust the table according to the contig annotation list
             table1 = table[grepl(GOY005_001_12_list, table$V1),]
             # Step2 : Grab the table with only genes identified in the prophyrin pathway
             table1$V13 <- strsplit(as.character(table1$V13), split = ",",fixed=FALSE) #split a column based on commas
             table1 <- tidyr::unnest(table1, cols = V13) #unlist a list of strings in a column made of a list
             table1 <- table1[- grep("ko", table1$V13),] # remove the ko pathways
             table1 <- table1[grep("map00860", table1$V13),] # grabs all rows that have the pyroph genes
             table1 <- data.frame(gene_id=table1$V1, KO=table1$V12)
             # Step3 : Merge the tables together by common column
             merge <-left_join(table1, KO, by = c("KO")) 
             # test <- merge(table1, KO, by = c("KO"), all.x = TRUE)
             name1 = data.frame(sample_name = name,
                                merge)
             
             print(name1)
             assign(name,name1, envir=.GlobalEnv)
           })
  } }
for (file in files_list){
  if(grepl("test.emapper.annotations_GCF_000020465.1", file)){
    lapply(file,
           function(x) {
             name = basename(x)
             name <- gsub("test.emapper.annotations_","",name)
             table <- read.table("test.emapper.annotations_GCF_000020465.1", sep="\t", header=FALSE, comment.char="#",
                                 na.strings=".", stringsAsFactors=FALSE,
                                 quote="", fill=FALSE) # this is to remove the EOF error/ warning that we get with read.delim
             table$V1 <- gsub("^GCF_000020465.1_", "", table$V1) # remove part of a string form a column in df
             # Step1 : Adjust the table according to the contig annotation list
             table1 = table[grepl(GCF_000020465.1_list, table$V1),]
             # Step2 : Grab the table with only genes identified in the prophyrin pathway
             table1$V13 <- strsplit(as.character(table1$V13), split = ",",fixed=FALSE) #split a column based on commas
             table1 <- tidyr::unnest(table1, cols = V13) #unlist a list of strings in a column made of a list
             table1 <- table1[- grep("ko", table1$V13),] # remove the ko pathways
             table1 <- table1[grep("map00860", table1$V13),] # grabs all rows that have the pyroph genes
             table1 <- data.frame(gene_id=table1$V1, KO=table1$V12)
             # Step3 : Merge the tables together by common column
             merge <-left_join(table1, KO, by = c("KO")) 
             # test <- merge(table1, KO, by = c("KO"), all.x = TRUE)
             name1 = data.frame(sample_name = name,
                                merge)
             
             print(name1)
             assign(name,name1, envir=.GlobalEnv)
           })
  } }
for (file in files_list){
  if(grepl("test.emapper.annotations_GCF_000015125.1", file)){
    lapply(file,
           function(x) {
             name = basename(x)
             name <- gsub("test.emapper.annotations_","",name)
             table <- read.table("test.emapper.annotations_GCF_000015125.1", sep="\t", header=FALSE, comment.char="#",
                                 na.strings=".", stringsAsFactors=FALSE,
                                 quote="", fill=FALSE) # this is to remove the EOF error/ warning that we get with read.delim
             table$V1 <- gsub("^GCF_000015125.1_", "", table$V1) # remove part of a string form a column in df
             # Step1 : Adjust the table according to the contig annotation list
             table1 = table[grepl(GCF_000015125.1_list, table$V1),]
             # Step2 : Grab the table with only genes identified in the prophyrin pathway
             table1$V13 <- strsplit(as.character(table1$V13), split = ",",fixed=FALSE) #split a column based on commas
             table1 <- tidyr::unnest(table1, cols = V13) #unlist a list of strings in a column made of a list
             table1 <- table1[- grep("ko", table1$V13),] # remove the ko pathways
             table1 <- table1[grep("map00860", table1$V13),] # grabs all rows that have the pyroph genes
             table1 <- data.frame(gene_id=table1$V1, KO=table1$V12)
             # Step3 : Merge the tables together by common column
             merge <-left_join(table1, KO, by = c("KO")) 
             # test <- merge(table1, KO, by = c("KO"), all.x = TRUE)
             name1 = data.frame(sample_name = name,
                                merge)
             
             print(name1)
             assign(name,name1, envir=.GlobalEnv)
           })
  } }
for (file in files_list){
  if(grepl("test.emapper.annotations_GCF_001509575.1", file)){
    lapply(file,
           function(x) {
             name = basename(x)
             name <- gsub("test.emapper.annotations_","",name)
             table <- read.table("test.emapper.annotations_GCF_001509575.1", sep="\t", header=FALSE, comment.char="#",
                                 na.strings=".", stringsAsFactors=FALSE,
                                 quote="", fill=FALSE) # this is to remove the EOF error/ warning that we get with read.delim
             table$V1 <- gsub("^GCF_001509575.1_", "", table$V1) # remove part of a string form a column in df
             # Step1 : Adjust the table according to the contig annotation list
             table1 = table[grepl(GCF_001509575.1_list, table$V1),]
             # Step2 : Grab the table with only genes identified in the prophyrin pathway
             table1$V13 <- strsplit(as.character(table1$V13), split = ",",fixed=FALSE) #split a column based on commas
             table1 <- tidyr::unnest(table1, cols = V13) #unlist a list of strings in a column made of a list
             table1 <- table1[- grep("ko", table1$V13),] # remove the ko pathways
             table1 <- table1[grep("map00860", table1$V13),] # grabs all rows that have the pyroph genes
             table1 <- data.frame(gene_id=table1$V1, KO=table1$V12)
             # Step3 : Merge the tables together by common column
             merge <-left_join(table1, KO, by = c("KO")) 
             # test <- merge(table1, KO, by = c("KO"), all.x = TRUE)
             name1 = data.frame(sample_name = name,
                                merge)
             
             print(name1)
             assign(name,name1, envir=.GlobalEnv)
           })
  } }
for (file in files_list){
  if(grepl("test.emapper.annotations_GCA_013335765.1", file)){
    lapply(file,
           function(x) {
             name = basename(x)
             name <- gsub("test.emapper.annotations_","",name)
             table <- read.table("test.emapper.annotations_GCA_013335765.1", sep="\t", header=FALSE, comment.char="#",
                                 na.strings=".", stringsAsFactors=FALSE,
                                 quote="", fill=FALSE) # this is to remove the EOF error/ warning that we get with read.delim
             table$V1 <- gsub("^GCA_013335765.1_", "", table$V1) # remove part of a string form a column in df
             # Step1 : Adjust the table according to the contig annotation list
             table1 = table[grepl(GCA_013335765.1_list, table$V1),]
             # Step2 : Grab the table with only genes identified in the prophyrin pathway
             table1$V13 <- strsplit(as.character(table1$V13), split = ",",fixed=FALSE) #split a column based on commas
             table1 <- tidyr::unnest(table1, cols = V13) #unlist a list of strings in a column made of a list
             table1 <- table1[- grep("ko", table1$V13),] # remove the ko pathways
             table1 <- table1[grep("map00860", table1$V13),] # grabs all rows that have the pyroph genes
             table1 <- data.frame(gene_id=table1$V1, KO=table1$V12)
             # Step3 : Merge the tables together by common column
             merge <-left_join(table1, KO, by = c("KO")) 
             # test <- merge(table1, KO, by = c("KO"), all.x = TRUE)
             name1 = data.frame(sample_name = name,
                                merge)
             
             print(name1)
             assign(name,name1, envir=.GlobalEnv)
           })
  } }
for (file in files_list){
  if(grepl("test.emapper.annotations_GCA_013335335.1", file)){
    lapply(file,
           function(x) {
             name = basename(x)
             name <- gsub("test.emapper.annotations_","",name)
             table <- read.table("test.emapper.annotations_GCA_013335335.1", sep="\t", header=FALSE, comment.char="#",
                                 na.strings=".", stringsAsFactors=FALSE,
                                 quote="", fill=FALSE) # this is to remove the EOF error/ warning that we get with read.delim
             table$V1 <- gsub("^GCA_013335335.1_", "", table$V1) # remove part of a string form a column in df
             # Step1 : Adjust the table according to the contig annotation list
             table1 = table[grepl(GCA_013335335.1_list, table$V1),]
             # Step2 : Grab the table with only genes identified in the prophyrin pathway
             table1$V13 <- strsplit(as.character(table1$V13), split = ",",fixed=FALSE) #split a column based on commas
             table1 <- tidyr::unnest(table1, cols = V13) #unlist a list of strings in a column made of a list
             table1 <- table1[- grep("ko", table1$V13),] # remove the ko pathways
             table1 <- table1[grep("map00860", table1$V13),] # grabs all rows that have the pyroph genes
             table1 <- data.frame(gene_id=table1$V1, KO=table1$V12)
             # Step3 : Merge the tables together by common column
             merge <-left_join(table1, KO, by = c("KO")) 
             # test <- merge(table1, KO, by = c("KO"), all.x = TRUE)
             name1 = data.frame(sample_name = name,
                                merge)
             
             print(name1)
             assign(name,name1, envir=.GlobalEnv)
           })
  } }

# Merge into a table
merge = rbind(EMN021,PES001,RIG001,TAF017,PLV001_002,PLV001_001,`11_GOY005`,`12_GOY005`,GOY006,
              GCF_000020465.1,GCF_000015125.1,GCF_001509575.1,GCA_013335765.1,GCA_013335335.1)
write.csv(merge, "porphyrine_genes_across_genomes_with_labels.csv")

######################################################
# Heatmap for the COG annotations - COG pathways - FigureS ####

files <- list.files(path="F:/review_folder/photosynthesis/eggnog",
                    pattern="test.emapper*", full.names=TRUE, recursive=FALSE)

# Convert to long format 
lapply(files, 
       function(x) {
         # table <- read.table(x, sep = "\t",header=FALSE ,
         #                     na.strings ="", stringsAsFactors= F ) 
         table <- read_excel(x)
         # remove the first two rows and assign the first as header and rename the column to locus tag
         table <- table[-c(1),]
         names(table) <- table[1,]
         table <- table[-1,]
         # retrieve the filename basename and remove the extension
         name = basename(x)
         name <- gsub("test.emapper.annotations_","",name)
         name <- gsub(".xlsx","",name)
         name
      
         name_table = data.frame(name = table$COG_category)
         name_table = data.frame(name = table(name_table))
         name_table = data.frame(name_table, name=c(rep(name,nrow(name_table))))
         colnames(name_table) = c("abbrev", "freq", "name")
         
         name_table_gather <- name_table %>%
           gather(-abbrev, -name, key = abb, value = fre) %>% #join by (-), 
           group_by(abbrev) %>%
           # drop_na(otunum) %>%
           ungroup()
         
         assign(name, name_table_gather, envir=.GlobalEnv) #write dataframes to workspace in r
       } )

# Merge the tables 
m1 = merge(x = EMN001_021, y =PES001_018, all = TRUE)
m2 = merge(x = RIG001_014, y =TAF017_RA, all = TRUE)
m3 = merge(x = PLV001_002, y =PLV001_001, all = TRUE)
m4 = merge(x = GOY005_001, y =GOY006_RA, all = TRUE)

m5 = merge(x = m1, y = m2, all = TRUE)
m6 = merge(x = m3, y = m4, all = TRUE)

m13 = merge(x = m5, y = m6, all = TRUE)

m8 = merge(x = GCA_013335335.1, y =GCA_013335765.1, all = TRUE)
m9 = merge(x = GCA_019163275.1, y =GCF_000015125.1, all = TRUE)
m10 = merge(x = GCF_000020465.1, y =GCF_001509575.1, all = TRUE)

m11 = merge(x = m8, y =m9, all = TRUE)
m12 = merge(x = m11, y =m10, all = TRUE)

m7 = merge(x = m12, y =m13, all = TRUE)

m7$abbrev <- factor(m7$abbrev, levels=unique(m7$abbrev)) 

# Reorder x-axis
name_order <- c('GOY006_RA',
                'GOY005_001', 
                'PES001_018', 
                'RIG001_014', 
                'EMN001_021', 
                'PLV001_001', 
                'TAF017_RA', 
                'PLV001_002',
                'GCF_000020465.1',
                'GCA_013335765.1',
                'GCA_013335335.1',
                'GCF_001509575.1',
                'GCA_019163275.1',
                'GCF_000015125.1') #this vector might be useful for other plots/analyses

# Rename the column and the values in the factor : Change the labels of the legend
levels(m7$abbrev)[levels(m7$abbrev)=="A"] <- "A:RNA processing and modification"
levels(m7$abbrev)[levels(m7$abbrev)=="B"] <- "B:Chromatin structure and dynamics"
levels(m7$abbrev)[levels(m7$abbrev)=="C"] <- "C:Energy production and conversion"
levels(m7$abbrev)[levels(m7$abbrev)=="D"] <- "D:Cell cycle control, cell division, chromosome partitioning"
levels(m7$abbrev)[levels(m7$abbrev)=="E"] <- "E:Amino acid transport and metabolism"
levels(m7$abbrev)[levels(m7$abbrev)=="F"] <- "F:Nucleotide transport and metabolism"
levels(m7$abbrev)[levels(m7$abbrev)=="G"] <- "G:Carbohydrate transport and metabolism"
levels(m7$abbrev)[levels(m7$abbrev)=="H"] <- "H:Coenzyme transport and metabolism"
levels(m7$abbrev)[levels(m7$abbrev)=="I"] <- "I:Lipid transport and metabolism"
levels(m7$abbrev)[levels(m7$abbrev)=="J"] <- "J:Translation, ribosomal structure and biogenesis"
levels(m7$abbrev)[levels(m7$abbrev)=="K"] <- "K:Transcription"
levels(m7$abbrev)[levels(m7$abbrev)=="L"] <- "L:Replication, recombination and repair"
levels(m7$abbrev)[levels(m7$abbrev)=="M"] <- "M:Cel wall/membrane/envelope biogenesis"
levels(m7$abbrev)[levels(m7$abbrev)=="N"] <- "N:Cell motility"
levels(m7$abbrev)[levels(m7$abbrev)=="O"] <- "O:Posttranslational modification, protein turnover, chaperones"
levels(m7$abbrev)[levels(m7$abbrev)=="P"] <- "P:Inorganic ion transport and metabolism"
levels(m7$abbrev)[levels(m7$abbrev)=="Q"] <- "Q:Secondary metabolites biosynthesis, transport and catabolim"
levels(m7$abbrev)[levels(m7$abbrev)=="R"] <- "R:General function prediction only"
levels(m7$abbrev)[levels(m7$abbrev)=="S"] <- "S:Function unknown"
levels(m7$abbrev)[levels(m7$abbrev)=="T"] <- "T:Signal transduction mechanisms"
levels(m7$abbrev)[levels(m7$abbrev)=="U"] <- "U:Intracellular trafficking, secretion and vesicular transport"
levels(m7$abbrev)[levels(m7$abbrev)=="V"] <- "V:Defense mechanisms"
levels(m7$abbrev)[levels(m7$abbrev)=="W"] <- "W:Extracellular structures"
levels(m7$abbrev)[levels(m7$abbrev)=="X"] <- "X:Mobilome, prophages, transposons"
levels(m7$abbrev)[levels(m7$abbrev)=="Y"] <- "Y:Nuclear structure"
levels(m7$abbrev)[levels(m7$abbrev)=="Z"] <- "Z:Cytoskeleton"

#change the colours of the heatmap
# Automatic
# coul <- colorRampPalette(brewer.pal(9, "Blues"))(9)
# coul <- colorRampPalette(brewer.pal(12, "Paired"))(47)
coul <- colorRampPalette(brewer.pal(12, "Paired"))(63)

# heatmap:::
pdf("COG_category_sample_proportional_counts_heatmap.pdf",width=20,height=20)
ggplot(m7, aes(y=factor(name, level = name_order), x=factor(abbrev))) +
  geom_tile(aes(fill=fre), colour="Black", lwd =0.1, linetype=1) +
  # geom_text(aes(label = fre), color = "black", size = 1) +
  scale_fill_gradient(low = "#F497B6", high = "#0C285E") +
  guides(fill = guide_colourbar(title = "Gene proportion %",ticks = TRUE)) +
  coord_fixed() +
  xlab("COG pathways") + ylab("Genomes") +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))
  # theme_bw()
dev.off()

######################################################
