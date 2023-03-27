library(tidyverse)
library(patchwork)

theme_set(theme_classic(base_size = 9))

# Load data
kegg_ko <- readxl::read_excel(snakemake@params[["dataset_s8"]],
                              sheet = 4, skip = 2)
photosyn <- read.table(snakemake@params[["ko_photo"]],
                       sep = ";", header=FALSE , na.strings ="", stringsAsFactors= F )


# Panel A:
panel_a <- kegg_ko %>% 
  rename(hierarchy_3rdl = `hierarchy 3rd level`) %>%
  group_by(hierarchy_3rdl,  specificity) %>%
  count() %>%
  group_by(hierarchy_3rdl) %>%
  filter(sum(n) > 1) %>%
  mutate(specificity = factor(specificity, levels = c("C. limicola", "Chlorobium MAG")),
         hierarchy_3rdl = factor(hierarchy_3rdl, levels = rev(unique(.$hierarchy_3rdl)))) %>%
  ggplot(aes(x = hierarchy_3rdl, y = specificity, fill = n)) +
  geom_tile() +
  coord_flip() +
  labs(x = "3rd-level KEGG BRITE hierarchy",
       y = "specific for",
       fill = "number of KEGG\northologs") +
  scale_fill_gradient(low = "white", high = "red") +
  theme(axis.text.y = element_text(size = 6))

# Panel B:
files <- list.files(path="tmp/eggnog",
                    pattern="test.emapper*", full.names=TRUE, recursive=FALSE)

# Convert to long format 
lapply(files, 
       function(x) {
         # table <- read.table(x, sep = "\t",header=FALSE ,
         #                     na.strings ="", stringsAsFactors= F ) 
         table <- readxl::read_excel(x)
         # remove the first two rows and assign the first as header and rename the column to locus tag
         table <- table[-c(1),]
         names(table) <- table[1,]
         table <- table[-1,]
         # retrieve the filename basename and remove the extension
         name = basename(x)
         name <- gsub("test.emapper.annotations_","",name)
         name <- gsub(".xlsx","",name)
         
         new_table <- table[ grepl(paste(photosyn$V1, collapse="|"), table$KEGG_ko),]
      
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
m7 =merge (m7, photosyn, by.x = 'ko', by.y ='V1')

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

panel_b <- ggplot(m7, aes(y=factor(name, level = name_order), x=factor(V2))) +
  geom_tile(aes(fill=freq), colour="Black", lwd =0.1, linetype=1) +
  scale_fill_gradient(low = "#F497B6", high = "#0C285E") +
  guides(fill = guide_colourbar(title = "Gene proportion %",ticks = TRUE)) +
  coord_fixed() +
  xlab("KO genes") + ylab("Genomes") +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))

# Patch
plt <- panel_a + panel_b +
  plot_layout(ncol = 1, heights = c(.75, .25)) +
  plot_annotation(tag_levels = "a")

# Save
ggsave(snakemake@output[["pdf"]],
       plot = plt, width = 160, height = 180, units = "mm", useDingbats = F)
ggsave(snakemake@output[["png"]],
       plot = plt, width = 160, height = 180, units = "mm", dpi = 300)
