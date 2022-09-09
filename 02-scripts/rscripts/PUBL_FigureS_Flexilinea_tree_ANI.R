library(tidyverse)
library(pheatmap)

heatmap_cols <- c("#F497B6","#E8749A","#D12E64", "#87153B","#5E0C27",
                  "#97E4F4","#74C9E8", "#2D92D1","#174E86","#0C285E")

flexilinea_ani <- readxl::read_excel(snakemake@params[["dataset_s8"]],
                                     sheet = 3, skip = 2)


flexilinea_order <- c("GCA_012798535.1", "GCA_012515195.1", "GCA_017432685.1",
                      "GCA_017433665.1", "GCF_001192795.1", "GCA_002305345.1", "GCA_002407845.1",
                      "GCA_002329625.1", "GCA_018056295.1", "GCA_022837205.1", "EMN001_017",
                      "PES001_009", "GOY005_018", "FUM002_008", "RIG001_011", "EMN001_010",
                      "FUM002_005", "PES001_006", "GOY005_016", "JAE008_005",
                      "JAE009_025", "VLC008_012", "GCF_001717545.1", "VLC001_002", "JAE015_012",
                      "JAE014_017", "GCA_905372085.1", "PYK005_007", "ESA006_002", "PYK004_001",
                      "JAE016_007", "MOA001_001", "VLC005_010", "ECO010_007", "TAF008_003",
                      "OAK003_004", "OFN001_004", "ECO002_009", "PES001_013", "RIG001_015",
                      "EMN001_001")

ani_matrix <- flexilinea_ani %>%
  select(-c(name, `NCBI taxonomy Id`)) %>%
  column_to_rownames(var = "sample1") %>%
  as.matrix()
ani_matrix <- ani_matrix[order(match(colnames(ani_matrix), flexilinea_order)),
                         order(match(colnames(ani_matrix), flexilinea_order))]


pdf(snakemake@output[["pdf"]], 6, 6)
pheatmap(ani_matrix, display_numbers = F,
         color = heatmap_cols,
         main = "",
         show_rownames = TRUE, show_colnames = TRUE, fontsize  = 5,
         cluster_cols = F,
         cluster_rows = F)
dev.off()
