library(rgeos)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(data.table)
library(tidyverse)
library(ggrepel)
library(patchwork)
library(ggh4x)
library(ggpubr)
library(scatterpie)
theme_set(theme_classic(base_size = 8))

# Load data
## Sample information with the geographic coordinates
sample_info <- readxl::read_excel(snakemake@params[["dataset_s1"]],
                                  sheet = 1, skip = 2)
## The fragment length distribution
fraglength <- readxl::read_excel(snakemake@params[["dataset_s1"]],
                                 sheet = 4, skip = 2) %>%
  filter(sample %in% c("EMN001", "PES001", "PLV001"))
## The summary of the de novo assembly by calN50
assembly <- readxl::read_excel(snakemake@params[["dataset_s2"]],
                               sheet = 1, skip = 2) %>%
  filter(sample %in% c("EMN001", "PES001", "PLV001") | str_sub(sample, 1, 3) %in% c("JAE", "VLC")) %>%
  mutate(sampletype = if_else(str_sub(sample, 1, 3) %in% c("JAE", "VLC"),
                              "modern", "ancient"))
## Overview list of the MAGs
mags <- readxl::read_excel(snakemake@params[["dataset_s3"]],
                           sheet = 3, skip = 2) %>%
  mutate(sample = str_sub(binID, 1, 6)) %>%
  filter(sample %in% c("EMN001", "PES001", "GOY005", "RIG001", "PLV001"),
         `checkM completeness [%]` >= 50,
         `checkM contamination [%]` < 10,
         `GUNC contamination [%]` < 10,
         `GUNC clade separation score` < 0.45,
         `ratio non-syn. to syn. minor alleles [%]` < 1) %>%
  mutate(sample = factor(sample, levels = c("EMN001", "PES001", "RIG001", "GOY005", "PLV001")))
nonred_mags <- readxl::read_excel(snakemake@params[["dataset_s3"]],
                                  sheet = 5, skip = 2)
## Damageprofiler overview for C. limicola and D. oralis for EMN001
damageprofiler <- readxl::read_excel(snakemake@params[["dataset_s3"]],
                                     sheet = 7, skip = 2) %>%
  filter(sample == "EMN001")

# Colour schemes
## EMN001, PES001, PLV001, and modern
preservation_col <- c("#C23616", "#0097E6", "#E1B12C")


# Panel a: maps
## Prepare the table for plotting the map
map_info <- sample_info %>%
    mutate(individualId = str_sub(`sampleId`, 1, 6),
           sampletype = if_else(`Analysis group` == "Neanderthal",
                                "Neanderthal", "modern human"),
           sampleperiod = recode(`Specimen period`,
             `18th century` = "historic",
             `19th century` = "historic",
             `1st-2nd century` = "historic",
             `Later Stone Age` = "Neolithic",
              `Middle Paleolithic` = "Paleolithic",
              `Modern` = "modern",
              `Upper Paleolithic` = "Paleolithic",
              `Upper Paleolithic (Gravettian)` = "Paleolithic",
              `Upper Paleolithic (Magdalenian)` = "Paleolithic",
           ),
           sampleperiod = factor(sampleperiod,
                                 levels = c("Paleolithic", "Mesolithic",
                                            "Neolithic", "historic", "modern"))) %>%
    select(individualId, sampletype,
           lat = Latitude, long = Longitude,
           sampleperiod)
## Plot full map including samples from Africa
world_map <- ggplot(data = ne_countries(scale = 110, returnclass = "sf")) +
  ## Plot coast line
  geom_sf(colour = "grey80",
          fill = "grey80",
          linetype = "solid") +
  ## Plot land masses with highlighting country borders
  geom_sf(data = ne_coastline(scale = 110, returnclass = "sf"),
          colour = "grey90",
          linetype = "solid",
          size = 0.3) +
  ## Restrict size of the map to Africa and Western Eurasia
  coord_sf(xlim = c(-15, 45),
           ylim = c(-35, 55),
           expand = F) +
  ## Plot the sample points
  geom_point(data = map_info,
             aes(x = long, y = lat, shape = sampletype,
                 fill = sampleperiod),
             colour = "black",
             size = 2) +
  scale_fill_manual(values = c("#B71C1C", "#EF5350", "#FFAB91",
                               "#FFD54F", "#AED581")) +
  scale_shape_manual(values = c(21, 24)) +
  ## Adjust the theme
  labs(fill = "sample period",
       shape = "sample type") +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom",
        legend.box.background = element_rect(colour = "black", size = 0.5,
                                             fill = "white"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.spacing.x = unit(0, "mm"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(c(1,1,0,0), "mm"),
        plot.background = element_rect(colour = "black", size = 0.5)) +
  guides(shape = guide_legend(nrow = 2, title.position = "top",
                              keyheight = unit(4, "mm")),
         fill = guide_legend(nrow = 2, title.position = "top",
                             keyheight = unit(4, "mm"),
                             override.aes = list(pch = 21)))
## Extract the full legend of the world map into a grob to place it later into
## the map of Europe
wm_legend <- get_legend(world_map)
## Plot map of samples with Chlorobiaceae samples
eur_map <- ggplot(data = ne_countries(scale = 110, returnclass = "sf")) +
  ## Plot coast line
  geom_sf(colour = "grey80",
          fill = "grey80",
          linetype = "solid") +
  ## Plot land masses with highlighting country borders
  geom_sf(data = ne_coastline(scale = 110, returnclass = "sf"),
          colour = "grey90",
          linetype = "solid",
          size = 0.3) +
  ## Restrict size of the map to Western Eurasia
  coord_sf(ylim = c(27, 55),
           xlim = c(-50, 30),
           expand = F) +
  ## Plot the sample points and labels
  geom_point(data = map_info %>%
                    filter(individualId %in% c("PES001", "GOY005", "GOY006",
                                               "EMN001", "PLV001", "RIG001",
                                               "TAF017")),
             aes(x = long, y = lat, shape = sampletype,
                 fill = sampleperiod),
             colour = "black",
             size = 2.5) +
  geom_label_repel(data = map_info %>%
                          filter(individualId %in% c("PES001", "GOY005", "GOY006",
                                                     "EMN001", "PLV001", "RIG001",
                                                     "TAF017")) %>%
                          distinct(),
             aes(x = long, y = lat, label = individualId),
             nudge_x = 1.5, nudge_y = -1.5,
             ylim = c(38, Inf),
             colour = "black",
             size = 2.5) +
  scale_fill_manual(values = c("#B71C1C", "#EF5350", "#FFAB91",
                               "#FFD54F", "#AED581")) +
  scale_shape_manual(values = c(21, 24)) +
  ## Adjust the theme
  labs(fill = "specimen period",
       shape = "sample type") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        plot.margin = unit(c(0,0,0,0), "mm"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
## Plot the map of the Chlorobiaceae samples and place the overview map and the
## legend as inleys
panel_a <- eur_map +
  annotation_custom(grob = ggplotGrob(world_map + theme(legend.position = "none")),
                    xmin = -37, xmax = -10,
                    ymin = -Inf, ymax = Inf) +
  annotation_custom(grob = wm_legend,
                    xmin = -10, xmax = Inf,
                    ymin = -Inf, ymax = 32)

# Panel b: The fragment length distribution of the samples EMN001, PES001, and PLV001
fraglength_lt <- select(fraglength, -c(`no. of DNA molecules`, '> 140bp')) %>%
  pivot_longer(-sample, names_to = "length", values_to = "prop") %>%
  mutate(length = as.double(str_replace(length, "bp", "")),
         sample = factor(sample, levels = c("EMN001", "PES001", "PLV001")))
panel_b <- fraglength_lt %>%
  ggplot(aes(x = length, y = prop, group = sample, colour = sample)) +
  geom_line(size = 0.8) +
  labs(x = "DNA molecule length [bp]",
       y = "fraction of molecules") +
  scale_x_continuous(breaks = c(30, 60, 90, 120, 140)) +
  scale_y_continuous(labels = scales::percent_format(),
                     limits = c(NA, 0.02)) +
  scale_colour_manual(values = preservation_col) +
  scale_fill_manual(values = preservation_col) +
  theme(legend.position = c(0.05, 0.2),
        legend.justification = "left",
        legend.margin = margin(0, 0, 0, 0),
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        panel.grid.major.y = element_line(colour = "grey50", linetype = 3, size = 0.5),
        axis.line.x = element_blank()) +
  guides(colour = guide_legend(nrow = 3, keyheight = unit(3, "mm")))

# Panel c: Nx distribution
panel_c <- assembly %>%
  select(sample, N0:N100) %>%
  pivot_longer(-sample, names_to = "decile", values_to = "length") %>%
  mutate(decile = as.double(str_replace(decile, "N", "")) - 5,
         length = length / 1e3,
         sample = if_else(sample %in% c("EMN001", "PES001", "PLV001"), sample, "modern"),
         sample = factor(sample, levels = c("modern", "EMN001", "PES001", "PLV001"))) %>%
  group_by(sample, decile) %>%
  summarise(length = mean(length)) %>%
  ggplot(aes(x = decile, y = length, group = sample, colour = sample)) +
  geom_step(aes(linetype = sample), size = 0.8) +
  labs(x = "Nx",
       y = "contig length [kb]") +
  scale_x_continuous(breaks = seq(0, 90, 10),
                     limits = c(-5, 95)) +
  scale_y_log10(limit = c(0.5, 500),
                breaks = c(1, 5, 10, 50, 100, 500),
                expand = c(0, NA),
                guide = "axis_logticks") +
  scale_colour_manual(values = c("grey50", preservation_col)) +
  scale_linetype_manual(values = c("twodash", rep("solid", 3))) +
  theme(legend.position = c(1, 0.95),
        legend.justification = "right",
        legend.margin = margin(0, 0, 0, 0),
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        axis.line.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey50", size=0.5, linetype = 3),
        axis.ticks.length.y = unit(3, "mm"),
        ggh4x.axis.ticks.length.minor = rel(0.5),
        ggh4x.axis.ticks.length.mini = rel(0.3)) +
  guides(colour = guide_legend(nrow = 2, keyheight = unit(4, "mm")))

# Panel d: Overview of the MAGs
## Prepare the data for plotting
mags_overview <- mags %>%
  mutate(representative = binID %in% nonred_mags$binID) %>%
  left_join(list(nonred_mags %>%
                 select(binID, `members of cluster`, `oral taxon`, `SGB type`) %>%
                 separate_rows(`members of cluster`, sep = ", ") %>%
                 rename(reprMAG = binID,
                        binID = `members of cluster`) %>%
                 mutate(binID = ifelse(is.na(binID), reprMAG, binID)),
                 nonred_mags %>%
                 filter(str_sub(binID, 1, 6) %in% c("EMN001", "PES001", "GOY005",
                                                    "PLV001", "RIG001")) %>%
                 select(binID, `oral taxon`, `SGB type`) %>%
                 mutate(reprMAG = binID)) %>%
            bind_rows(),
            by = "binID") %>%
  distinct() %>%
  select(binID, sample, MIMAG,
         completeness = `checkM completeness [%]`,
         representative,
         reprMAG,
         `oral taxon`,
         `SGB type`) %>%
  mutate(sample = factor(sample, levels = c("EMN001", "PES001", "GOY005", "PLV001", "RIG001")),
         chlorobium = binID %in% c("EMN001_021", "PLV001_002", "GOY005_001",
                                   "PES001_018", "PLV001_001", "RIG001_014")) %>%
  arrange(sample, desc(completeness)) %>%
  mutate(binID = factor(binID, levels = .$binID)) %>%
  group_by(sample) %>%
  mutate(r = row_number(desc(completeness))) %>%
  ungroup() %>%
  mutate(x = (r - 1) %/% 8,
         y = abs((r - 1) %% 8 - 8),
         type = case_when(chlorobium ~ "Chlorobium",
                          `oral taxon` != "none" ~ "oral taxon",
                          `SGB type` == "kSGB" ~ "known SGB",
                          `SGB type` == "uSGB" ~ "unknown SGB"),
         type = factor(type, levels = c("Chlorobium", "oral taxon", "known SGB",
                                        "unknown SGB")),
         quality = factor(if_else(MIMAG == "high", "HQ", "MQ"), levels = c("HQ", "MQ")),
         id = as.numeric(sample))

sample_names <- tibble(x = as.numeric(mags_overview$sample),
                       y = as.character(mags_overview$sample)) %>%
                distinct() %>%
                deframe()

## Convert the data to a numeric-only dataframe for scatterpie
mags_overview_numeric <- mags_overview %>%
  select(id, x, y, r, type, completeness) %>%
  mutate(type = LETTERS[as.numeric(type)],
         E = 100 - completeness) %>%
  pivot_wider(names_from = "type", values_from = "completeness", values_fill = 0)
## Plot
panel_d <- ggplot() +
  geom_scatterpie(data = mags_overview_numeric,
                  aes(x = x, y = y, group = r, r = 0.4),
                  cols = LETTERS[1:5], colour = NA) +
  geom_segment(data = mags_overview %>%
                      filter(quality == "HQ"),
               aes(x = x - 0.2, xend = x + 0.2,
                   y = y - 0.5, yend = y - 0.5),
               colour = "#9C27B0") +
  geom_hline(yintercept = 0, size = 0.5, colour = "grey50") +
  coord_equal() +
  facet_wrap(~ id,
             labeller = as_labeller(sample_names),
             nrow = 1) +
  labs(x = "samples with HQ Chlorobium MAGs",
       y = "metagenome-assembled genomes") +
  scale_fill_manual(values = c("#E57373", "#42A5F5", "#2C3E50", "#795548", "white"),
                    labels = c("Chlorobium MAG", "oral taxa", "other known MAG", "other unknown MAG", "")) +
  theme(legend.position = "top",
        legend.margin = margin(0, 0, 0, 0),
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key.height = unit(3, "mm"),
        legend.key.width = unit(3, "mm"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = 7),
        strip.background = element_blank()) +
  guides(fill = guide_legend(nrow = 2))

# Panel E
## Prepare data
damage_profile <- damageprofiler %>%
  select(-c(sample, MAG)) %>%
  rename(taxon = genus, Pos = `position from 5' end of read`) %>%
  pivot_longer(-c(taxon, Pos), names_to = "substitution", values_to = "freq") %>%
  mutate(taxon = factor(taxon, levels = c("Chlorobium", "Flexilinea")),
         substitution = case_when(
                          substitution %in% c("C>T", "G>A") ~ substitution,
                          substitution %in% c("T>C", "A>G") ~ "transitions",
                          T ~ "transversions"
                        ),
         substitution = factor(substitution, levels = c("C>T", "G>A",
                                                        "transitions", "transversions"))) %>%
  group_by(taxon, substitution, Pos) %>%
  summarise(freq = mean(freq)) %>%
  ungroup()

## Plot
panel_e <- damage_profile %>%
  ggplot(aes(x = Pos, y = freq, group = substitution)) +
  geom_line(aes(colour = substitution, alpha = substitution), size = 0.8) +
  geom_text(data = select(damage_profile, taxon, substitution) %>%
                   filter(substitution == "C>T") %>%
                   distinct(),
            aes(x = 12.5, y = 0.27, label = taxon),
            colour = "black", size = 2.2, fontface = "italic") +
  facet_wrap(~ taxon, nrow = 1) +
  labs(x = "position from the 5' end of the read",
       y = "frequency") +
  scale_y_continuous(labels = scales::percent_format(),
                     limits = c(NA, 0.38)) +
  scale_colour_manual(values = c("#F44336", "#2196F3", "#CE93D8", "#A5D6A7")) +
  scale_alpha_manual(values = c(1, 1, 0.8, 0.8)) +
  theme(legend.position = c(0.5, 0.95),
        legend.justification = "center",
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) +
  guides(colour = guide_legend(nrow = 2, keyheight = unit(3, "mm")))

# Stitch together
plt <- panel_a + panel_b + panel_c + panel_d + panel_e +
  plot_layout(design = "AAADDD
                        BBCCEE",
              heights = c(0.6, 0.4)) +
  plot_annotation(tag_levels = "a")

# Save
ggsave(snakemake@output[["pdf"]],
       plot = plt, width = 180, height = 140, units = "mm", useDingbats = F)
ggsave(snakemake@output[["png"]],
       plot = plt, width = 180, height = 140, units = "mm", dpi = 300)
