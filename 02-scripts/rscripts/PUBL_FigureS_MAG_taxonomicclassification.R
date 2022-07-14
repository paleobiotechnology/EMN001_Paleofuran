library(tidyverse)

theme_set(theme_classic(base_size = 8))

# Auxilliary functions
## https://github.com/rmcelreath/rethinking/blob/master/R/colors.r
col.desat <- function( acol , amt=0.3 ) {
    acol <- col2rgb(acol)
    ahsv <- rgb2hsv(acol)
    ahsv[2] <- ahsv[2] * amt
    hsv( ahsv[1] , ahsv[2] , ahsv[3] )
}

## https://stackoverflow.com/a/8197703
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Set colour scheme
col_mags <- c(gg_color_hue(2)[1], col.desat(gg_color_hue(2)[1]),
              gg_color_hue(2)[2], col.desat(gg_color_hue(2)[2]))

# Load data
nonred_mags <- readxl::read_excel(snakemake@params[["dataset_s3"]],
                                  sheet = 5, skip = 2)
mag_clusters <- readxl::read_excel(snakemake@params[["dataset_s3"]],
                                  sheet = 6, skip = 2)

# Prepare data
taxon_assignment <- nonred_mags %>%
  select(`binID`, `primary cluster`, `secondary cluster`, taxon = `GTDB classification`, oral = `oral taxon`) %>%
  left_join(mag_clusters, by = c("binID" = "representative MAG")) %>%
  mutate(taxon = str_match(taxon, ".+;(o__.+)")[,2],
         family = str_match(taxon, ".+;(f__.+);g__.*")[,2],
         genus = str_match(taxon, ".+;(g__.+);s__.*")[,2],
         species = str_match(taxon, ".+;(s__.+)")[,2])

taxon_assignment_summary <- taxon_assignment %>%
  select(`secondary cluster`, family, genus, oral, contains("Q")) %>%
  group_by(genus) %>%
  summarise(across(contains("Q"), sum),
            oral = sum(oral != "none"),
            n_cluster = n()) %>%
  mutate(total_ancient = `ancient - HQ` + `ancient - MQ`,
         total_modern = `modern - HQ` + `modern - MQ`,
         total = total_ancient + total_modern)

# Plot
plt <- taxon_assignment_summary %>%
  mutate(label = if_else(is.na(genus), str_c("Unassigned: ", n_cluster),
                         str_c(str_replace(genus, "g__", ""), ": ", n_cluster))) %>%
  arrange(desc(`ancient - HQ`), desc(`modern - HQ`)) %>%
  mutate(label = factor(label, levels = .$label),
         oral = if_else(oral > 0, "oral", "other")) %>%
  arrange(label, oral) %>%
  mutate(cumsum_samples = cumsum(oral == "oral")) %>%
  mutate(oral = case_when(
           cumsum_samples <= 26 & oral == "oral" ~ "oral - 1",
           cumsum_samples > 26 & oral == "oral" ~ "oral - 2",
           T ~ "other"
         )) %>%
  arrange(oral, desc(as.integer(label))) %>%
  mutate(label = factor(label, levels = .$label)) %>%
  select(label, contains("Q"), oral) %>%
  pivot_longer(contains("Q"), names_to = "MAG type", values_to = "n") %>%
  ggplot(aes(x = label, y = n)) +
  geom_col(aes(group = label, fill = `MAG type`)) +
  coord_flip() +
  facet_wrap(~ oral, nrow = 1, scales = "free_y",
             labeller = as_labeller(c("oral - 1" = "oral", "oral - 2" = "oral", "other" = "other"))) +
  labs(x = "genus",
       y = "number of MAGs",
       fill = "MAG type") +
  scale_fill_manual(values = col_mags) +
  theme(legend.position = "top",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid.major.x = element_line(colour = "grey50", size = 0.5, linetype = 3),
        axis.title = element_text(size = 8),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size = 6),
        strip.background = element_blank(),
        strip.text = element_text(size = 8)) +
  guides(fill = guide_legend(keywidth = unit(3, "mm"), keyheight = unit(3, "mm")))

# Save
ggsave(snakemake@output[["pdf"]],
       plot = plt, width = 160, height = 120, units = "mm", useDingbats = F)
ggsave(snakemake@output[["png"]],
       plot = plt, width = 160, height = 120, units = "mm", dpi = 300)

