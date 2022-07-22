library(tidyverse)

theme_set(theme_classic(base_size = 8))

# Load data
damage <- readxl::read_excel(snakemake@params[["dataset_s3"]],
                             sheet = 7, skip = 2)

damage_profile <- damage %>%
  select(-MAG) %>%
  rename(Pos = `position from 5' end of read`,
         taxon = genus) %>%
  pivot_longer(-c(sample, taxon, Pos), names_to = "substitution", values_to = "freq") %>%
  mutate(taxon = factor(taxon, levels = c("Chlorobium", "Flexilinea")),
         substitution = case_when(
                          substitution %in% c("C>T", "G>A") ~ substitution,
                          substitution %in% c("T>C", "A>G") ~ "transitions",
                          T ~ "transversions"
                        ),
         substitution = factor(substitution, levels = c("C>T", "G>A",
                                                        "transitions", "transversions"))) %>%
  group_by(sample, taxon, substitution, Pos) %>%
  summarise(freq = mean(freq)) %>%
  ungroup() %>%
  mutate(sample = factor(sample, levels = c("EMN001", "PES001", "GOY005", "PLV001", "RIG001")))

## Plot
plt <- damage_profile %>%
  ggplot(aes(x = Pos, y = freq, group = substitution)) +
  geom_line(aes(colour = substitution, alpha = substitution), size = 0.8) +
  facet_grid(sample ~ taxon) +
  labs(x = "position from the 5' end of the read",
       y = "frequency") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_colour_manual(values = c("#F44336", "#2196F3", "#CE93D8", "#A5D6A7")) +
  scale_alpha_manual(values = c(1, 1, 0.8, 0.8)) +
  theme(legend.position = "top",
        legend.justification = "center",
        legend.title = element_blank(),
        panel.grid.major.y = element_line(size = 0.5, colour = "grey50", linetype = 3),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(face = "italic")) +
  guides(colour = guide_legend(nrow = 1, keyheight = unit(3, "mm")))

# Save
ggsave(snakemake@output[["pdf"]],
       plot = plt, width = 160, height = 140, units = "mm", useDingbats = F)
ggsave(snakemake@output[["png"]],
       plot = plt, width = 160, height = 140, units = "mm", dpi = 300)
