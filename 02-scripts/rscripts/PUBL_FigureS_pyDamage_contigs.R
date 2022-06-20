library(data.table)
library(tidyverse)
library(patchwork)
library(ggrepel)

theme_set(theme_classic(base_size = 9))

# Load data
sample_info <- readxl::read_excel(snakemake@params[["dataset_s1"]],
                                  sheet = 1, skip = 2)
pydamage <- readxl::read_excel(snakemake@params[["dataset_s2"]],
                                sheet = 3, skip = 2)

# Plot
## Panel a: fraction of contigs with coverage >= 5-fold
pyd_panela <- pydamage %>%
  select(sample, frac = `% of contigs with coverage >= 5x`) %>%
  mutate(frac = frac / 100) %>%
  left_join(sample_info %>%
            mutate(sample = str_sub(sampleId, 1, 6)) %>%
            select(sample, sampletype = `Common name`) %>%
            distinct(),
            by = "sample") %>%
  ggplot(aes(x = 1, y = frac)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = sampletype), size = 2.25,
              pch = 21, width = 0.25, alpha = .8) +
  coord_flip() +
  labs(x = "",
       y = "fraction of contigs with coverage >= 5-fold",
       fill = "sample type") +
  scale_y_continuous(labels = scales::percent_format(),
                     breaks = seq(0, 1, 0.25),
                     limit = c(0, NA)) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top")

## Panel b: fraction of contigs with q-value < 0.05 with p_damage_model <= 0.6
pyd_panelb <-  pydamage %>%
  select(sample, frac = `% of contigs with q-value < 0.05`) %>%
  mutate(frac = frac / 100) %>%
  left_join(sample_info %>%
            mutate(sample = str_sub(sampleId, 1, 6)) %>%
            select(sample, sampletype = `Common name`),
            by = "sample") %>%
  ggplot(aes(x = 1, y = frac)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = sampletype), size = 2.25,
              pch = 21, width = 0.25, alpha = .8) +
  coord_flip() +
  labs(x = "",
       y = "fraction of contigs with q-value < 0.05",
       fill = "sample type") +
  scale_y_continuous(labels = scales::percent_format(),
                     breaks = seq(0, 1, 0.25),
                     limits = c(0, NA)) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top")

## Panel c: prediction accuracy for the samples EMN001, PES001, and PLV001
predacc_df <- pydamage %>%
  select(sample, contains("<= PA")) %>%
  pivot_longer(-sample, names_to = "predacc_bin", values_to = "frac") %>%
  mutate(predacc_bin = factor(predacc_bin, levels = unique(.$predacc_bin))) %>%
  group_by(sample) %>%
  arrange(predacc_bin) %>%
  mutate(cumsum_frac = cumsum(frac)) %>%
  filter(sample %in% c("EMN001", "PLV001", "PES001")) %>%
  left_join(sample_info %>%
              mutate(sample = str_sub(sampleId, 1, 6)) %>%
              select(sample, sampletype = `Common name`),
              by = "sample") %>%
  ungroup()

pyd_panelc <- predacc_df %>%
  ggplot(aes(x = predacc_bin, y = cumsum_frac)) +
  geom_line(aes(group = sample), size = 0.5, colour = "grey20", lty = 3) +
  geom_point(aes(fill = sampletype, group = sample), size = 2.5, pch = 21) +
  geom_text_repel(data = filter(predacc_df,
                                predacc_bin == "90% <= PA < 95%"),
                  aes(label = sample),
                  size = 2.5,
                  nudge_y = 0.3,
                  xlim = c(NA, 14),
                  ylim = c(0.25, NA)) +
  labs(x = "predicted accuracy (PA) bins",
       y = "cumulative fraction of contigs") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Stitch together
plt <- guide_area() + pyd_panela + pyd_panelb + pyd_panelc +
  plot_layout(width = c(0.5, 0.5), heights = c(0.1, 0.2, 0.7),
              guides = "collect", design = c("AA
                                              BC
                                              DD")) +
  plot_annotation(tag_levels = "a")

# Save
ggsave(snakemake@output[["pdf"]],
       plot = plt, width = 160, height = 120, units = "mm", useDingbats = F)
ggsave(snakemake@output[["png"]],
       plot = plt, width = 160, height = 120, units = "mm", dpi = 300)
