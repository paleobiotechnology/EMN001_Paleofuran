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
  mutate(frac = frac / 100,
         sample = if_else(sample %in% c("EMN001", "PES001", "PLV001"),
                          sample, "other"),
         sample = factor(sample, levels = c("EMN001", "PES001", "PLV001", "other"))) %>%
  ggplot(aes(x = 1, y = frac)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = sample), size = 2.25,
              pch = 21, width = 0.25, alpha = .8) +
  coord_flip() +
  labs(x = "",
       y = "fraction of contigs with coverage >= 5-fold",
       fill = "sample type") +
  scale_y_continuous(labels = scales::percent_format(),
                     breaks = seq(0, 1, 0.25),
                     limit = c(0, NA)) +
  scale_fill_manual(values = c("#C23616", "#0097E6", "#E1B12C", "grey50")) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top",
        legend.title = element_blank()) +
  guides(fill = guide_legend(nrow = 1))

## Panel b: fraction of contigs with q-value < 0.05 with p_damage_model <= 0.6
pyd_panelb <-  pydamage %>%
  select(sample, frac = `% of contigs with q-value < 0.05`) %>%
  mutate(frac = frac / 100,
         sample = if_else(sample %in% c("EMN001", "PES001", "PLV001"),
                          sample, "other"),
         sample = factor(sample, levels = c("EMN001", "PES001", "PLV001", "other"))) %>%
  ggplot(aes(x = 1, y = frac)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = sample), size = 2.25,
              pch = 21, width = 0.25, alpha = .8) +
  coord_flip() +
  labs(x = "",
       y = "fraction of contigs with q-value < 0.05",
       fill = "sample type") +
  scale_y_continuous(labels = scales::percent_format(),
                     breaks = seq(0, 1, 0.25),
                     limits = c(0, NA)) +
  scale_fill_manual(values = c("#C23616", "#0097E6", "#E1B12C", "grey50")) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top",
        legend.title = element_blank()) +
  guides(fill = guide_legend(nrow = 1))

## Panel c: prediction accuracy for the samples EMN001, PES001, and PLV001
predacc_df <- pydamage %>%
  select(sample, contains("<= PA")) %>%
  pivot_longer(-sample, names_to = "predacc_bin", values_to = "frac") %>%
  mutate(predacc_bin = factor(predacc_bin, levels = unique(.$predacc_bin))) %>%
  group_by(sample) %>%
  arrange(predacc_bin) %>%
  mutate(cumsum_frac = cumsum(frac)) %>%
  filter(sample %in% c("EMN001", "PLV001", "PES001")) %>%
  mutate(sample = factor(sample, levels = c("EMN001", "PES001", "PLV001"))) %>%
  ungroup()

pyd_panelc <- predacc_df %>%
  ggplot(aes(x = predacc_bin, y = cumsum_frac)) +
  geom_line(aes(group = sample), size = 0.5, colour = "grey20", lty = 3) +
  geom_point(aes(fill = sample, group = sample), size = 2, pch = 21, alpha = 0.8) +
  labs(x = "predicted accuracy (PA) bins",
       y = "cumulative fraction of contigs") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("#C23616", "#0097E6", "#E1B12C")) +
  theme(legend.position = "none",
        panel.grid.major.y = element_line(colour = "grey50", linetype = 3, size = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  guides(fill = "none")

# Stitch together
plt <- guide_area() + pyd_panela + pyd_panelb + pyd_panelc +
  plot_layout(width = c(0.5, 0.5), heights = c(0.1, 0.3, 0.6),
              guides = "collect", design = c("AA
                                              BC
                                              DD")) +
  plot_annotation(tag_levels = "a") & theme(legend.position = "top")

# Save
ggsave(snakemake@output[["pdf"]],
       plot = plt, width = 160, height = 120, units = "mm", useDingbats = F)
ggsave(snakemake@output[["png"]],
       plot = plt, width = 160, height = 120, units = "mm", dpi = 300)
