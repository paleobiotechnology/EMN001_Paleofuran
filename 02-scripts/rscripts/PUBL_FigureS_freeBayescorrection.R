library(data.table)
library(tidyverse)
library(patchwork)

# Load data
excluded_samples <- c("DLV001", "DLV002", "SPM001", "SPM002", "TAF016")
substitutions <- fread(snakemake@params[["subst"]]) %>%
  filter(!(sample %in% excluded_samples))
mafs <- fread(snakemake@params[["maf"]], header = T) %>%
  filter(!(sample %in% excluded_samples))

# Plot
## Panel a: ratio transition vs. transversions
ratio_ts_tv <- substitutions %>%
  pivot_longer(-sample, names_to = "subst", values_to = "counts") %>%
  mutate(substtype = if_else(subst %in% c("AG", "CT", "GA", "TC"),
                             "transition", "transversion")) %>%
  group_by(sample, substtype) %>%
  summarise(counts = sum(counts)) %>%
  pivot_wider(names_from = "substtype", values_from = "counts") %>%
  mutate(total = transition + transversion,
         across(transition:transversion, ~ . / total)) %>%
  arrange(transition) %>%
  mutate(sample = factor(sample, levels = .$sample)) %>%
  select(-total) %>%
  pivot_longer(-sample, names_to = "substtype", values_to = "freq") %>%
  mutate(substtype = factor(substtype, levels = c("transition", "transversion")))

panel_a <- ratio_ts_tv %>%
  ggplot(aes(x = sample, y = freq, fill = substtype)) +
  geom_col(position = position_stack(reverse = T)) +
  geom_hline(yintercept = 0.5, colour = "grey50", lty = 1, size = 0.5) +
  coord_flip() +
  labs(y = "substitution frequency",
       fill = "") +
  scale_y_continuous(expand = c(0.005, 0),
                     labels = scales::percent_format()) +
  scale_fill_manual(values = c("#CE93D8", "#A5D6A7")) +
  theme_classic(base_size = 8) +
  theme(legend.position = "top",
        legend.key.size = unit(0.6, "line"),
        legend.margin=margin(0, 0, 0, 0),
        legend.box.margin=margin(-3, -3, -3, -3),
        axis.text.y = element_text(size = 6),
        axis.line.y = element_blank()) +
  guides(fill = guide_legend(nrow = 1))

## Panel b: ratio T>C and A>G to C>T and G>A
ratio_damage_muts <- substitutions %>%
  select(sample, AG, CT, GA, TC) %>%
  mutate(`C>T & G>A` = CT + GA,
         `T>C & A>G` = AG + TC,
         total = `C>T & G>A` + `T>C & A>G`,
         across(`C>T & G>A`:`T>C & A>G`, ~ . / total)) %>%
  mutate(sample = factor(sample, levels = levels(ratio_ts_tv$sample))) %>%
  select(sample, `C>T & G>A`, `T>C & A>G`) %>%
  pivot_longer(-sample, names_to = "substtype", values_to = "freq") %>%
  mutate(substtype = factor(substtype, levels = c("C>T & G>A", "T>C & A>G")))

panel_b <- ratio_damage_muts %>%
  ggplot(aes(x = sample, y = freq, fill = substtype)) +
  geom_col(position = position_stack(reverse = T)) +
  geom_hline(yintercept = 0.5, colour = "grey50", lty = 1, size = 0.5) +
  coord_flip() +
  labs(y = "substitution frequency",
       fill = "") +
  scale_y_continuous(expand = c(0.005, 0),
                     labels = scales::percent_format()) +
  scale_fill_manual(values = c("#80DEEA", "#FFAB91")) +
  theme_classic(base_size = 8) +
  theme(legend.position = "top",
        legend.key.size = unit(0.6, "line"),
        legend.margin=margin(0, 0, 0, 0),
        legend.box.margin=margin(-3, -3, -3, -3),
        axis.text.y = element_text(size = 6),
        axis.line.y = element_blank(),
        axis.title.y = element_blank()) +
  guides(fill = guide_legend(nrow = 1))

# Panel c: histogram of the MAF
panel_c <- mafs %>%
  pivot_longer(-sample, names_to = "bin", values_to = "n") %>%
  mutate(decile_bin = floor((as.numeric(bin) + 1) / 2)) %>%
  group_by(decile_bin) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  filter(decile_bin > 3) %>%
  mutate(frac = n / sum(n),
         label = c("]30%, 40%]",
                   "]40%, 50%]",
                   "]50%, 60%]",
                   "]60%, 70%]",
                   "]70%, 80%]",
                   "]80%, 90%]",
                   "]90%, 100%]")) %>%
  ggplot(aes(x = label, y = frac)) +
  geom_col() +
  labs(x = "alternative allele frequency",
       y = "frequency") +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(size = 7))


# Stitch together
plt <- panel_a + panel_b + panel_c +
  plot_layout(heights = c(0.75, 0.25), design = "AB
                                                 CC") +
  plot_annotation(tag_levels = "a")

# Save
ggsave(snakemake@output[["pdf"]],
       plot = plt, width = 160, height = 130, units = "mm", useDingbats = F)
ggsave(snakemake@output[["png"]],
       plot = plt, width = 160, height = 130, units = "mm", dpi = 300)
