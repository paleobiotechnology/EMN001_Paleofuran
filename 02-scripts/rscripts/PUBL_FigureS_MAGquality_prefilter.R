library(tidyverse)
library(patchwork)

theme_set(theme_classic(base_size = 8))

# Load data
mags_preflt <- readxl::read_excel(snakemake@params[["dataset_s3"]],
                                  sheet = 2, skip = 2) %>%
  mutate(sample = str_sub(binID, 1, 6),
         sampletype = if_else(str_detect(binID, "(JAE|VLC)"), "modern", "ancient"),
         MIMAG = if_else(`checkM completeness [%]` >= 90 & `checkM contamination [%]` < 5,
                         "HQ", "MQ"))

# Panel a: overview of the number of HQ and MQ MAGs per sample

panel_a <- mags_preflt %>%
  group_by(sample, sampletype, MIMAG) %>%
  count() %>%
  ggplot(aes(x = sampletype, y = n, group = sampletype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = sampletype), size = 1.75, pch = 21, width = .25) +
  facet_wrap(~ MIMAG, nrow = 1) +
  labs(x = "MAG quality following the MIMAG",
       y = "number of MAGs",
       fill = "sample type") +
  theme(legend.position = "top",
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        panel.grid.major.y = element_line(size = 0.5, colour = "grey50", linetype = 3),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 8),
        strip.background = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 2.5, alpha = 1)))

# Panel b: BUSCO completeness estimates
busco_stats <- mags_preflt %>%
  select(binID, sampletype, MIMAG, contains("BUSCO generic")) %>%
  pivot_longer(contains("BUSCO generic"), names_to = "metric", values_to = "frac") %>%
  filter(!is.na(frac)) %>%
  mutate(metric = str_match(metric, "BUSCO (specific|generic) ([a-z]+) .+")[,3],
         metric = factor(metric, levels = c("complete", "fragmented", "missing")),
         frac = frac / 100) 

panel_b <- busco_stats %>%
  ggplot(aes(x = sampletype, y = frac, group = sampletype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = sampletype), size = 1.75, pch = 21, width = .25, alpha = .3) +
  facet_grid(MIMAG ~ metric) +
  labs(x = "BUSCO's completeness estimate",
       y = "fraction of the marker genes",
       fill = "sample type") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(legend.position = "top",
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        panel.grid.major.y = element_line(size = 0.5, colour = "grey50", linetype = 3),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(size = 8),
        strip.background = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 2.5, alpha = 1)))

# GUNC related analyses
gunc_stats <- mags_preflt %>%
  select(binID, sampletype, MIMAG,
         contamination = `GUNC contamination [%]`,
         `surplus clades` = `GUNC effective no. of surplus clades`,
         CSS = `GUNC clade separation score`) %>%
  mutate(contamination = contamination / 100) %>%
  pivot_longer(contamination:CSS, names_to = "metric", values_to = "value")

# Panel c: GUNC surplus clades
panel_c <- gunc_stats %>%
  filter(metric == "surplus clades") %>%
  ggplot(aes(x = sampletype, y = value, group = sampletype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = sampletype), size = 1.75, pch = 21, width = .25, alpha = .7) +
  facet_wrap(~ MIMAG, nrow = 1) +
  labs(x = "MAG quality following the MIMAG",
       y = "GUNC - eff. no. of surplus clades",
       fill = "sample type") +
  theme(legend.position = "top",
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        panel.grid.major.y = element_line(size = 0.5, colour = "grey50", linetype = 3),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 7),
        strip.text = element_text(size = 8),
        strip.background = element_blank()) +
guides(fill = guide_legend(override.aes = list(size = 2.5, alpha = 1)))

# Panel d: GUNC contamination
panel_d <- gunc_stats %>%
  filter(metric == "contamination") %>%
  ggplot(aes(x = sampletype, y = value, group = sampletype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0.05, size = 0.5, colour = "red", lty = 2) +
  geom_jitter(aes(fill = sampletype), size = 1.75, pch = 21, width = .25, alpha = .7) +
  facet_wrap(~ MIMAG, nrow = 1) +
  labs(x = "MAG quality following the MIMAG",
       y = "GUNC - contamination",
       fill = "sample type") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(legend.position = "top",
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        panel.grid.major.y = element_line(size = 0.5, colour = "grey50", linetype = 3),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(size = 8),
        strip.background = element_blank()) +
guides(fill = guide_legend(override.aes = list(size = 2.5, alpha = 1)))

# Panel e: GUNC CSS
panel_e <- gunc_stats %>%
  filter(metric == "CSS") %>%
  ggplot(aes(x = sampletype, y = value, group = sampletype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0.45, size = 0.5, colour = "red", lty = 2) +
  geom_jitter(aes(fill = sampletype), size = 1.75, pch = 21, width = .25, alpha = .7) +
  facet_wrap(~ MIMAG, nrow = 1) +
  labs(x = "MAG quality following the MIMAG",
       y = "GUNC - CSS",
       fill = "sample type") +
  theme(legend.position = "top",
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        panel.grid.major.y = element_line(size = 0.5, colour = "grey50", linetype = 3),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(size = 8),
        strip.background = element_blank()) +
guides(fill = guide_legend(override.aes = list(size = 2.5, alpha = 1)))

# Stitch together
plt <- guide_area() + panel_a + panel_b + panel_c + panel_d + panel_e +
  plot_layout(heights = c(0.05, 0.316, 0.316, 0.316), guides = "collect",
              design = "AA
                        BC
                        DC
                        EF") +
  plot_annotation(tag_levels = "a")

# Save
ggsave(snakemake@output[["pdf"]],
       plot = plt, width = 160, height = 160, units = "mm", useDingbats = F)
ggsave(snakemake@output[["png"]],
       plot = plt, width = 160, height = 160, units = "mm", dpi = 300)
