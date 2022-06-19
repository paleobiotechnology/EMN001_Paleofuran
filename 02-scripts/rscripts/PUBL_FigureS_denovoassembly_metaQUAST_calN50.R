library(data.table)
library(tidyverse)
library(patchwork)
library(ggh4x)

theme_set(theme_classic(base_size = 9))

# Load data
assembly <- readxl::read_excel(snakemake@params[["assembly"]],
                                sheet = 1, skip = 2) %>%
  mutate(sampletype = if_else(str_sub(sample, 1, 3) %in% c("JAE", "VLC"),
                              "modern", "ancient"))

# Plot
## Panel a: fraction of contigs with a certain minimal length
panel_a_plt <- assembly %>%
  mutate(across(contains("# of contigs >="), ~ . / `# of contigs`)) %>%
  select(sample, contains("# of contigs >=")) %>%
  pivot_longer(-sample, names_to = "minL", values_to = "nContigs") %>%
  mutate(minL = str_replace(str_replace(minL, "# of contigs ", ""), ",000 bp", "kb"),
         sampletype = if_else(str_sub(sample, 1, 3) %in% c("JAE", "VLC"),
                              "modern", "ancient")) %>%
  mutate(minL = factor(minL, levels = unique(.$minL)),
         nContigs = nContigs + 0.00001) %>%
  ggplot(aes(x = sampletype, y = nContigs)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = sampletype), size = 2, pch = 21, width = .25) +
  facet_wrap(~ minL, nrow = 1) +
  labs(x = "minimal contig length ",
       y = "fraction of the contigs",
       fill = "sample type") +
  scale_y_log10(labels = c("0%", "0.1%", "1%", "10%", "50%", "100%"),
                breaks = c(0.00001, 0.001, 0.01, 0.1, 0.5, 1)) +
  theme(legend.position = "top",
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        panel.grid.major.y = element_line(size = 0.5, colour = "grey50", linetype = 3),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 9)) +
  guides(fill = guide_legend(override.aes = list(size = 3)))

## Panel b: Nx plot
panel_b_plt <- assembly %>%
  select(sample, starts_with("N")) %>%
  pivot_longer(-sample, names_to = "decile", values_to = "length") %>%
  mutate(sampletype = if_else(str_sub(sample, 1, 3) %in% c("JAE", "VLC"),
                              "modern", "ancient"),
         decile = factor(decile, levels = unique(.$decile)),
         length = length / 1e3) %>%
  filter(decile != "N100") %>%
  ggplot(aes(x = sampletype, y = length, group = sampletype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = sampletype), size = 2, pch = 21, width = .25) +
  facet_wrap(~ decile, nrow = 1) +
  labs(x = "contig deciles",
       y = "contig length [kb]",
       fill = "sample type") +
  scale_y_log10(limit = c(0.5, NA),
                breaks = c(1, 5, 10, 50, 100, 500, 1000),
                guide = "axis_logticks") +
  theme(legend.position = "top",
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        panel.grid.major.y = element_line(colour = "grey50", size=0.5, linetype = 3),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.length.y = unit(3, "mm"),
        ggh4x.axis.ticks.length.minor = rel(0.5),
        ggh4x.axis.ticks.length.mini = rel(0.3),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 9)) +
  guides(fill = guide_legend(override.aes = list(size = 3)))

## Stitch together
plt <- guide_area() + panel_a_plt + panel_b_plt +
  plot_layout(ncol = 1, guides = "collect", heights = c(0.05, 0.475, 0.475)) +
  plot_annotation(tag_levels = "a")

# Save
ggsave(snakemake@output[["pdf"]],
       plot = plt, width = 160, height = 160, units = "mm", useDingbats = F)
ggsave(snakemake@output[["png"]],
       plot = plt, width = 160, height = 160, units = "mm", dpi = 300)
