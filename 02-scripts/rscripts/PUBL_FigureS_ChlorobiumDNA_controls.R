library(data.table)
library(tidyverse)
library(patchwork)
library(ggh4x)

theme_set(theme_classic(base_size = 9))

# Load data
dentalcalculus <- fread(snakemake@params[["dentalcalculus"]])
controls <- readxl::read_excel(snakemake@params[["dataset_s7"]],
                               sheet = 3, skip = 2)

# Calculate relative abundance
relab <- bind_rows(dentalcalculus %>%
  mutate(relab = alignedReads / totalReads,
         sampletype = "dental calculus") %>%
  select(sample, relab, sampletype),
  controls %>%
  select(sample, relab = `% aligned`) %>%
  mutate(sampletype = case_when(
      sample == "ElMiron" ~ "toe bone",
      str_sub(sample, 1, 3) == "EMN" ~ "sediments",
      T ~ "lab controls"),
    relab = relab / 100)) %>%
  mutate(sampletype = factor(sampletype,
                             levels = c("dental calculus",
                                        "sediments",
                                        "toe bone",
                                        "lab controls")))

# Plot
plt <- ggplot(relab,
              aes(x = sampletype, y = relab, group = sampletype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = sampletype), size = 2.25, pch = 21, width = .25, alpha = .7) +
  labs(x = "",
       y = "relative abundance",
       fill = "sample type") +
  scale_y_log10(breaks = c(0.0001, 0.01, 0.1),
                labels = c("0.1%", "1%", "10%"),
                guide = "axis_logticks") +
  theme(legend.position = "none",
        panel.grid.major.y = element_line(size = 0.5, colour = "grey50", linetype = 3),
        axis.line = element_blank(),
        axis.ticks.length.y = unit(3, "mm"),
        ggh4x.axis.ticks.length.minor = rel(0.5),
        ggh4x.axis.ticks.length.mini = rel(0.3),
        strip.text = element_text(size = 8),
        strip.background = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size = 2.5, alpha = 1)))

# Save
ggsave(snakemake@output[["pdf"]],
       plot = plt, width = 160, height = 100, units = "mm", useDingbats = F)
ggsave(snakemake@output[["png"]],
       plot = plt, width = 160, height = 100, units = "mm", dpi = 300)
