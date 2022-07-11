library(tidyverse)
library(patchwork)

theme_set(theme_classic(base_size = 8))

# Load data
mags_preflt <- readxl::read_excel(snakemake@params[["dataset_s3"]],
                                  sheet = 2, skip = 2) %>%
  mutate(sampletype = if_else(str_detect(binID, "(JAE|VLC)"), "modern", "ancient"),
         MIMAG = if_else(`checkM completeness [%]` >= 90 & `checkM contamination [%]` < 5,
                         "HQ", "MQ"))
mags_postflt <- readxl::read_excel(snakemake@params[["dataset_s3"]],
                                   sheet = 3, skip = 2) %>%
  mutate(sampletype = if_else(str_detect(binID, "(JAE|VLC)"), "modern", "ancient"))

# Panel a: number of contigs
figure5_panela <- left_join(mags_preflt %>%
                            select(binID, sampletype, pre = `no. of contigs`),
                            mags_postflt %>%
                            select(binID, sampletype, post = `no. of contigs`),
                            by = c("binID", "sampletype")) %>%
  mutate(diff = (post - pre) / pre) %>%
  ggplot(aes(x = sampletype, y = diff, group = sampletype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = sampletype), size = 1.75, pch = 21, width = .25, alpha = .7) +
  geom_hline(yintercept = 0, size = 0.5, lty = 1, colour = "black") +
  labs(x = "",
       y = "fraction of contigs removed",
       fill = "sample type") +
  scale_y_continuous(labels = scales::percent_format()) +
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

# Panel b: checkM contamination
figure5_panelb <- left_join(mags_preflt %>%
                            select(binID, sampletype, pre = `checkM contamination [%]`),
                            mags_postflt %>%
                            select(binID, sampletype, post = `checkM contamination [%]`),
                            by = c("binID", "sampletype")) %>%
  mutate(diff = (post - pre) / 100) %>%
  ggplot(aes(x = sampletype, y = diff, group = sampletype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = sampletype), size = 1.75, pch = 21, width = .25, alpha = .7) +
  geom_hline(yintercept = 0, size = 0.5, lty = 1, colour = "black") +
  labs(x = "",
       y = expression(paste(Delta, " checkM contamination")),
       fill = "sample type") +
  scale_y_continuous(labels = scales::percent_format()) +
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

# Panel c: GUNC contamination
figure5_panelc <- left_join(mags_preflt %>%
                            select(binID, sampletype, pre = `GUNC contamination [%]`),
                            mags_postflt %>%
                            select(binID, sampletype, post = `GUNC contamination [%]`),
                            by = c("binID", "sampletype")) %>%
  mutate(diff = (post - pre) / 100) %>%
  ggplot(aes(x = sampletype, y = diff, group = sampletype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = sampletype), size = 1.75, pch = 21, width = .25, alpha = .7) +
  geom_hline(yintercept = 0, size = 0.5, lty = 1, colour = "black") +
  labs(x = "",
       y = expression(paste(Delta, " GUNC contamination")),
       fill = "sample type") +
  scale_y_continuous(labels = scales::percent_format()) +
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

# Panel d: GUNC CSS
figure5_paneld <- left_join(mags_preflt %>%
                            select(binID, sampletype, pre = `GUNC clade separation score`),
                            mags_postflt %>%
                            select(binID, sampletype, post = `GUNC clade separation score`),
                            by = c("binID", "sampletype")) %>%
  mutate(diff = post - pre) %>%
  ggplot(aes(x = sampletype, y = diff, group = sampletype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = sampletype), size = 1.75, pch = 21, width = .25, alpha = .7) +
  geom_hline(yintercept = 0, size = 0.5, lty = 1, colour = "black") +
  labs(x = "",
       y = expression(paste(Delta, " GUNC CSS")),
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

# Panel e: genome size
figure5_panele <- left_join(mags_preflt %>%
                            select(binID, sampletype, pre = `genome size [Mb]`),
                            mags_postflt %>%
                            select(binID, sampletype, post = `genome size [Mb]`),
                            by = c("binID", "sampletype")) %>%
  mutate(diff = post - pre) %>%
  ggplot(aes(x = sampletype, y = diff, group = sampletype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = sampletype), size = 1.75, pch = 21, width = .25, alpha = .7) +
  geom_hline(yintercept = 0, size = 0.5, lty = 1, colour = "black") +
  labs(x = "",
       y = expression(paste(Delta, " genome size [Mb]")),
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

# Panel f: checkM completeness
figure5_panelf <- left_join(mags_preflt %>%
                            select(binID, sampletype, pre = `checkM completeness [%]`),
                            mags_postflt %>%
                            select(binID, sampletype, post = `checkM completeness [%]`),
                            by = c("binID", "sampletype")) %>%
  mutate(diff = (post - pre) / 100) %>%
  ggplot(aes(x = sampletype, y = diff, group = sampletype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(fill = sampletype), size = 1.75, pch = 21, width = .25, alpha = .7) +
  geom_hline(yintercept = 0, size = 0.5, lty = 1, colour = "black") +
  labs(x = "",
       y = expression(paste(Delta, " checkM completeness")),
       fill = "sample type") +
  scale_y_continuous(labels = scales::percent_format()) +
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

# Panel g: polymorphism rate
figure5_panelg <- mags_postflt %>%
  filter(MIMAG != "low") %>%
  mutate(MIMAG = recode(MIMAG, `high` = "HQ", `medium` = "MQ")) %>%
  select(binID, sampletype, MIMAG, ratio = `ratio non-syn. to syn. minor alleles [%]`) %>%
  mutate(ratio = ratio / 100) %>%
  ggplot(aes(x = sampletype, y = ratio, group = sampletype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0.01, size = 0.5, colour = "red", lty = 2) +
  geom_jitter(aes(fill = sampletype), size = 2, pch = 21, width = .25, alpha = .7) +
  facet_wrap(~ MIMAG, nrow = 1) +
  labs(x = "",
       y = "ratio non-syn. to syn. minor alleles",
       fill = "sample type") +
  scale_y_continuous(labels = scales::percent_format(), breaks = c(0, 0.01, 0.02)) +
  theme(legend.position = "top",
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        panel.grid.major.y = element_line(size = 0.5, colour = "grey50", linetype = 3),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 7),
        strip.text = element_text(size = 8),
        strip.background = element_blank()) +
guides(fill = guide_legend(override.aes = list(size = 2.5, alpha = 1)))

# Panel h: number of MAGs per MIMAG group
figure5_panelh <- left_join(mags_preflt %>%
                            select(binID, sampletype, pre = MIMAG) %>%
                            mutate(pre = recode(pre, `HQ` = "high", `MQ` = "medium")),
                            mags_postflt %>%
                            select(binID, sampletype, post = MIMAG),
                            by = c("binID", "sampletype")) %>%
  pivot_longer(pre:post, names_to = "filterstep", values_to = "MIMAG") %>%
  mutate(filterstep = factor(filterstep, levels = c("pre", "post")),
         MIMAG = factor(MIMAG, levels = c("high", "medium", "low"))) %>%
  group_by(sampletype, filterstep, MIMAG) %>%
  count() %>%
  pivot_wider(names_from = "filterstep", values_from = "n", values_fill = 0) %>%
  mutate(diff = post - pre) %>%
  ggplot(aes(x = MIMAG, y = diff)) +
  geom_col(aes(fill = MIMAG)) +
  facet_wrap(~ sampletype, nrow = 1) +
  scale_fill_brewer(palette = "Reds", direction = -1) +
  labs(x = "",
       y = expression(paste(Delta, "number of MAGs")),
       fill = "MIMAG category") +
  theme(legend.position = "top",
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.key.size = unit(4, "mm"),
        panel.grid.major.y = element_line(size = 0.5, colour = "grey50", linetype = 3),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 7),
        strip.text = element_text(size = 8),
        strip.background = element_blank()) +
  guides(fill = guide_legend(nrow = 1))

# Stitch together
plt <- guide_area() + figure5_panela + figure5_panelb + figure5_panelc +
       figure5_paneld + figure5_panele + figure5_panelf +
       figure5_panelh + figure5_panelg +
  plot_layout(heights = c(0.1, 0.45, 0.45), guides = "collect",
              design = c("AAAA
                          BCDE
                          FGHI")) +
  plot_annotation(tag_levels = "a") & theme(legend.position = "top")

# Save
ggsave(snakemake@output[["pdf"]],
       plot = plt, width = 160, height = 160, units = "mm", useDingbats = F)
ggsave(snakemake@output[["png"]],
       plot = plt, width = 160, height = 160, units = "mm", dpi = 300)

