library(Matrix)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(dplyr)
library(stringr)
library(reshape2)
library(patchwork)

source(file.path("functions", "FileIO_HelperFunctions.R"))
source(file.path("functions", "Step14_Analysis_HelperFunctions.R"))
source(file.path("functions", "Step15_Plotting_HelperFunctions.R"))

options(scipen = 999)

granularity <- c("broad_class")

bulk_datasets <- c("Mayo", "MSBB", "ROSMAP")

best_errors_list_broad <- Load_BestErrors("broad_class")
best_errors_list_sub <- Load_BestErrors("sub_class")

# Unpack variables into environment for readability
list2env(best_errors_list_broad, globalenv())
list2env(best_errors_list_sub, globalenv())

tissues_use <- c("TCX", "FP", "ACC") # When subsetting
tissues_use_full <- c("Mayo TCX", "MSBB FP", "ROSMAP ACC") # When subsetting


# Color setup ------------------------------------------------------------------

algs <- unique(best_errors_broad$algorithm)
regs <- unique(best_errors_broad$regression_method)

# Algorithms
algorithm_colors <- RColorBrewer::brewer.pal(length(algs), "Set2")
names(algorithm_colors) <- sort(algs)

# Regression methods
# F45533CC is several shades lighter than viridis::turbo(1, alpha = 0.8, begin = 0.9, end = 0.9)
regression_colors <- c(viridis::turbo(length(regs)-1, alpha = 0.8, begin = 0.15, end = 0.5),
                       "#F45533CC")

# Tissues (modified viridis turbo color scheme)
tiss <- colSums(table(best_errors_broad$tissue, best_errors_broad$test_data_name) > 1)
tissue_colors <- c(viridis::turbo(tiss[["Mayo"]], begin = 0.1, end = 0.2),
                   viridis::turbo(tiss[["MSBB"]], begin = 0.35, end = 0.6),
                   viridis::turbo(tiss[["ROSMAP"]], begin = 0.7, end = 0.9))
names(tissue_colors) <- sort(unique(best_errors_broad$tissue_full))

# MSBB colors need to be darker and a little more differentiated -- original
# colors are 20EAABFF, 67FD68FF, AEFA37FF, E1DD37FF
tissue_colors[grepl("MSBB", names(tissue_colors))] <- c("#00CA8BFF", "#37CD38FF",
                                                        "#7ECA07FF", "#C1BD17FF")

# Slightly more pastel than default looks better
tissue_fill_colors <- str_replace(tissue_colors, "FF$", "88")
names(tissue_fill_colors) <- names(tissue_colors)


# Figure 1 variations ----------------------------------------------------------

pdf(file.path(dir_figures,
              str_glue("draft_figure1_variations.pdf")),
    width=7.5, height = 10)


# Figure 1B --------------------------------------------------------------------

baselines_plot_tissue <- baselines_broad %>%
  subset(reference_data_name != "All zeros") %>%
  Create_BoxStats(c("tissue_full", "error_type")) %>%
  group_by(error_type) %>%
  mutate(best_val = if (unique(error_type) == "Correlation") max_val else min_val) %>%
  ungroup()

errs_box_noalg <- best_errors_broad %>%
  Create_BoxStats(c("tissue_full", "error_type",
                    "normalization", "regression_method"))

cap <- "Show correlation only, for all tissues. RMSE and MAPE graphs would be in supplementary."
plt1Bv1 <- ggplot(subset(errs_box_noalg, error_type == "Correlation"),
                  aes(x = normalization, fill = regression_method,
                      y = median_val, ymin = min_val, ymax = max_val)) +
  geom_hline(aes(yintercept = best_val, color = "Baseline"),
             data = subset(baselines_plot_tissue, error_type == "Correlation"),
             linetype = "twodash") +
  geom_crossbar(position = position_dodge2(padding = 0),
                fatten = 1.5, width = 0.75) +
  theme_bw() +
  scale_fill_manual(values = regression_colors) +
  scale_color_manual(values = "darkgray") +
  facet_wrap(~tissue_full, nrow = 3) +
  labs(fill = "Regression Method", color = "", caption = cap) +
  xlab("Normalization Strategy") +
  ylab("Correlation") +
  ggtitle("Figure 1B v1: Correlation: normalization vs regression")

cap <- paste("Original, but cut down: show 3 representative tissues for correlation\n",
             "only and kick the rest to supplementary.The sub class graph looks\n",
             "nearly identical so I don't think it's necessary to show both.")
plt1Bv2 <- ggplot(subset(errs_box_noalg, error_type == "Correlation" &
                           tissue_full %in% tissues_use_full),
                  aes(x = normalization, fill = regression_method,
                      y = median_val, ymin = min_val, ymax = max_val)) +
  geom_hline(aes(yintercept = best_val, color = "Baseline"),
             data = subset(baselines_plot_tissue, error_type == "Correlation" &
                             tissue_full %in% tissues_use_full),
             linetype = "twodash") +
  geom_crossbar(position = position_dodge2(padding = 0),
                fatten = 1.5, width = 0.75) +
  theme_bw() +
  scale_fill_manual(values = regression_colors) +
  scale_color_manual(values = "slategray") +
  facet_wrap(~tissue_full, nrow = 1) +
  labs(fill = "Regression Method", color = "", caption = cap) +
  xlab("Normalization Strategy") +
  ylab("Correlation") +
  ggtitle("Figure 1B v2: Correlation: normalization vs regression")
print(plt1Bv1 / plt1Bv2 + plot_layout(heights = c(3, 1)))

cap <- paste("Showing all 3 error metrics for select tissues. I'm not sure if it\n",
             "will be confusing that for correlation, higher values are better but\n",
             "for the other two, lower values are better, and the baseline value\n",
             "is set to the highest/lowest baseline accordingly. We could also\n",
             "put more separation between the three plots so they look like 3\n",
             "separate plots instead of the current facet_grid. TPM RMSE for Mayo\n",
             "TCX is really high, which messes with the scale of those graphs.\n",
             "Not sure what a good solution to that is since if we're using\n",
             "representative tissues, the only other option for Mayo is CBE.")

plt1Bv3 <- ggplot(subset(errs_box_noalg,
                            tissue_full %in% tissues_use_full),
                  aes(x = normalization, fill = regression_method,
                      y = median_val, ymin = min_val, ymax = max_val)) +
  geom_hline(aes(yintercept = best_val, color = "Baseline"),
             data = subset(baselines_plot_tissue,
                             tissue_full %in% tissues_use_full),
             linetype = "twodash") +
  geom_crossbar(position = position_dodge2(padding = 0),
                fatten = 1.5, width = 0.75) +
  theme_bw() +
  scale_fill_manual(values = regression_colors) +
  scale_color_manual(values = "slategray") +
  facet_grid(error_type ~ tissue_full, scales = "free", switch = "y") +
  labs(fill = "Regression Method", color = "", caption = cap) +
  xlab("Normalization Strategy") +
  ylab(NULL) +
  ggtitle("Figure 1B v3: Error metrics: normalization vs regression") +
  theme(strip.background.y = element_blank(),
        strip.placement.y = "outside")
print(plt1Bv3 / plot_spacer() + plot_layout(heights = c(3, 1)))


# Figure 1C --------------------------------------------------------------------

cap <- paste("Took the top 3 scoring estimates for each tissue and error metric\n",
             "(9 total per tissue) and counted how many times each norm/regression\n",
             "appeared. Legend would need to be hand-fixed for this graph. These can\n",
             "also be changed to 3x3 grids. These graphs look very similar except a\n",
             "few tissues. Do we want to show both in the main figure or put sub\n",
             "class in supplemental?")

plt1Ca <- ggplot(Count_Grouped(ranked_df_broad, c("tissue", "normalization", "regression_method")),
                 aes(x = normalization, y = regression_method, color = count, size = count)) +
  geom_count() + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.width = unit(10, "pt")) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5) +
  facet_wrap(~tissue, nrow = 1) +
  labs(color = "Count") +
  xlab(NULL) +
  ylab(NULL) +
  guides(size = "none") +
  ggtitle("Figure 1C: Best normalization/regression per tissue", subtitle = "Broad class")

plt1Cb <- ggplot(Count_Grouped(ranked_df_sub, c("tissue", "normalization", "regression_method")),
                 aes(x = normalization, y = regression_method, color = count, size = count)) +
  geom_count() + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.width = unit(10, "pt")) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5) +
  facet_wrap(~tissue, nrow = 1) +
  labs(color = "Count", caption = cap) +
  xlab(NULL) +
  ylab(NULL) +
  guides(size = "none") +
  ggtitle(NULL, subtitle = "Sub class")
print(plt1Ca / plt1Cb / plot_spacer() + plot_layout(heights = c(1, 1, 2)))

# Figure 1D --------------------------------------------------------------------

cap <- paste("Took the top 3 scoring estimates for each tissue and error metric\n",
             "(9 total per tissue) and counted how many times each reference\n",
             "appeared. This can also be changed to a 3x3 grid. Broad class and\n",
             "sub class look nearly identical, should sub class go in supplemental?")
plt1Dv1a <- ggplot(Count_Grouped(ranked_df_broad, c("tissue", "reference_data_name")),
                   aes(x = tissue, y = reference_data_name, color = count, size = count)) +
  geom_count() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.width = unit(10, "pt")) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5) +
  labs(color = "Count") +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Figure 1D: Best reference", subtitle = "Broad class") +
  guides(size = "none")

plt1Dv1b <- ggplot(Count_Grouped(ranked_df_sub, c("tissue", "reference_data_name")),
                   aes(x = tissue, y = reference_data_name, color = count, size = count)) +
  geom_count() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.width = unit(10, "pt")) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5) +
  labs(color = "Count", caption = cap) +
  xlab("Tissue") +
  ylab(NULL) +
  ggtitle(NULL, subtitle = "Sub class") +
  guides(size = "none")

print(plt1Dv1a + plt1Dv1b)

dev.off()


# Draft Main Figure 1 ----------------------------------------------------------

row1 <- ggplot() + ggtitle("A: Placeholder for diagram of process")
row2 <- plt1Bv3 + ggtitle("B") + labs(caption = NULL) +
  theme(legend.title = element_text(size = 10))
row3 <- plt1Ca + ggtitle("C", subtitle = NULL)
row4 <- (plt1Dv1a + ggtitle("D", subtitle = NULL) + labs(caption = NULL)) +
  (ggplot() + ggtitle("E: Placeholder for ??"))

#row3 <- plt1Ca + ggtitle("C") + labs(caption = NULL) +
#  theme(legend.title = element_text(size = 10))
#row4 <- (plt1Dv3 + ggtitle("D", subtitle = "Unfiltered") + labs(caption = NULL) +
#           ylim(0, 0.85)) +
#  (plt2v4 + ggtitle(NULL, subtitle = "Filtered") + ylim(0, 0.85)) +
#  (plt4v2 + ggtitle("E") + labs(caption = NULL) +
#     theme(legend.title = element_text(size = 10))) +
#  plot_layout(widths = c(1, 1, 2))
#row5 <- plt5v2 + ggtitle("F") + labs(caption = NULL) +
#  theme(legend.title = element_text(size = 10))

pdf(file.path(dir_figures,
              str_glue("draft_figure1.pdf")),
    width=7.5, height = 10)

print((row1) / (row2) / (row3) / (row4) + plot_layout(heights = c(2, 2, 1, 1)) +
        plot_annotation("Draft main figure 1 -- please ignore weird formatting / sizing for now"))

dev.off()


# Figure 2 variations ----------------------------------------------------------

pdf(file.path(dir_figures,
              str_glue("draft_figure2_variations.pdf")),
    width=7.5, height = 10)


# Figure 2A --------------------------------------------------------------------

cap <- paste("Took the top 3 scoring estimates for each tissue and error metric\n",
             "(9 total per tissue) and counted how many times each algorithm\n",
             "appeared. This can also be changed to a 3x3 grid.")

plt2Av1a <- ggplot(Count_Grouped(ranked_df_broad, c("tissue", "algorithm", "error_type")),
                   aes(x = error_type, y = algorithm, color = count, size = count)) +
  geom_count() + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.width = unit(10, "pt")) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5) +
  facet_wrap(~tissue, nrow = 1) +
  labs(color = "Count") +
  xlab(NULL) +
  ylab(NULL) +
  guides(size = "none") +
  ggtitle("Figure 2A v1: Best algorithm per tissue for each error metric",
          subtitle = "Broad class")

# Sub class is missing several algorithms that dropped out, fill these in
alg_info <- Count_Grouped(ranked_df_broad, c("tissue", "algorithm", "error_type")) %>%
  select(tissue, algorithm)
alg_sub <- Count_Grouped(ranked_df_sub, c("tissue", "algorithm", "error_type")) %>%
  merge(alg_info, by = c("tissue", "algorithm"), all = TRUE)
alg_sub$error_type[is.na(alg_sub$error_type)] <- "Correlation"

plt2Av1b <- ggplot(alg_sub,
                   aes(x = error_type, y = algorithm, color = count, size = count)) +
  geom_count() + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.width = unit(10, "pt")) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5) +
  facet_wrap(~tissue, nrow = 1) +
  labs(color = "Count", caption = cap) +
  xlab(NULL) +
  ylab(NULL) +
  guides(size = "none") +
  ggtitle(NULL, subtitle = "Sub class")

cap <- paste("Original: get rid of error metric to make the graph simpler.")
plt2Av2a <- ggplot(Count_Grouped(ranked_df_broad, c("tissue", "algorithm")),
                   aes(x = tissue, y = algorithm, color = count, size = count)) +
  geom_count() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.width = unit(10, "pt")) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5) +
  xlab(NULL) +
  ylab(NULL) +
  guides(color = "none", size = "none") +
  ggtitle("Figure 2A v2: Best algorithm per tissue", subtitle = "Broad class")

# Sub class is missing several algorithms that dropped out, fill these in
alg_info <- Count_Grouped(ranked_df_broad, c("tissue", "algorithm")) %>%
  select(tissue, algorithm)
alg_sub <- Count_Grouped(ranked_df_sub, c("tissue", "algorithm")) %>%
  merge(alg_info, by = c("tissue", "algorithm"), all = TRUE)

plt2Av2b <- ggplot(alg_sub,
                   aes(x = tissue, y = algorithm, color = count, size = count)) +
  geom_count() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.width = unit(10, "pt")) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5) +
  labs(color = "Count", caption = cap) +
  xlab(NULL) +
  ylab(NULL) +
  guides(size = "none") +
  ggtitle(NULL, subtitle = "Sub class")

plt2Av2 <- plt2Av2a + plt2Av2b

print(plt2Av1a / plt2Av1b / (plt2Av2 + plot_layout(widths = c(1,1)))) # + plot_layout(heights = c(2, 2, 1)))


# Figure 2B --------------------------------------------------------------------

baselines_plot <- baselines_broad %>%
  subset(reference_data_name != "All zeros") %>%
  Create_BoxStats(grouping_cols = c("tissue_full", "normalization",
                                    "regression_method", "error_type"))

errs_better <- merge(best_errors_broad, baselines_plot,
                     by = c("tissue_full", "normalization", "regression_method",
                            "error_type")) %>%
  group_by(error_type) %>%
  mutate(better = if(unique(error_type) == "Correlation") (value > max_val) else (value < min_val))

better_stats <- errs_better %>%
  group_by(tissue_full, algorithm) %>%
  dplyr::summarize(count = n(),
                   pct_better_than_baseline = sum(better) / n(),
                   .groups = "drop")

baselines_plot_filt <- baselines_top_broad %>%
  subset(reference_data_name != "All zeros") %>%
  mutate(data_transform = paste(normalization, "+", regression_method)) %>%
  merge(best_dt_broad, by = c("tissue", "algorithm", "data_transform"), all = FALSE) %>%
  Create_BoxStats(grouping_cols = c("tissue_full", "normalization",
                                    "regression_method", "error_type"))

errs_better_filt <- best_errors_top_broad %>%
  mutate(data_transform = paste(normalization, "+", regression_method)) %>%
  merge(best_dt_broad, by = c("tissue", "algorithm", "data_transform"), all = FALSE) %>%
  merge(baselines_plot_filt,
        by = c("tissue_full", "error_type")) %>%
  group_by(error_type) %>%
  mutate(better = if(unique(error_type) == "Correlation") (value > max_val) else (value < min_val))

better_stats_filt <- errs_better_filt %>%
  group_by(tissue_full, algorithm) %>%
  dplyr::summarize(count = n(),
                   pct_better_than_baseline = sum(better) / n(),
                   .groups = "drop")

baselines_plot_filt_sub <- baselines_top_sub %>%
  subset(reference_data_name != "All zeros") %>%
  mutate(data_transform = paste(normalization, "+", regression_method)) %>%
  merge(best_dt_broad, by = c("tissue", "algorithm", "data_transform"), all = FALSE) %>%
  Create_BoxStats(grouping_cols = c("tissue_full", "normalization",
                                    "regression_method", "error_type"))

errs_better_filt_sub <- best_errors_top_sub %>%
  mutate(data_transform = paste(normalization, "+", regression_method)) %>%
  merge(best_dt_sub, by = c("tissue", "algorithm", "data_transform"), all = FALSE) %>%
  merge(baselines_plot_filt_sub,
        by = c("tissue_full", "error_type")) %>%
  group_by(error_type) %>%
  mutate(better = if(unique(error_type) == "Correlation") (value > max_val) else (value < min_val))

# Sub class is missing some algorithms, this fills them in
alg_info <- errs_better_filt %>%
  ungroup() %>%
  select(tissue, algorithm) %>%
  distinct()

errs_better_filt_sub <- errs_better_filt_sub %>%
  merge(alg_info, by = c("tissue", "algorithm"), all = TRUE)
errs_better_filt_sub$better[is.na(errs_better_filt_sub$better)] <- FALSE

cap <- paste("This graph is filtered to the top 3 errors for each algorithm, for\n",
             "the best norm/regression for each tissue. Breaking out by\n",
             "error metric doesn't provide better information. This could also\n",
             "be a supplementary figure if no room. This graph is for broad class.")
plt2Bv1 <- ggplot(better_stats_filt,
                  aes(x = algorithm, y = pct_better_than_baseline, fill = algorithm)) +
  geom_col() +
  facet_wrap(~tissue_full) +
  scale_fill_manual(values = algorithm_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(caption = cap) +
  xlab(NULL) +
  ylab("Percent") +
  guides(fill = "none") +
  ggtitle("Figure 2B v1: Percent of test results better than baseline", subtitle = "Filtered")

cap <- paste("Showing unfiltered (all norms/regressions/references) vs filtered,\n",
             "using 3 representative tissues and putting the rest in supplemental\n",
             "figures. This graph is for broad class. I'm not sure if it's interesting\n",
             "to show unfiltered values or not.")
plt2Bv2a <- ggplot(subset(better_stats, tissue_full %in% tissues_use_full),
                   aes(x = algorithm, y = pct_better_than_baseline, fill = algorithm)) +
  geom_col() +
  facet_wrap(~tissue_full, ncol = 1) +
  scale_fill_manual(values = algorithm_colors) +
  ylim(0,1.05) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab(NULL) +
  ylab("Percent") +
  guides(fill = "none") +
  ggtitle("Figure 2B v2", subtitle = "Unfiltered")

plt2Bv2b <- ggplot(subset(better_stats_filt, tissue_full %in% tissues_use_full),
                   aes(x = algorithm, y = pct_better_than_baseline, fill = algorithm)) +
  geom_col() +
  facet_wrap(~tissue_full, ncol = 1) +
  scale_fill_manual(values = algorithm_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 1.05) +
  labs(caption = cap) +
  xlab(NULL) +
  ylab("Percent") +
  guides(fill = "none") +
  ggtitle(NULL, subtitle = "Filtered")

plt2Bv2 <- plt2Bv2a + plt2Bv2b

better_stats_filt2 <- errs_better_filt %>%
  group_by(algorithm) %>%
  dplyr::summarize(count = n(),
                   pct_better_than_baseline = sum(better) / n(),
                   .groups = "drop")

better_stats_filt_sub2 <- errs_better_filt_sub %>%
  group_by(algorithm) %>%
  dplyr::summarize(count = n(),
                   pct_better_than_baseline = sum(better) / n(),
                   .groups = "drop")
better_stats_filt_sub2$pct_better_than_baseline[better_stats_filt_sub2$pct_better_than_baseline == 0] <- 0.01

cap <- paste("The original that appears in the draft main figure, all data is\n",
             "collapsed across all tissues, and we can put the tissue-specific\n",
             "breakout in supplemental.")
plt2Bv3a <- ggplot(better_stats_filt2,
                   aes(x = algorithm, y = pct_better_than_baseline, fill = algorithm)) +
  geom_col() +
  scale_fill_manual(values = algorithm_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 1.05) +
  xlab(NULL) +
  ylab("Percent") +
  guides(fill = "none") +
  ggtitle("Figure 2B v3", subtitle = "Broad class")

plt2Bv3b <- ggplot(better_stats_filt_sub2,
                   aes(x = algorithm, y = pct_better_than_baseline, fill = algorithm)) +
  geom_col() +
  scale_fill_manual(values = algorithm_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, 1.05) +
  labs(caption = cap) +
  xlab(NULL) +
  ylab("Percent") +
  ggtitle(NULL, subtitle = "Sub class") +
  guides(fill = "none")

plt2Bv3 <- plt2Bv3a + plt2Bv3b

print((plt2Bv1 + plt2Bv2 + plot_layout(widths = c(3, 2))) /
        (plt2Bv3 + plot_spacer() + plot_layout(widths = c(1, 1, 3))) + plot_layout(heights = c(3, 1)))


# Figure 2? --------------------------------------------------------------------

# TODO get the actual mean inhibitory:excitatory ratio, not the percent of samples
inh_mean <- quality_stats_broad %>%
  group_by(tissue, algorithm) %>%
  summarize(pct_bad_inhibitory_ratio = mean(pct_bad_inhibitory_ratio),
            .groups = "drop")

plt <- ggplot(inh_mean, aes(x = algorithm, y = pct_bad_inhibitory_ratio,
                            fill = algorithm)) +
  geom_col() + theme_bw() +
  facet_wrap(~tissue, nrow = 3) +
  scale_fill_manual(values = algorithm_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Percent high inhibitory to excitatory ratio") +
  xlab(NULL) +
  ylab("Average percent of samples") +
  guides(fill = "none")
print(plt)

inh_mean2 <- quality_stats_broad %>%
  group_by(tissue, tissue_full) %>%
  summarize(pct_bad_inhibitory_ratio = mean(pct_bad_inhibitory_ratio),
            .groups = "drop")

plt <- ggplot(inh_mean2, aes(x = tissue, y = pct_bad_inhibitory_ratio,
                             fill = tissue_full)) +
  geom_col() + theme_bw() +
  scale_fill_manual(values = tissue_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Percent high inhibitory to excitatory ratio") +
  xlab(NULL) +
  ylab("Average percent of samples") +
  guides(fill = "none")
print(plt)


# Figure 2? --------------------------------------------------------------------

# TODO fix: this is pulling in pct valid estimates from the best params, not
# from all params

# Fill in stats for reference data sets that didn't have any valid param sets
# and therefore aren't in this data frame
file_params <- expand.grid(tissue = unique(quality_stats_broad$tissue),
                           reference_data_name = unique(quality_stats_broad$reference_data_name),
                           normalization = unique(quality_stats_broad$normalization),
                           regression_method = unique(quality_stats_broad$regression_method),
                           algorithm = unique(quality_stats_broad$algorithm))

# CibersortX does not have tmm normalization. It also has two possible
# reference_input_types, but since one of them (cibersortx) only has 2 possible
# params in the file, we're ignoring that one to avoid skewing the numbers badly.
cibersort_only <- file_params %>%
  subset(algorithm == "CibersortX" & normalization %in% c("CPM", "TPM"))

# Music doesn't have tmm or tpm normalization
music_only <- file_params %>%
  subset(algorithm == "MuSiC" & normalization == "CPM")

file_params <- file_params %>%
  subset(!(algorithm %in% c("CibersortX", "MuSiC"))) %>%
  rbind(cibersort_only, music_only)

qstats_broad <- quality_stats_broad %>%
  select(-param_id, -pct_bad_inhibitory_ratio) %>%
  distinct() %>%
  subset(reference_input_type != "cibersortx") %>%
  merge(file_params, by = colnames(file_params), all = TRUE)

# These values are NA where there was missing data, set to 0
qstats_broad$pct_valid_results[is.na(qstats_broad$pct_valid_results)] <- 0

qstats_mean_broad <- qstats_broad %>%
  group_by(algorithm) %>%
  summarize(pct_valid_results = mean(pct_valid_results),
            .groups = "drop")

qstats_sub <- quality_stats_sub %>%
  select(-param_id, -pct_bad_inhibitory_ratio) %>%
  distinct() %>%
  subset(reference_input_type != "cibersortx") %>%
  merge(file_params, by = colnames(file_params), all = TRUE)

# These values are NA where there was missing data, set to 0
qstats_sub$pct_valid_results[is.na(qstats_sub$pct_valid_results)] <- 0

qstats_mean_sub <- qstats_sub %>%
  group_by(algorithm) %>%
  summarize(pct_valid_results = mean(pct_valid_results),
            .groups = "drop")

plt <- ggplot(qstats_mean_broad,
              aes(x = algorithm, y = pct_valid_results, fill = algorithm)) +
  geom_col() + theme_bw() +
  #facet_wrap(~tissue, nrow = 3) +
  scale_fill_manual(values = algorithm_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Percent valid estimates")
print(plt)

plt <- ggplot(qstats_mean_sub,
              aes(x = algorithm, y = pct_valid_results, fill = algorithm)) +
  geom_col() + theme_bw() +
  #facet_wrap(~tissue, nrow = 3) +
  scale_fill_manual(values = algorithm_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Percent valid estimates")
print(plt)


# Figure 2? --------------------------------------------------------------------

bulk_dataset <- "Mayo"
algorithm <- "Music"
err_files <- Get_ErrorFiles(bulk_dataset, algorithm, granularity)

marker_summary <- lapply(err_files, function(EF) {
  err_list <- readRDS(EF)

  if (length(err_list) == 0) {
    next
  }

  errs_all <- err_list$means$all_signature %>%
    subset(tissue != "All") %>%
    merge(err_list$params, by.x = "param_id", by.y = "row.names") %>%
    mutate(marker_combo = paste(marker_type, marker_subtype, marker_input_type),
           marker_combo = str_replace_all(marker_combo, " None", ""),
           algorithm = algorithm) %>%
    Paper_Renames() %>%
    merge(best_dt_broad, by = colnames(best_dt_broad))

  if (nrow(errs_all) == 0) {
    return(NULL)
  }

  errs_all$n_marker_type <- sapply(errs_all$n_markers, function(N) {
    if (N <= 1) {
      return("percent")
    }
    return("fixed")
  })

  errs_summary <- errs_all %>%
    group_by(tissue, tissue_full, data_transform, marker_combo, marker_order,
             n_marker_type, total_markers_used) %>%
    dplyr::summarize(Correlation = max(Correlation),
                     RMSE = min(RMSE),
                     MAPE = min(MAPE),
                     .groups = "drop")
  return(errs_summary)
})

marker_summary <- do.call(rbind, marker_summary)

markers_plot <- best_errors_broad



dev.off()

# Draft Main Figure 2 ----------------------------------------------------------

row1 <- (plt2Av2a + ggtitle("A")) + (plt2Av2b + labs(caption = NULL))
row2 <- (plt2Bv3a + ggtitle("B")) + (plt2Bv3b + labs(caption = NULL)) +
  ggplot() + ggtitle("C: Placeholder for ??") + plot_layout(widths = c(1, 1, 3))
row3 <- ggplot() + ggtitle("D: Placeholder for line graph", subtitle = "Best errors for each algorithm")
row4 <- ggplot() + ggtitle("E: Placeholder for valid results graph") +
  ggplot() + ggtitle("F: Placeholder for Exc:Inh ratio graph")

pdf(file.path(dir_figures,
              str_glue("draft_figure2.pdf")),
    width=7.5, height = 10)

print((row1) / (row2) / (row3) / (row4) + #plot_layout(heights = c(2, 2, 1, 1)) +
        plot_annotation("Draft main figure 2 -- please ignore weird formatting / sizing for now"))

dev.off()


# TODO left off here -- I'm thinking Fig 1 is general stats, Fig 2 is comparison
# of algorithms across broad/sub class, and Fig 3 is significance

# Figure 3? --------------------------------------------------------------------

significance <- readRDS(file.path(dir_analysis,
                                  str_glue("significance_lists_{granularity}.rds")))

sig_toplevel <- Paper_Renames(do.call(rbind, significance$significance_props_toplevel))

# Fill in missing data with NA
params <- expand.grid(celltype = unique(sig_toplevel$celltype),
                      tissue = unique(sig_toplevel$tissue),
                      algorithm = unique(sig_toplevel$algorithm),
                      normalization = unique(sig_toplevel$normalization),
                      regression_method = unique(sig_toplevel$regression_method))

sig_toplevel <- merge(sig_toplevel, params,
                      by = colnames(params),
                      all = TRUE)
sig_toplevel$p_adj_thresh[is.na(sig_toplevel$p_adj_thresh)] <- 1

sig_toplevel$data_transform <- paste(sig_toplevel$normalization, "+",
                                     sig_toplevel$regression_method)

if (granularity == "sub_class") {
  ct_order <- c("Astrocyte", "Endothelial", paste0("Exc.", 1:10),
                paste0("Inh.", 1:7), "Microglia", "Oligodendrocyte", "OPC",
                "Pericyte", "VLMC")
  sig_toplevel$celltype <- factor(sig_toplevel$celltype, levels = ct_order)
}

sig_final <- merge(sig_toplevel, best_dt,
                   by = c("tissue", "algorithm", "data_transform"),
                   all.x = FALSE)

# Cap log2_fc values to +/- 1 so color scaling is better. Need to account for
# Inf and -Inf values. Set non-significant log2 values to NA.
sig_final$log2_fc[sig_final$fc == 0] <- NA
sig_final$log2_fc[is.infinite(sig_final$log2_fc)] <- 1
sig_final$log2_fc[sig_final$log2_fc > 1] <- 1
sig_final$log2_fc[sig_final$log2_fc < -1] <- -1

sig_final$log2_fc[sig_final$p_adj_thresh >= 0.01] <- NA

plot_limit <- c(-1, 1)
n_cols <- if (granularity == "broad_class") 4 else 6

cap <- paste("Squares with color are significant at p < 0.01. Intensity of color\n",
             "is log2 fold change. I can also make this graph for any combination\n",
             "of normalization/regression for supplementary, in case we want to\n",
             "show that the same patterns show up no matter how we transform the data.\n",
             "V2 of this graph is on the next page.")
plt5v1 <- ggplot(sig_final, aes(x = tissue, y = algorithm, fill = log2_fc)) +
  geom_tile(color = "black") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.width = unit(10, "pt")) +
  facet_wrap(~celltype, ncol = n_cols) +
  coord_fixed() +
  scale_fill_distiller(palette = "RdBu", limit = plot_limit, na.value = "white") +
  xlab(NULL) + ylab(NULL) +
  labs(fill = "log2-FC", caption = cap) +
  ggtitle("Figure 5v1: Significant changes at p < 0.01 (using best data for each tissue)")
print(plt5v1)

ok <- subset(sig_final, is.finite(log2_fc))
sig_final_ok <- subset(sig_final, celltype %in% unique(ok$celltype))

cap <- paste("Alternately only show cell types that have significant differences,\n",
             "or we could cut down even further to cell types that have > N differences.")
plt5v2 <- ggplot(sig_final_ok, aes(x = tissue, y = algorithm, fill = log2_fc)) +
  geom_tile(color = "black") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.width = unit(10, "pt")) +
  facet_wrap(~celltype, ncol = 6) +
  coord_fixed() +
  scale_fill_distiller(palette = "RdBu", limit = plot_limit, na.value = "white") +
  xlab(NULL) + ylab(NULL) +
  labs(fill = "log2-FC", caption = cap) +
  ggtitle("Figure 5v2: Significant changes at p < 0.01 (using best data for each tissue)")
print(plt5v2)

dev.off()

# TODO PDF for supplemental

# Supplementary Figure 1a and 1b -----------------------------------------------

baselines_plot_tissue <- baselines_random %>%
  Create_BoxStats(c("tissue", "error_type"))

errs_box_noalg <- best_errors_broad %>%
  Create_BoxStats(c("tissue", "error_type",
                    "normalization", "regression_method"))

cap <- paste("Not sure if it's better to have everything on the same scale or\n",
             "have each plot on its own scale due to TCX scaling.")
plts1a <- ggplot(subset(errs_box_noalg, error_type == "rMSE"),
                 aes(x = normalization, fill = regression_method,
                     y = median_val, ymin = min_val, ymax = max_val)) +
  geom_hline(aes(yintercept = min_val, color = "Baseline"),
             data = subset(baselines_plot_tissue, error_type == "rMSE"),
             linetype = "twodash") +
  geom_crossbar(position = position_dodge2(padding = 0),
                fatten = 1.5, width = 0.75) +
  theme_bw() +
  scale_fill_manual(values = regression_colors) +
  scale_color_manual(values = "darkgray") +
  facet_wrap(~tissue, nrow = 3, scales = "fixed") +
  labs(fill = "Regression Method", color = "", caption = cap) +
  xlab("Normalization Strategy") +
  ylab("RMSE") +
  ggtitle("Supplementary Figure S1a: RMSE: normalization vs regression")
print(plts1a / plot_spacer() + plot_layout(heights = c(3, 1)))

plts1b <- ggplot(subset(errs_box_noalg, error_type == "mAPE"),
                 aes(x = normalization, fill = regression_method,
                     y = median_val, ymin = min_val, ymax = max_val)) +
  geom_hline(aes(yintercept = min_val, color = "Baseline"),
             data = subset(baselines_plot_tissue, error_type == "mAPE"),
             linetype = "twodash") +
  geom_crossbar(position = position_dodge2(padding = 0),
                fatten = 1.5, width = 0.75) +
  theme_bw() +
  scale_fill_manual(values = regression_colors) +
  scale_color_manual(values = "darkgray") +
  facet_wrap(~tissue, nrow = 3, scales = "fixed") +
  labs(fill = "Regression Method", color = "") +
  xlab("Normalization Strategy") +
  ylab("MAPE") +
  ggtitle("Supplementary Figure S1b: MAPE: normalization vs regression")
print(plts1b / plot_spacer() + plot_layout(heights = c(3, 1)))


# Supplementary Figure 1c ------------------------------------------------------

baselines_plot <- Create_BoxStats(baselines_random,
                                  grouping_cols = c("test_data_name", "algorithm",
                                                    "normalization", "regression_method",
                                                    "error_type"))
errs_box <- Create_BoxStats(best_errors_broad,
                            grouping_cols = c("test_data_name", "algorithm",
                                              "normalization", "regression_method",
                                              "error_type")) %>%
  rbind(baselines_plot)

cap <- paste("Errors broken out by algorithm, including baseline. For readability,\n",
             "Tissues have been collapsed to study level.")
plts1cv1 <- ggplot(subset(errs_box, error_type == "cor"),
                 aes(x = normalization, fill = regression_method,
                     y = median_val, ymin = min_val, ymax = max_val)) +
  geom_crossbar(position = position_dodge2(padding = 0),
                fatten = 1.5, width = 0.75) +
  theme_bw() +
  scale_fill_manual(values = regression_colors) +
  scale_color_manual(values = "darkgray") +
  facet_grid(test_data_name ~ algorithm, scales = "fixed") +
  labs(fill = "Regression Method", color = "", caption = cap) +
  xlab("Normalization Strategy") +
  ylab("Correlation") +
  ggtitle("Supplementary Figure S1c v1: Correlation: normalization vs regression by algorithm")
print(plts1cv1 / plot_spacer() + plot_layout(heights = c(3, 1)))

baselines_plot_tissue <- Create_BoxStats(baselines_random,
                                         grouping_cols = c("tissue", "algorithm",
                                                           "normalization", "regression_method",
                                                           "error_type"))
errs_box_best <- Create_BoxStats(best_errors_broad,
                                 grouping_cols = c("tissue", "algorithm",
                                                   "normalization", "regression_method",
                                                   "error_type")) %>%
  rbind(baselines_plot_tissue) %>%
  mutate(tissue_display = tissue,
         tissue = str_replace(tissue, ".*\ ", ""),
         data_transform = paste(normalization, "+", regression_method)) %>%
  merge(best_dt, by = c("tissue", "algorithm", "data_transform"), all.x = FALSE)

cap <- paste("Or we can just use the best norm/regressions for each tissue and\n",
             "put a dashed line for baseline.")
plts1cv2 <- ggplot(subset(errs_box_best, error_type == "cor" & algorithm != "Baseline"),
                   aes(x = algorithm, fill = algorithm,
                       y = median_val, ymin = min_val, ymax = max_val)) +
  geom_hline(aes(yintercept = max_val, color = "Baseline"),
             data = subset(errs_box_best, error_type == "cor" & algorithm == "Baseline"),
             linetype = "twodash") +
  geom_crossbar(fatten = 1.5, width = 0.5) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = algorithm_colors) +
  scale_color_manual(values = "darkgray") +
  facet_wrap(~tissue_display, nrow = 3, scales = "fixed") +
  labs(fill = "Regression Method", color = "", caption = cap) +
  xlab("Normalization Strategy") +
  ylab("Correlation") +
  ggtitle("Supplementary Figure S1c v2")
#print(plts1cv2 / plot_spacer() + plot_layout(heights = c(3, 1)))

cap <- paste("Or we can use the best norm/regressions and put it all in one\n",
             "plot, but I need to fix this data to pull in more than one \n", "
             baseline value per tissue. Repeat for RMSE and MAPE.")
plts1cv3 <- ggplot(subset(errs_box_best, error_type == "cor"), # & algorithm != "Baseline"),
                 aes(x = tissue_display, fill = algorithm,
                     y = median_val, ymin = min_val, ymax = max_val)) +
  geom_crossbar(position = position_dodge2(padding = 0),
                fatten = 1.5, width = 0.75) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c(algorithm_colors, "Baseline" = "darkgray")) +
  scale_color_manual(values = "darkgray") +
  labs(fill = "Regression Method", color = "", caption = cap) +
  xlab("Normalization Strategy") +
  ylab("Correlation") +
  ggtitle("Supplementary Figure S1c v3")
print(plts1cv2 / plts1cv3)


# Supplmentary Figure 2a and 2b ------------------------------------------------

# Quality stats, TODO


# Supplementary Figure 3a, 3b --------------------------------------------------

plts3a <- ggplot(Count_Grouped(ranked_df_broad, c("algorithm", "error_type")),
                 aes(x = algorithm, y = error_type, color = count, size = count)) +
  geom_count() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5) +
  labs(color = "Count", size = "Count") +
  xlab("Algorithm") +
  ylab("Error Metric") +
  ggtitle("Supplementary Figure S3a", subtitle = "Best algorithm per error metric") +
  guides(size = "none")
#print(plt)

cap <- "Not sure how useful/interesting these graphs are."
plts3b <- ggplot(Count_Grouped(ranked_df_broad, c("error_type", "normalization", "regression_method")),
                 aes(x = normalization, y = regression_method, color = count, size = count)) +
  geom_count() + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5) +
  facet_wrap(~error_type, nrow = 1) +
  labs(color = "Count", size = "Count", caption = cap) +
  xlab("Normalization") +
  ylab("Regression Method") +
  ggtitle("Supplementary Figure S3b", subtitle = "Best normalization/regression per error type") +
  guides(size = "none")
#print(plt)
print((plts3a + plts3b) / plot_spacer())

dev.off()

