library(Matrix)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(dplyr)
library(stringr)
library(reshape2)
library(patchwork)

source(file.path("functions", "FileIO_HelperFunctions.R"))
source(file.path("functions", "Step15_Analysis_HelperFunctions.R"))
source(file.path("functions", "Step15_Plotting_HelperFunctions.R"))

options(scipen = 999)

bulk_datasets <- c("Mayo", "MSBB", "ROSMAP")

best_errors_list_broad <- Load_BestErrors("broad_class")
best_errors_list_sub <- Load_BestErrors("sub_class")
quality_stats_broad <- Load_QualityStats(bulk_datasets, "broad_class")
quality_stats_sub <- Load_QualityStats(bulk_datasets, "sub_class")

# Unpack error variables into environment for readability
list2env(best_errors_list_broad, globalenv())
list2env(best_errors_list_sub, globalenv())

tissues_use <- c("TCX", "PHG", "ACC") # When subsetting
tissues_use_full <- c("Mayo TCX", "MSBB PHG", "ROSMAP ACC") # When subsetting


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
  Create_BoxStats(c("tissue_full", "error_metric")) %>%
  group_by(error_metric) %>%
  mutate(best_val = if (unique(error_metric) == "Correlation") max_val else min_val) %>%
  ungroup()

errs_box_noalg <- best_errors_broad %>%
  Create_BoxStats(c("tissue_full", "error_metric",
                    "normalization", "regression_method"))

cap <- "Show correlation only, for all tissues. RMSE and MAPE graphs would be in supplementary."
plt1Bv1 <- ggplot(subset(errs_box_noalg, error_metric == "Correlation"),
                  aes(x = normalization, fill = regression_method,
                      y = median_val, ymin = min_val, ymax = max_val)) +
  geom_hline(aes(yintercept = best_val, color = "Baseline"),
             data = subset(baselines_plot_tissue, error_metric == "Correlation"),
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
plt1Bv2 <- ggplot(subset(errs_box_noalg, error_metric == "Correlation" &
                           tissue_full %in% tissues_use_full),
                  aes(x = normalization, fill = regression_method,
                      y = median_val, ymin = min_val, ymax = max_val)) +
  geom_hline(aes(yintercept = best_val, color = "Baseline"),
             data = subset(baselines_plot_tissue, error_metric == "Correlation" &
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
  facet_grid(error_metric ~ tissue_full, scales = "free", switch = "y") +
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
             "few tissues, so these could possibly go in supplemental if the main\n",
             "figure is pooled data.")

ranked_df_all <- rbind(top3_by_tissue_broad, top3_by_tissue_sub)

plt1Cav1 <- ggplot(Count_Grouped(top3_by_tissue_broad, c("tissue", "normalization", "regression_method")),
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
  ggtitle("Figure 1C v1: Best normalization/regression per tissue", subtitle = "Broad class")

plt1Cbv1 <- ggplot(Count_Grouped(ranked_df_sub, c("tissue", "normalization", "regression_method")),
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

plt1Cv1 <- plt1Cav1 / plt1Cbv1

cap <- paste("Original. Here the counts for broad and sub class are pooled into one graph.\n",
             "Need to hand-fix the legend to be whole numbers.")
plt1Cv2 <- ggplot(Count_Grouped(ranked_df_all, c("tissue", "normalization", "regression_method")),
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
  ggtitle("Plot 1C v2: Pooled broad and sub class")

print(plt1Cv1 / plt1Cv2 / plot_spacer() + plot_layout(heights = c(1, 1, 1, 1)))


# Figure 1D --------------------------------------------------------------------

cap <- paste("Took the top 3 scoring estimates for each tissue and error metric\n",
             "(9 total per tissue) and counted how many times each reference\n",
             "appeared. This can also be changed to a 3x3 grid. Broad class and\n",
             "sub class look nearly identical, should sub class go in supplemental?")
plt1Dv1a <- ggplot(Count_Grouped(top3_by_tissue_broad, c("tissue", "reference_data_name")),
                   aes(x = tissue, y = reference_data_name, color = count, size = count)) +
  geom_count() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.width = unit(10, "pt")) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5) +
  labs(color = "Count") +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Figure 1D v1: Best reference", subtitle = "Broad class") +
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

ranked_all <- rbind(top3_by_tissue_broad, ranked_df_sub)

cap <- "Pooled results for broad and sub class"
plt1Dv2 <- ggplot(Count_Grouped(ranked_all, c("tissue", "reference_data_name")),
                  aes(x = tissue, y = reference_data_name, color = count, size = count)) +
  geom_count() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.width = unit(10, "pt")) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5) +
  labs(color = "Count", caption = cap) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Figure 1D v2: Best reference", subtitle = "Broad class") +
  guides(size = "none")

print(plt1Dv1a + plt1Dv1b)
print(plt1Dv2)


# Figure 1E --------------------------------------------------------------------

# TODO calculate this in Step 15
total_valid_broad <- quality_stats_broad$n_valid_by_algorithm %>%
  subset(algorithm != "Baseline") %>%
  group_by(algorithm) %>%
  summarize(n_valid = sum(n_valid),
            n_possible = sum(n_possible),
            pct_valid = n_valid / n_possible,
            granularity = "Broad class")

total_valid_sub <- quality_stats_all_sub %>%
  subset(algorithm != "Baseline") %>%
  group_by(algorithm) %>%
  summarize(n_valid_results = sum(n_valid_results),
            n_possible_results = sum(n_possible_results),
            pct_valid_results = n_valid_results / n_possible_results,
            granularity = "Sub class")

# So something shows up for DeconRNASeq
total_valid_sub$pct_valid_results[total_valid_sub$pct_valid_results < 0.01] <- 0.01

cap = "Graphs separated by class"
plt1Ev1a <- ggplot(total_valid_broad,
                   aes(x = algorithm, y = pct_valid, fill = algorithm)) +
  geom_col() + theme_bw() +
  scale_fill_manual(values = algorithm_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Figure 1E v1: Percent results passing QC", subtitle = "Broad class") +
  xlab(NULL) +
  ylab("Percent") +
  guides(fill = "none")

plt1Ev1b <- ggplot(total_valid_sub,
                   aes(x = algorithm, y = pct_valid_results, fill = algorithm)) +
  geom_col() + theme_bw() +
  scale_fill_manual(values = algorithm_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(caption = cap) +
  ggtitle(NULL, subtitle = "Sub class") +
  xlab(NULL) +
  ylab("Percent") +
  guides(fill = "none")

plt1Ev1 <- plt1Ev1a + plt1Ev1b

total_valid_combined <- rbind(total_valid_broad, total_valid_sub)

cap <- paste("Original: Both classes in the same graph. Breaking out\n",
             "by tissue or norm/regression doesn't provide more\n",
             "information")
plt1Ev2 <- ggplot(total_valid_combined,
                  aes(x = algorithm, y = pct_valid_results, fill = granularity)) +
  geom_col(position = "dodge2") + theme_bw() +
  scale_fill_manual(values = regression_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(fill = NULL, caption = cap) +
  ggtitle("Figure 1E v2: Percent results passing QC") +
  xlab(NULL) +
  ylab("Percent")

print(plt1Ev1 / (plt1Ev2 + plot_spacer()) / plot_spacer() + plot_layout(heights = c(1, 1, 3)))

dev.off()


# Draft Main Figure 1 ----------------------------------------------------------

row1 <- ggplot() + ggtitle("A: Placeholder for diagram of process")
row2 <- plt1Bv3 + ggtitle("B") + labs(caption = NULL) +
  theme(legend.title = element_text(size = 10))
row3 <- plt1Cv2 + ggtitle("C", subtitle = NULL) + labs(caption = NULL)
row4 <- (plt1Dv2 + ggtitle("D", subtitle = NULL) + labs(caption = NULL)) +
  (plt1Ev2 + ggtitle("E") + labs(caption = NULL))

cap <- paste(
  "Figure 1. Deconvolution performance across all tests. A) Diagram of our",
  "deconvolution benchmarking pipeline.",
  "B) The average \ncorrelation, RMSE, and MAPE for representative tissues of the",
  "best-performing estimates, across different normalization and \nregression strategies.",
  "C) Count of how many times each normalization and regression pairing was used",
  "in the top three best-\nperforming estimates by each of the three error metrics,",
  "for each tissue (maximum of 18 estimates per tissue, 9 each for broad \nclass and",
  "sub class cell type estimates).",
  "D) Count of how many times each single cell reference was used as input for the",
  "top \nthree best-performing estimates per error metric and tissue.",
  "E) Percentage of estimates from each algorithm that passed quality \ncontrol, across",
  "all tests of broad and sub class cell types.",
  "<not sure if we need to define acronyms here> CPM = counts per \nmillion, TMM =",
  "trimmed mean of M values, TPM = transcripts per million. CBE = cerebellum,",
  "TCX = temporal cortex, FP = frontal \npole, IFG = inferior frontal gyrus, PHG =",
  "parahippocampal gyrus, STG = superior temporal gyrus, ACC = anterior cingulate",
  "cortex, \nDLPFC = dorsolateral prefrontal cortex, PCC = posterior cingulate cortex."
)

main_fig1 <- (row1) / (row2) / (row3) / (row4) + plot_layout(heights = c(2, 2, 1, 1)) +
        plot_annotation("Draft main figure 1 -- please ignore weird formatting / sizing for now",
                        caption = cap)


# Figure 2 variations ----------------------------------------------------------

pdf(file.path(dir_figures,
              str_glue("draft_figure2_variations.pdf")),
    width=7.5, height = 10)


# Figure 2A --------------------------------------------------------------------

cap <- paste("Took the top 3 scoring estimates for each tissue and error metric\n",
             "(9 total per tissue) and counted how many times each algorithm\n",
             "appeared. This can also be changed to a 3x3 grid.")

plt2Av1a <- ggplot(Count_Grouped(top3_by_tissue_broad, c("tissue", "algorithm", "error_type")),
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
alg_info <- Count_Grouped(top3_by_tissue_broad, c("tissue", "algorithm", "error_type")) %>%
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

cap <- paste("Original: get rid of error metric to make the graph simpler.\n",
             "Thematically this fits better in Figure 1 but there wasn't room.")
plt2Av2a <- ggplot(Count_Grouped(top3_by_tissue_broad, c("tissue", "algorithm")),
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
alg_info <- Count_Grouped(top3_by_tissue_broad, c("tissue", "algorithm")) %>%
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

# TODO none of this is right

baselines_plot <- baselines_broad %>%
  subset(reference_data_name != "All zeros") %>%
  Create_BoxStats(grouping_cols = c("tissue_full", "normalization",
                                    "regression_method", "error_metric"))

errs_better <- merge(best_errors_broad, baselines_plot) %>%
  mutate(better = case_when(error_metric == "Correlation" ~ value > max_val,
                            TRUE ~ value < min_val))

better_stats <- errs_better %>%
  group_by(tissue_full, algorithm) %>%
  dplyr::summarize(count = n(),
                   pct_better_than_baseline = sum(better) / count,
                   .groups = "drop")

baselines_plot_filt <- baselines_broad %>%
  subset(reference_data_name != "All zeros") %>%
  merge(best_dt_broad) %>%
  Create_BoxStats(grouping_cols = c("tissue_full", "normalization",
                                    "regression_method", "error_metric"))

errs_better_filt <- best_errors_top_broad %>%
  merge(best_dt_broad) %>%
  merge(baselines_plot_filt) %>%
  mutate(better = case_when(error_metric == "Correlation" ~ value > max_val,
                            TRUE ~ value < min_val))

better_stats_filt <- errs_better_filt %>%
  group_by(tissue_full, algorithm) %>%
  dplyr::summarize(count = n(),
                   pct_better_than_baseline = sum(better) / count,
                   .groups = "drop")

baselines_plot_filt_sub <- baselines_top_sub %>%
  #subset(reference_data_name != "All zeros") %>%
  merge(best_dt_broad) %>%
  Create_BoxStats(grouping_cols = c("tissue_full", "normalization",
                                    "regression_method", "error_metric"))

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

cap <- paste("Showing 3 representative tissues and putting \nthe rest in supplemental",
             "figures. This graph \nis for broad class. We could also wedge\n",
             "sub class next to it.")

plt2Bv2 <- ggplot(subset(better_stats_filt, tissue_full %in% tissues_use_full),
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
  ggtitle(NULL)

better_stats_filt2 <- errs_better_filt %>%
  group_by(algorithm) %>%
  dplyr::summarize(count = n(),
                   pct_better_than_baseline = sum(better) / count,
                   granularity = "Broad class",
                   .groups = "drop")

better_stats_filt_sub2 <- errs_better_filt_sub %>%
  group_by(algorithm) %>%
  dplyr::summarize(count = n(),
                   pct_better_than_baseline = sum(better) / n(),
                   granularity = "Sub class",
                   .groups = "drop")
better_stats_filt_sub2$pct_better_than_baseline[better_stats_filt_sub2$pct_better_than_baseline == 0] <- 0.01

cap <- paste("All data is collapsed across all tissues, and we \ncan put the",
             "tissue-specific breakout \nin supplemental.")
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

better_stats_filt_combined <- rbind(better_stats_filt2, better_stats_filt_sub2)

cap <- "Original. Classes combined in the same graph."
plt2Bv4 <-ggplot(better_stats_filt_combined,
                 aes(x = algorithm, y = pct_better_than_baseline, fill = granularity)) +
  geom_col(position = "dodge2") +
  scale_fill_manual(values = regression_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(fill = NULL, caption = cap) +
  ylim(0, 1.05) +
  xlab(NULL) +
  ylab("Percent") +
  ggtitle("Figure 2B v4")

print((plt2Bv1 + plt2Bv2 + plot_layout(widths = c(3, 2))) /
        (plt2Bv3 + plt2Bv4 + plot_layout(widths = c(1, 1, 3))) + plot_layout(heights = c(3, 1)))


# Figure 2? --------------------------------------------------------------------

print(ggplot() + ggtitle("Placeholder: two other figures in progress"))

# TODO redo when median values are ready
inh_mean <- quality_stats_top_broad %>%
  merge(best_dt_broad, by = c("tissue", "data_transform", "algorithm"), all = FALSE) %>%
  group_by(tissue, algorithm) %>%
  summarize(median_ratio = median(mean_ratio),
            mean_ratio = mean(mean_ratio),
            .groups = "drop")

# Cap to remove large values. We'll manually put a break in or something...
inh_mean$mean_ratio[inh_mean$mean_ratio > 20] <- 20

cap <- paste("Mean ratio of excitatory to inhibitory neurons for broad class\n",
             "estimates, for the top 3 estimates / best norm/regression. Values\n",
             "larger than 20 have been capped so the graph is readable: Mean\n",
             "values range from 2.37 to 2.76e17 :( and 11 of the bars in this\n",
             "graph would go above 100 without capping. Median values are similar.\n",
             "DWLS is the biggest offender, followed by MuSiC.\n",
             "This particular data isn't encouraging and shows that inhibitory\n",
             "neurons are being way under-estimated and likely have values close\n",
             "to zero in a lot of cases.")
plt <- ggplot(inh_mean, aes(x = algorithm, y = mean_ratio,
                            fill = algorithm)) +
  geom_col() + theme_bw() +
  facet_wrap(~tissue, nrow = 3) +
  scale_fill_manual(values = algorithm_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(caption = cap) +
  ggtitle("Mean ratio of Excitatory:Inhibitory neurons") +
  xlab(NULL) +
  ylab("Ratio") +
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

bulk_dataset <- "Mayo"
algorithm <- "Music"
err_files <- Get_ErrorFiles(bulk_dataset, algorithm, granularity)

marker_summary <- lapply(err_files, function(EF) {
  err_list <- readRDS(EF)

  if (length(err_list) == 0) {
    next
  }

  errs_all <- err_list$means %>%
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

# TODO

dev.off()


# Draft Main Figure 2 ----------------------------------------------------------

row1 <- (plt2Av2a + ggtitle("A")) + (plt2Av2b + labs(caption = NULL))
row2 <- (plt2Bv4 + ggtitle("B") + labs(caption = NULL)) +
  ggplot() + ggtitle("C: Placeholder for ??") + plot_layout(widths = c(1,1))
row3 <- ggplot() + ggtitle("D: Placeholder for line graph",
                           subtitle = "Best errors for each algorithm vs tissue") +
  ggplot() + ggtitle("E: Placeholder for Exc:Inh ratio graph")

cap <- paste(
  "Figure 2. Comparison of algorithm quality. A) Count of how many times each algorithm",
  "appears in the top three best-performing \nestimates per error metric and tissue",
  "(maximum of 9 estimates per tissue, 3 error metrics with 3 estimates each). B)",
  "Percent of \nestimates from each algorithm that out-perform the baseline error rate.",
  "Data is filtered to the top three best-performing estimates \nper algorithm per",
  "tissue, one estimate per error metric. C) ?? we can think of something to put here.",
  "D) comparison of the best \ncorrelation between algorithms across tissues. E) Median",
  "ratio of excitatory to inhibitory neurons in the top three estimates \nfrom each",
  "algorithm. <We might have room for another row of small graphs if necessary.>"
)

main_fig2 <- (row1) / (row2) / (row3) +
        plot_annotation("Draft main figure 2 -- please ignore weird formatting / sizing for now",
                        caption = cap)


# TODO I'm thinking Fig 1 is general stats, Fig 2 is comparison of algorithms
# across broad/sub class, and Fig 3 is significance + some bar graphs of
# individual estimates

# Figure 3 variations ----------------------------------------------------------

pdf(file.path(dir_figures,
              str_glue("draft_figure3_variations.pdf")),
    width=7.5, height = 10)


# Figure 3A --------------------------------------------------------------------

significance_broad <- Load_Significance(quality_stats_broad$significance_toplevel,
                                        best_dt_broad,
                                        p_sig = 0.01, log2_cap = 1)

plot_limit <- c(-1, 1)
n_cols <- 4 #if (granularity == "broad_class") 4 else 6

cap <- paste("Squares with color are significant at p < 0.01. Intensity of color",
             "is log2 fold change. This graph uses the \nbest norm/regression for",
             "each tissue, with significance calculated on the average of the top 3",
             "estimates \nfor each algorithm. I can also make this graph for any combination",
             "of normalization/regression for \nsupplementary, in case we want to",
             "show that the same patterns show up across data transforms.")
plt3Av1 <- ggplot(significance_broad, aes(x = tissue, y = algorithm, fill = log2_fc)) +
  geom_tile(color = "black") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.width = unit(10, "pt")) +
  facet_wrap(~celltype, ncol = n_cols) +
  coord_fixed() +
  scale_fill_distiller(palette = "RdBu", limit = plot_limit, na.value = "white") +
  xlab(NULL) + ylab(NULL) +
  labs(fill = "log2-FC", caption = cap) +
  ggtitle("Figure 3A v1: Broad class significant changes at p < 0.01")

ok <- subset(significance_broad, is.finite(log2_fc))
sig_final_ok <- subset(significance_broad, celltype %in% unique(ok$celltype))

cap <- paste("Only cell types with at least one significant value. Alternately we could\n",
             "cut down even further to cell types that have >= 3 differences, as in v3.")
plt3Av2 <- ggplot(sig_final_ok, aes(x = tissue, y = algorithm, fill = log2_fc)) +
  geom_tile(color = "black") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.width = unit(10, "pt")) +
  facet_wrap(~celltype, ncol = 6) +
  coord_fixed() +
  scale_fill_distiller(palette = "RdBu", limit = plot_limit, na.value = "white") +
  xlab(NULL) + ylab(NULL) +
  labs(fill = "log2-FC", caption = cap) +
  ggtitle("Figure 3A v2")

ok2 <- significance_broad %>%
  group_by(celltype) %>%
  summarize(n_significant = sum(is.finite(log2_fc))) %>%
  subset(n_significant >= 3)

sig_final_ok2 <- subset(significance_broad, celltype %in% ok2$celltype)

# Main figure but not graphed here
plt3Av3 <- ggplot(sig_final_ok2, aes(x = tissue, y = algorithm, fill = log2_fc)) +
  geom_tile(color = "black") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.width = unit(10, "pt")) +
  facet_wrap(~celltype, ncol = 4) +
  coord_fixed() +
  scale_fill_distiller(palette = "RdBu", limit = plot_limit, na.value = "white") +
  xlab(NULL) + ylab(NULL) +
  labs(fill = "log2-FC") +
  ggtitle("Figure 3A v3")

print(plt3Av1 / plt3Av2 / plt3Av3)


# Figure 3B --------------------------------------------------------------------

significance_sub <- Load_Significance("sub_class", p_sig = 0.01, log2_cap = 1)

plot_limit <- c(-1, 1)
n_cols <- 6

cap <- paste("Same idea as broad class.")
plt3Bv1 <- ggplot(significance_sub, aes(x = tissue, y = algorithm, fill = log2_fc)) +
  geom_tile(color = "black") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.width = unit(10, "pt")) +
  facet_wrap(~celltype, ncol = n_cols) +
  coord_fixed() +
  scale_fill_distiller(palette = "RdBu", limit = plot_limit, na.value = "white") +
  xlab(NULL) + ylab(NULL) +
  labs(fill = "log2-FC", caption = cap) +
  ggtitle("Figure 3B v1: Sub class significant changes at p < 0.01")

ok <- subset(significance_sub, is.finite(log2_fc))
sig_final_ok <- subset(significance_sub, celltype %in% unique(ok$celltype))

# Not shown
cap <- paste("Only cell types with at least one significant value.")
plt3Bv2 <- ggplot(sig_final_ok, aes(x = tissue, y = algorithm, fill = log2_fc)) +
  geom_tile(color = "black") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.width = unit(10, "pt")) +
  facet_wrap(~celltype, ncol = 6) +
  coord_fixed() +
  scale_fill_distiller(palette = "RdBu", limit = plot_limit, na.value = "white") +
  xlab(NULL) + ylab(NULL) +
  labs(fill = "log2-FC", caption = cap) +
  ggtitle("Figure 3A v2")

ok2 <- significance_sub %>%
  group_by(celltype) %>%
  summarize(n_significant = sum(is.finite(log2_fc))) %>%
  subset(n_significant >= 3)

sig_final_ok2 <- subset(significance_sub, celltype %in% ok2$celltype)

cap <- paste("Cell types with >= 3 significant values. If we show everything with any\n",
             "significant values, we get 16 of 24 cell types which is a lot.")
plt3Bv3 <- ggplot(sig_final_ok2, aes(x = tissue, y = algorithm, fill = log2_fc)) +
  geom_tile(color = "black") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.width = unit(10, "pt")) +
  facet_wrap(~celltype, ncol = 5) +
  coord_fixed() +
  scale_fill_distiller(palette = "RdBu", limit = plot_limit, na.value = "white") +
  labs(caption = cap) +
  xlab(NULL) + ylab(NULL) +
  labs(fill = "log2-FC") +
  ggtitle("Figure 3B v2")

print(plt3Bv1 / plt3Bv3)

dev.off()


# Draft Main Figure 3 ----------------------------------------------------------

row1 <- (plt3Av3 + ggtitle("A") + labs(caption = NULL))
row2 <- (plt3Bv3 + ggtitle("B") + labs(caption = NULL))
row3 <- ggplot() + ggtitle("C: Placeholder for some individual estimate stuff")

cap <- paste(
  "Figure 3. Consensus between algorithms shows significant AD-related changes in",
  "multiple cell types and tissues. A, B) Significant \ndifferences in broad class (A)",
  "and sub class (B) cell types between AD and control samples. Squares are colored in",
  "if the difference \nis significant at p < 0.01. Only cell types with more than three",
  "significant values are shown. C) Range of individual estimates for",
  "\n<select algorithms>. This graph will show that they're mostly reasonable estimates."
)

main_fig3 <- (row1) / (row2) / (row3) +
        plot_annotation("Draft main figure 3 -- please ignore weird formatting / sizing for now",
                        caption = cap)

pdf(file.path(dir_figures,
              str_glue("draft_main_figures.pdf")),
    width=7.5, height = 10)

print(main_fig1)
print(main_fig2)
print(main_fig3)

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

