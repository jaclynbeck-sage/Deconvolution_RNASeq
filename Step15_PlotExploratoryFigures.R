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

# TODO remove the Mayo TCX param with > 800 rMSE? it's all from TPM

options(scipen = 999)

datasets <- c("cain", "lau", "leng", "mathys", "seaRef")

algorithms <- c("CibersortX", "DeconRNASeq", "Dtangle", "DWLS", "HSPE", "Music",
                "Scaden", "Baseline")

granularity <- c("broad_class")

bulk_datasets <- c("Mayo", "MSBB", "ROSMAP")

best_errors_list <- readRDS(file.path(dir_analysis,
                                      str_glue("best_errors_{granularity}.rds")))
best_errors <- best_errors_list$best_errors_all

best_errs_plot <- Paper_Renames(best_errors) %>%
  mutate(tissue = tissue_full)

# There are up to 3 param_ids per set of input parameters. Get the best of each
# error metric for each set
errs_melt <- best_errs_plot  %>%
  group_by_at(c("tissue", Get_ParameterColumnNames())) %>%
  dplyr::summarize(Correlation = max(Correlation),
                   RMSE = min(RMSE),
                   MAPE = min(MAPE),
                   .groups = "drop") %>%
  melt(variable.name = "error_type")

baselines_melt <- subset(errs_melt, algorithm == "Baseline")
errs_melt <- subset(errs_melt, algorithm != "Baseline")

# Color setup ------------------------------------------------------------------

refs <- unique(errs_melt$reference_data_name)
algs <- unique(errs_melt$algorithm)
norms <- unique(str_replace(errs_melt$normalization, "log_", ""))
regs <- unique(errs_melt$regression_method)

# Reference data sets
reference_colors <- RColorBrewer::brewer.pal(length(refs), "Set3")
names(reference_colors) <- sort(refs)

# Algorithms
algorithm_colors <- RColorBrewer::brewer.pal(length(algs), "Set2")
names(algorithm_colors) <- sort(algs)

# Normalizations
normalization_colors <- RColorBrewer::brewer.pal(length(norms), "Set1")

# Regression methods
regression_colors <- viridis::turbo(length(regs), begin = 0.1, end = 0.9)#RColorBrewer::brewer.pal(length(regs), "Accent")

# Tissues (modified viridis turbo color scheme)
tiss <- colSums(table(errs_melt$tissue, errs_melt$test_data_name) > 0)
tissue_colors <- c(viridis::turbo(tiss[["Mayo"]], begin = 0.1, end = 0.2),
                   viridis::turbo(tiss[["MSBB"]], begin = 0.35, end = 0.6),
                   viridis::turbo(tiss[["ROSMAP"]], begin = 0.7, end = 0.9))
names(tissue_colors) <- sort(unique(errs_melt$tissue))

# MSBB colors need to be darker and a little more differentiated -- original
# colors are 20EAABFF, 67FD68FF, AEFA37FF, E1DD37FF
tissue_colors[grepl("MSBB", names(tissue_colors))] <- c("#00CA8BFF", "#37CD38FF",
                                                        "#7ECA07FF", "#C1BD17FF")

# Slightly more pastel than default looks better
tissue_fill_colors <- str_replace(tissue_colors, "FF$", "88")
names(tissue_fill_colors) <- names(tissue_colors)

# Bulk data sets -- match some of the tissue colors
bulk_colors <- tissue_fill_colors[c(1, 4, 9)]
names(bulk_colors) <- sort(bulk_datasets)

baselines_random <- subset(baselines_melt, reference_data_name != "All zeros")
baselines_zeros <- subset(baselines_melt,
                          reference_data_name == "All zeros" & error_type == "RMSE")

# Best baseline scores for each bulk data set, regardless of normalization / regression
baselines_plot <- baselines_random %>%
  Create_BoxStats(c("test_data_name", "error_type"))

zeros_plot <- baselines_zeros %>%
  Create_BoxStats(c("test_data_name", "error_type"))

pdf(file.path(dir_figures,
              str_glue("error_plots_{granularity}_summary.pdf")),
    width=10, height = 12)

# Spread of scores regardless of normalization / regression --------------------

errs_box <- Create_BoxStats(errs_melt,
                            c("tissue", "algorithm", "reference_data_name",
                              "test_data_name", "error_type"))

groups <- c("tissue", "algorithm", "test_data_name", "error_type")
errs_box_noref <- Create_BoxStats(errs_melt,
                                  c("tissue", "algorithm", "test_data_name",
                                    "error_type"))

# correlation
plt1 <- Plot_FacetBoxPlot(subset(errs_box, error_type == "Correlation"),
                          fill = "tissue",
                          x_axis = "reference_data_name",
                          fill_colors = tissue_colors,
                          facet_vars = c("test_data_name", "algorithm")) +
  geom_hline(aes(yintercept = max_val),
             data = subset(baselines_plot, error_type == "Correlation"),
             linetype = "twodash")
print(plt1 + plot_annotation(title = "Spread of correlation regardless of normalization / regression"))

plt2 <- Plot_FacetBoxPlot(subset(errs_box_noref, error_type == "Correlation"),
                          fill = "tissue",
                          x_axis = "algorithm",
                          fill_colors = tissue_colors,
                          facet_vars = "test_data_name") +
  geom_hline(aes(yintercept = max_val),
             data = subset(baselines_plot, error_type == "Correlation"),
             linetype = "twodash")
print(plt2 + plot_annotation(title = "Spread of correlation regardless of reference, normalization, or regression"))

# Plt 1 but divided by tissue so there is a baseline line for each tissue
# separately
baselines_plot_tissue <- baselines_random %>%
  Create_BoxStats(c("tissue", "error_type"))

plt3 <- Plot_FacetBoxPlot(subset(errs_box, error_type == "Correlation"),
                          fill = "tissue",
                          x_axis = "reference_data_name",
                          fill_colors = tissue_colors,
                          facet_vars = c("tissue", "algorithm")) +
  geom_hline(aes(yintercept = max_val),
             data = subset(baselines_plot_tissue, error_type == "Correlation"),
             linetype = "twodash")
print(plt3 + plot_annotation(title = "Spread of correlation by tissue, regardless of normalization, or regression"))


# rMSE -- this is the only one where the "zeros" baseline is relevant
plt4a <- Plot_FacetBoxPlot(subset(errs_box, error_type == "RMSE"),
                           fill = "tissue",
                           x_axis = "reference_data_name",
                           fill_colors = tissue_colors,
                           facet_vars = c("test_data_name", "algorithm")) +
  geom_hline(aes(yintercept = min_val),
             data = subset(baselines_plot, error_type == "RMSE"),
             linetype = "twodash") +
  geom_hline(aes(yintercept = min_val), data = zeros_plot, color = "red",
             linetype = "twodash")

print(plt4a + plot_annotation(title = "Spread of RMSE regardless of normalization, or regression"))

plt4b <- Plot_FacetBoxPlot(subset(errs_box_noref, error_type == "RMSE"),
                           fill = "tissue",
                           x_axis = "algorithm",
                           fill_colors = tissue_colors,
                           facet_vars = "test_data_name") +
  geom_hline(aes(yintercept = min_val),
             data = subset(baselines_plot, error_type == "RMSE"),
             linetype = "twodash") +
  geom_hline(aes(yintercept = min_val), data = zeros_plot, color = "red",
             linetype = "twodash")

print(plt4b + plot_annotation(title = "Spread of RMSE regardless of reference, normalization, or regression"))

# Plt 4a but divided by tissue so there is a baseline line for each tissue
# separately.
plt4c <- Plot_FacetBoxPlot(subset(errs_box, error_type == "RMSE"),
                           fill = "tissue",
                           x_axis = "reference_data_name",
                           fill_colors = tissue_colors,
                           facet_vars = c("tissue", "algorithm")) +
  geom_hline(aes(yintercept = min_val),
             data = subset(baselines_plot_tissue, error_type == "RMSE"),
             linetype = "twodash")
print(plt4c + plot_annotation(title = "Spread of RMSE by tissue, regardless of normalization, or regression"))


# mAPE
plt5a <- Plot_FacetBoxPlot(subset(errs_box, error_type == "MAPE"),
                           fill = "tissue",
                           x_axis = "reference_data_name",
                           fill_colors = tissue_colors,
                           facet_vars = c("test_data_name", "algorithm")) +
  geom_hline(aes(yintercept = min_val),
             data = subset(baselines_plot, error_type == "MAPE"),
             linetype = "twodash")

print(plt5a + plot_annotation(title = "Spread of MAPE regardless of normalization, or regression"))

plt5b <- Plot_FacetBoxPlot(subset(errs_box_noref, error_type == "MAPE"),
                           fill = "tissue",
                           x_axis = "algorithm",
                           fill_colors = tissue_colors,
                           facet_vars = "test_data_name") +
  geom_hline(aes(yintercept = min_val),
             data = subset(baselines_plot, error_type == "MAPE"),
             linetype = "twodash")
print(plt5b + plot_annotation(title = "Spread of MAPE regardless of reference, normalization, or regression"))

# Plt 5a but divided by tissue so there is a baseline line for each tissue
# separately.
plt5c <- Plot_FacetBoxPlot(subset(errs_box, error_type == "MAPE"),
                           fill = "tissue",
                           x_axis = "reference_data_name",
                           fill_colors = tissue_colors,
                           facet_vars = c("tissue", "algorithm")) +
  geom_hline(aes(yintercept = min_val),
             data = subset(baselines_plot_tissue, error_type == "MAPE"),
             linetype = "twodash")
print(plt5c + plot_annotation(title = "Spread of MAPE by tissue, regardless of normalization, or regression"))


# Normalization vs regression --------------------------------------------------

errs_box_norm <- errs_melt %>%
  Create_BoxStats(c("tissue", "test_data_name", "algorithm", "error_type",
                    "normalization", "regression_method"))

baselines_plot_norm <- baselines_random %>%
  Create_BoxStats(c("tissue", "normalization", "error_type"))

# TODO these plots are kind of large, need to figure out how to collapse into
# meaningful information
plt6 <- Plot_FacetBoxPlot(subset(errs_box_norm, error_type == "Correlation"),
                          x_axis = "normalization",
                          fill = "regression_method",
                          fill_colors = regression_colors,
                          facet_vars = c("tissue", "algorithm")) +
  geom_hline(aes(yintercept = max_val, color = normalization),
             data = subset(baselines_plot_norm, error_type == "Correlation"),
             linetype = "twodash")
print(plt6 + plot_annotation(title = "Correlation: normalization vs regression"))

plt7 <- Plot_FacetBoxPlot(subset(errs_box_norm, error_type == "RMSE"),
                          x_axis = "normalization",
                          fill = "regression_method",
                          fill_colors = regression_colors,
                          facet_vars = c("tissue", "algorithm")) +
  geom_hline(aes(yintercept = min_val, color = normalization),
             data = subset(baselines_plot_norm, error_type == "RMSE"),
             linetype = "twodash")
print(plt7 + plot_annotation(title = "RMSE: normalization vs regression"))

plt8 <- Plot_FacetBoxPlot(subset(errs_box_norm, error_type == "MAPE"),
                          x_axis = "normalization",
                          fill = "regression_method",
                          fill_colors = regression_colors,
                          facet_vars = c("tissue", "algorithm")) +
  geom_hline(aes(yintercept = min_val, color = normalization),
             data = subset(baselines_plot_norm, error_type == "MAPE"),
             linetype = "twodash")
print(plt8 + plot_annotation(title = "MAPE: normalization vs regression"))

# Algorithms collapsed
errs_box_noalg <- errs_melt %>%
  Create_BoxStats(c("tissue", "error_type",
                    "normalization", "regression_method"))

plt12 <- Plot_FacetBoxPlot(subset(errs_box_noalg, error_type == "Correlation"),
                           x_axis = "normalization",
                           fill = "regression_method",
                           fill_colors = regression_colors,
                           facet_vars = c("tissue")) +
  geom_hline(aes(yintercept = max_val, color = normalization),
             data = subset(baselines_plot_norm, error_type == "Correlation"),
             linetype = "twodash")
print(plt12 + plot_annotation(title = "Correlation: normalization vs regression"))

plt13 <- Plot_FacetBoxPlot(subset(errs_box_noalg, error_type == "RMSE"),
                           x_axis = "normalization",
                           fill = "regression_method",
                           fill_colors = regression_colors,
                           facet_vars = c("tissue")) +
  geom_hline(aes(yintercept = min_val, color = normalization),
             data = subset(baselines_plot_norm, error_type == "RMSE"),
             linetype = "twodash")
print(plt13 + plot_annotation(title = "RMSE: normalization vs regression"))

plt14 <- Plot_FacetBoxPlot(subset(errs_box_noalg, error_type == "MAPE"),
                           x_axis = "normalization",
                           fill = "regression_method",
                           fill_colors = regression_colors,
                           facet_vars = c("tissue")) +
  geom_hline(aes(yintercept = min_val, color = normalization),
             data = subset(baselines_plot_norm, error_type == "MAPE"),
             linetype = "twodash")
print(plt14 + plot_annotation(title = "MAPE: normalization vs regression"))

for (err_metric in c("Correlation", "RMSE", "MAPE")) {
  errs_sub <- subset(errs_melt, error_type == err_metric) %>%
    Create_BoxStats(c("test_data_name", "normalization", "regression_method"))

  plt <- Plot_FacetBoxPlot(errs_sub,
                           x_axis = "normalization",
                           facet_vars = "test_data_name",
                           fill = "regression_method",
                           fill_colors = regression_colors)
  print(plt + plot_annotation(title = str_glue("{err_metric}: normalization vs regression")))

  # One plot for each bulk data set
  for (bulk_name in unique(errs_sub$test_data_name)) {

    errs_sub2 <- subset(errs_box_norm,
                        test_data_name == bulk_name & error_type == err_metric)

    plt <- Plot_FacetBoxPlot(errs_sub2,
                             x_axis = "normalization",
                             fill = "regression_method",
                             facet_vars = c("tissue", "algorithm"),
                             fill_colors = regression_colors) +
      ggtitle(paste(bulk_name, err_metric))

    print(plt)
  }
}


# Baseline only ----------------------------------------------------------------

baselines_plot <- Create_BoxStats(baselines_random,
                                  grouping_cols = c("tissue", "normalization",
                                                    "regression_method",
                                                    "error_type"))

plt12 <- Plot_FacetBoxPlot(subset(baselines_plot, error_type == "Correlation"),
                           x_axis = "normalization",
                           fill = NULL, #"regression_method",
                           fill_colors = NULL,#regression_colors,
                           facet_vars = "tissue",
                           color = "regression_method") +
  scale_color_manual(values = regression_colors)
print(plt12 + plot_annotation(title = "Baseline correlation: normalization vs regression"))

plt13 <- Plot_FacetBoxPlot(subset(baselines_plot, error_type == "RMSE"),
                           x_axis = "normalization",
                           fill = NULL, #"regression_method",
                           fill_colors = NULL, #regression_colors,
                           facet_vars = "tissue",
                           color = "regression_method") +
  scale_color_manual(values = regression_colors)
print(plt13 + plot_annotation(title = "Baseline RMSE: normalization vs regression"))

plt14 <- Plot_FacetBoxPlot(subset(baselines_plot, error_type == "MAPE"),
                           x_axis = "normalization",
                           fill = NULL, #"regression_method",
                           fill_colors = NULL, #regression_colors,
                           facet_vars = "tissue",
                           color = "regression_method") +
  scale_color_manual(values = regression_colors)
print(plt14 + plot_annotation(title = "Baseline MAPE: normalization vs regression"))


# Errors better than baseline only ---------------------------------------------

errs_better <- merge(errs_melt, baselines_plot,
                     by = c("tissue", "normalization", "regression_method", "error_type")) %>%
  group_by(error_type) %>%
  mutate(better = if(unique(error_type) == "Correlation") (value >= max_val) else (value <= min_val)) #%>%
#subset(better == TRUE) %>%
#select(-better, -median_val, -upper_quartile, -lower_quartile)

better_stats <- errs_better %>%
  group_by(tissue, algorithm) %>%
  dplyr::summarize(count = n(),
                   pct_better_than_baseline = sum(better) / n(),
                   .groups = "drop")

plt <- ggplot(better_stats, aes(x = algorithm, y = pct_better_than_baseline, fill = algorithm)) +
  geom_col() +
  facet_wrap(~tissue) +
  scale_fill_manual(values = algorithm_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(plt + plot_annotation(title = "Percent of parameters better than baseline"))

# Quality stats ----------------------------------------------------------------

# TODO these fields don't exist in best_errs_plot, they need to come from the
# quality stats script
quality_stats <- subset(best_errs_plot, algorithm != "Baseline") %>%
  select(-Correlation, -RMSE, -MAPE) %>%
  distinct()

## Bad inhibitory ratio, norm vs regression ------------------------------------

qbox_stats <- quality_stats %>%
  group_by(test_data_name, tissue, algorithm) %>%
  dplyr::summarize(max_val = max(pct_bad_inhibitory_ratio),
                   min_val = min(pct_bad_inhibitory_ratio),
                   median_val = median(pct_bad_inhibitory_ratio),
                   upper_quartile = quantile(pct_bad_inhibitory_ratio, probs = 0.75, na.rm = TRUE),
                   lower_quartile = quantile(pct_bad_inhibitory_ratio, probs = 0.25, na.rm = TRUE),
                   .groups = "drop")

qbox_stats_norm <- quality_stats %>%
  group_by(test_data_name, tissue, algorithm, normalization, regression_method) %>%
  dplyr::summarize(max_val = max(pct_bad_inhibitory_ratio),
                   min_val = min(pct_bad_inhibitory_ratio),
                   median_val = median(pct_bad_inhibitory_ratio),
                   upper_quartile = quantile(pct_bad_inhibitory_ratio, probs = 0.75, na.rm = TRUE),
                   lower_quartile = quantile(pct_bad_inhibitory_ratio, probs = 0.25, na.rm = TRUE),
                   .groups = "drop")

plt <- ggplot(quality_stats, aes(x = normalization, y = pct_bad_inhibitory_ratio,
                                 fill = regression_method)) +
  geom_boxplot() + theme_bw() + #facet_wrap(~tissue)
  facet_grid(rows = vars(tissue), cols = vars(algorithm)) +
  scale_fill_manual(values = regression_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(plt + plot_annotation(title = "Percent bad inhibitory ratio"))

plt15 <- Plot_FacetBoxPlot(qbox_stats_norm,
                           x_axis = "normalization",
                           fill = "regression_method",
                           fill_colors = regression_colors,
                           facet_vars = c("tissue", "algorithm"))
print(plt15 + plot_annotation(title = "Percent bad inhibitory ratio"))

plt16 <- Plot_FacetBoxPlot(qbox_stats,
                           x_axis = "algorithm",
                           fill = "algorithm",
                           fill_colors = algorithm_colors,
                           facet_vars = "tissue")
print(plt16 + plot_annotation(title = "Percent bad inhibitory ratio"))


# params passing QC, norm vs regression

# Fill in stats for reference data sets that didn't have any valid param sets
# and therefore aren't in this data frame
file_params <- expand.grid(tissue = unique(quality_stats$tissue),
                           reference_data_name = datasets,
                           normalization = unique(quality_stats$normalization),
                           regression_method = unique(quality_stats$regression_method),
                           algorithm = unique(quality_stats$algorithm))

tissue_data <- quality_stats %>%
  select(tissue, test_data_name) %>%
  distinct()

file_params <- merge(file_params, tissue_data, by = "tissue")

# CibersortX does not have tmm normalization and is the only algorithm with
# two possible reference_input_types
cibersort_only <- file_params %>%
  subset(algorithm == "CibersortX" & normalization %in% c("counts/cpm", "tpm")) %>%
  merge(data.frame(algorithm = "CibersortX",
                   reference_input_type = c("cibersortx", "signature")),
        by = "algorithm")

# Music and Scaden do not have tmm or tpm normalization
file_params <- file_params %>%
  subset(algorithm != "CibersortX") %>%
  subset(!(algorithm %in% c("Music", "Scaden")) |
           normalization == "counts/cpm")

qstats_all <- quality_stats %>%
  select(-param_id, -pct_bad_inhibitory_ratio) %>%
  distinct() %>%
  merge(file_params, by = colnames(file_params), all = TRUE) %>%
  merge(cibersort_only, by = colnames(cibersort_only), all = TRUE)

# These values are NA where there was missing data, set to 0
qstats_all$pct_valid_results[is.na(qstats_all$pct_valid_results)] <- 0

plt <- ggplot(qstats_all, aes(x = normalization, y = pct_valid_results,
                              fill = regression_method)) +
  geom_boxplot() + theme_bw() + #facet_wrap(~tissue)
  facet_grid(rows = vars(tissue), cols = vars(algorithm)) +
  scale_fill_manual(values = regression_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(plt + plot_annotation(title = "Percent valid output"))

qbox_stats2 <- qstats_all %>% group_by(test_data_name, algorithm) %>%
  dplyr::summarize(max_val = max(pct_valid_results),
                   min_val = min(pct_valid_results),
                   median_val = median(pct_valid_results),
                   upper_quartile = quantile(pct_valid_results, probs = 0.75, na.rm = TRUE),
                   lower_quartile = quantile(pct_valid_results, probs = 0.25, na.rm = TRUE),
                   .groups = "drop")

plt <- Plot_FacetBoxPlot(qbox_stats2, x_axis = "algorithm", fill = "algorithm",
                         fill_colors = algorithm_colors,
                         facet_vars = "test_data_name")
print(plt + plot_annotation(title = "Percent valid output"))

# Top 3 param sets for each tissue ---------------------------------------------

ranked <- readRDS(file.path(dir_analysis, str_glue("ranked_errors_{granularity}.rds")))

ranked_df <- ranked$ranked_errors_all %>%
  mutate(normalization = str_replace(normalization, "counts", "cpm"),
         normalization = str_replace(normalization, "cpm", "counts/cpm"),
         normalization = str_replace(normalization, "log_", ""))

CountDf <- function(groups) {
  ranked_df %>%
    group_by_at(groups) %>%
    dplyr::summarize(count = n(), .groups = "drop")
}

plt <- ggplot(CountDf(c("tissue", "signature")),
              aes(x = tissue, y = signature, color = count, size = count)) +
  geom_count() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5)
print(plt + plot_annotation(title = "Best signature per tissue (unfiltered)"))

plt <- ggplot(CountDf(c("tissue", "reference_data_name")),
              aes(x = tissue, y = reference_data_name, color = count, size = count)) +
  geom_count() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5)
print(plt + plot_annotation(title = "Best reference per tissue (unfiltered)"))

plt <- ggplot(CountDf(c("tissue", "data_transform")),
              aes(x = tissue, y = data_transform, color = count, size = count)) +
  geom_count() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5)
print(plt + plot_annotation(title = "Best data transform per tissue (unfiltered)"))

plt <- ggplot(CountDf(c("tissue", "algorithm")),
              aes(x = tissue, y = algorithm, color = count, size = count)) +
  geom_count() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5)
print(plt + plot_annotation(title = "Best algorithm per tissue (unfiltered)"))

plt <- ggplot(CountDf(c("type", "algorithm")),
              aes(x = type, y = algorithm, color = count, size = count)) +
  geom_count() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5)
print(plt + plot_annotation(title = "Best algorithm per error metric (unfiltered)"))

plt <- ggplot(CountDf(c("type", "data_transform")),
              aes(x = type, y = data_transform, color = count, size = count)) +
  geom_count() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5)
print(plt + plot_annotation(title = "Best data transform per error metric (unfiltered)"))

plt <- ggplot(CountDf(c("normalization", "regression_method")),
              aes(x = normalization, y = regression_method, color = count, size = count)) +
  geom_count() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5)
print(plt + plot_annotation(title = "Best normalization vs regression method (unfiltered)"))

plt <- ggplot(CountDf(c("tissue", "regression_method")),
              aes(x = tissue, y = regression_method, color = count, size = count)) +
  geom_count() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5)
print(plt + plot_annotation(title = "Best regression method per tissue (unfiltered)"))

plt <- ggplot(CountDf(c("tissue", "normalization")),
              aes(x = tissue, y = normalization, color = count, size = count)) +
  geom_count() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5)
print(plt + plot_annotation(title = "Best normalization per tissue (unfiltered)"))

plt <- ggplot(CountDf(c("tissue", "normalization", "regression_method")),
              aes(x = normalization, y = regression_method, color = count, size = count)) +
  geom_count() + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5) +
  facet_wrap(~tissue, nrow = 1)
print(plt + plot_annotation(title = "Best normalization/regression per tissue (unfiltered)"))


# Same thing but limited to the "best" signatures used for analysis

ranked_df <- ranked$ranked_errors_best_signatures %>%
  mutate(normalization = str_replace(normalization, "counts", "cpm"),
         normalization = str_replace(normalization, "cpm", "counts/cpm"),
         normalization = str_replace(normalization, "log_", ""))

CountDf <- function(groups) {
  ranked_df %>%
    group_by_at(groups) %>%
    dplyr::summarize(count = n(), .groups = "drop")
}

plt <- ggplot(CountDf(c("tissue", "signature")),
              aes(x = tissue, y = signature, color = count, size = count)) +
  geom_count() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5)
print(plt + plot_annotation(title = "Best signature per tissue"))

# This should be the same as best signature since signature == reference_data_name
plt <- ggplot(CountDf(c("tissue", "reference_data_name")),
              aes(x = tissue, y = reference_data_name, color = count, size = count)) +
  geom_count() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5)
print(plt + plot_annotation(title = "Best reference per tissue"))

plt <- ggplot(CountDf(c("tissue", "data_transform")),
              aes(x = tissue, y = data_transform, color = count, size = count)) +
  geom_count() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5)
print(plt + plot_annotation(title = "Best data transform per tissue"))

plt <- ggplot(CountDf(c("tissue", "algorithm")),
              aes(x = tissue, y = algorithm, color = count, size = count)) +
  geom_count() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5)
print(plt + plot_annotation(title = "Best algorithm per tissue"))

plt <- ggplot(CountDf(c("type", "algorithm")),
              aes(x = type, y = algorithm, color = count, size = count)) +
  geom_count() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5)
print(plt + plot_annotation(title = "Best algorithm per error metric"))

plt <- ggplot(CountDf(c("type", "data_transform")),
              aes(x = type, y = data_transform, color = count, size = count)) +
  geom_count() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5)
print(plt + plot_annotation(title = "Best data transform per error metric"))

plt <- ggplot(CountDf(c("normalization", "regression_method")),
              aes(x = normalization, y = regression_method, color = count, size = count)) +
  geom_count() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5)
print(plt + plot_annotation(title = "Best normalization vs regression"))

plt <- ggplot(CountDf(c("tissue", "regression_method")),
              aes(x = tissue, y = regression_method, color = count, size = count)) +
  geom_count() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5)
print(plt + plot_annotation(title = "Best regression method per tissue"))

plt <- ggplot(CountDf(c("tissue", "normalization")),
              aes(x = tissue, y = normalization, color = count, size = count)) +
  geom_count() + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5)
print(plt + plot_annotation(title = "Best normalization per tissue"))

plt <- ggplot(CountDf(c("tissue", "normalization", "regression_method")),
              aes(x = normalization, y = regression_method, color = count, size = count)) +
  geom_count() + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5) +
  facet_wrap(~tissue, nrow = 1)
print(plt + plot_annotation(title = "Best normalization/regression per tissue"))


dev.off()

# Detailed plots per reference data set ----------------------------------------

file_params <- best_errors %>%
  select(reference_data_name, test_data_name, algorithm,
         normalization, regression_method) %>%
  distinct()

best_params <- best_errs_plot %>%
  select(-cor, -rMSE, -mAPE, -pct_valid_results, -pct_bad_inhibitory_ratio) %>%
  mutate(tissue = str_replace(tissue, paste0(test_data_name, " "), ""))

significance_list <- readRDS(file.path(dir_analysis,
                                       str_glue("significance_lists_{granularity}.rds")))
significance <- do.call(rbind, significance_list$significance_props_single) %>%
  subset(algorithm != "Baseline") %>%
  mutate(param_id = str_replace(avg_id, paste0(as.character(tissue), "_"), "")) %>%
  subset(param_id %in% best_params$param_id)

for (bulk_dataset in bulk_datasets) {
  bulk_se <- Load_BulkData(bulk_dataset)
  metadata <- as.data.frame(colData(bulk_se))

  best_params_sub <- subset(best_params, test_data_name == bulk_dataset)
  best_ests <- Get_AllBestEstimatesAsDf(bulk_dataset, granularity, metadata,
                                        best_params = best_params_sub)

  for (dataset in datasets) {
    pdf(file.path(dir_figures,
                  str_glue("error_plots_{bulk_dataset}_{dataset}_{granularity}_detailed.pdf")),
        width=10, height = 12)

    best_params_dataset <- subset(best_params_sub, reference_data_name == dataset)
    best_ests_dataset <- subset(best_ests,
                                param_id %in% best_params_dataset$param_id) %>%
      merge(best_params_dataset, by = c("tissue", "param_id", "algorithm"),
            all.x = FALSE)

    best_ests_dataset <- best_ests_dataset %>%
      group_by(tissue, algorithm, normalization, regression_method) %>%
      mutate(title = paste(algorithm, "params",
                           as.numeric(factor(param_id)))) %>%
      ungroup()

    best_errs_dataset <- best_errs_plot %>%
      select(-pct_valid_results, -pct_bad_inhibitory_ratio) %>%
      subset(reference_data_name == dataset & test_data_name == bulk_dataset) %>%
      mutate(tissue = str_replace(tissue, paste0(test_data_name, " "), "")) %>%
      group_by(tissue, algorithm, normalization, regression_method) %>%
      mutate(title = paste(algorithm, "params",
                           as.numeric(factor(param_id)))) %>%
      ungroup()

    ## Compare error vs normalization and regression ---------------------------

    errs_box <- melt(best_errs_dataset, variable.name = "error_type") %>%
      Create_BoxStats(c("tissue", "algorithm", "normalization",
                        "regression_method", "error_type"))

    best_vals <- errs_box %>%
      group_by(error_type, tissue, normalization, algorithm) %>%
      summarize(best_max = max(max_val),
                best_min = min(min_val),
                .groups = "keep") %>%
      mutate(best_val = if (error_type == "cor") best_max else best_min)

    for (err_type in c("cor", "rMSE", "mAPE")) {
      plots <- lapply(unique(errs_box$tissue), function(tiss) {
        plt <- Plot_FacetBoxPlot(subset(errs_box, error_type == err_type & tissue == tiss),
                                 x_axis = "normalization",
                                 fill = "regression_method",
                                 facet_vars = c("tissue", "algorithm"),
                                 fill_colors = regression_colors) +
          geom_hline(aes(yintercept = best_val, color = normalization),
                     data = subset(best_vals, error_type == err_type & tissue == tiss),
                     linetype = "twodash")
        return(plt)
      })

      plt <- plots[[1]]
      for (N in 2:length(plots)) {
        plt <- plt / plots[[N]]
      }
      plot_id <- str_glue("{bulk_dataset} x {dataset} Errors (Normalization vs Regression)")
      subtitle <- str_glue("Best {err_type}. Dashed lines are the best score for that normalization.")
      print(plt + plot_annotation(title = plot_id, subtitle = subtitle))
    }


    ## Separate out by normalization / regression ------------------------------

    params_plot <- best_params_dataset %>%
      select(normalization, regression_method) %>%
      distinct()

    for (R in 1:nrow(params_plot)) {
      norm <- params_plot$normalization[R]
      regression <- params_plot$regression_method[R]

      plot_id <- str_glue(
        paste0("{bulk_dataset} x {dataset} (",
               "Normalization: {norm}, Regression: {regression})")
      )

      best_ests_params <- subset(best_ests_dataset,
                                 normalization == norm &
                                   regression_method == regression)

      best_errs_params <- subset(best_errs_dataset,
                                 normalization == norm &
                                   regression_method == regression) %>%
        melt(variable.name = "error_type")


      ### Error comparisons ----------------------------------------------------

      errs_box <- Create_BoxStats(best_errs_params,
                                  c("tissue", "param_id", "title", "algorithm",
                                    "error_type"))

      best_vals <- errs_box %>%
        group_by(tissue, error_type) %>%
        summarize(best_max = max(max_val),
                  best_min = min(min_val),
                  .groups = "drop") %>%
        group_by(tissue, error_type) %>%
        mutate(best_val = if (error_type == "cor") best_max else best_min)

      for (err_type in c("cor", "rMSE", "mAPE")) {
        plots <- lapply(unique(best_errs_params$tissue), function(tiss) {
          plt <- ggplot(subset(best_errs_params, error_type == err_type & tissue == tiss),
                        aes(x = algorithm, y = value, fill = algorithm)) +
            geom_dotplot(binaxis = "y", stackdir = "center") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            facet_wrap(~tissue) +
            geom_hline(aes(yintercept = best_val),
                       data = subset(best_vals, error_type == err_type & tissue == tiss),
                       linetype = "twodash")
          return(plt)
        })

        plt <- plots[[1]]
        for (N in 2:length(plots)) {
          plt <- plt / plots[[N]]
        }
        subtitle <- str_glue("Best {err_type} (each dot is a single parameter set)")
        print(plt + plot_annotation(title = plot_id, subtitle = subtitle))
      }

      ### Estimate comparisons -------------------------------------------------

      # TODO this is percent cells, not percent RNA
      merscope <- data.frame(celltype = c("Astrocyte", "Excitatory", "Inhibitory",
                                          "Microglia", "Oligodendrocyte", "OPC",
                                          "Vascular"),
                             pct = c(10, 34, 7, 6, 32.5, 3.5, 11.5)/100)
      #merscope = subset(merscope, celltype %in% c("Astrocyte", "Oligodendrocyte"))
      #merscope = rbind(merscope, data.frame(celltype = c("Endothelial", "Pericyte"),
      #                                      pct = c(8, 3.5)/100))

      # These graphs need to be broken down by tissue so there is enough room
      # for all the cell types on the page
      for (tiss in unique(best_ests_params$tissue)) {
        best_ests_tissue <- subset(best_ests_params, tissue == tiss)
        n_cols <- if (granularity == "broad_class") 4 else 6
        plot_id_tissue <- str_glue(
          paste0("{bulk_dataset} {tiss} x {dataset} (",
                 "Normalization: {norm}, Regression: {regression})")
        )

        plt <- ggplot(best_ests_tissue, aes(x = algorithm, y = percent, color = title)) +
          geom_boxplot(width = 0.5) + theme_bw() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          theme(legend.position = "bottom") +
          facet_wrap(~celltype, scales = "fixed", ncol = n_cols) +
          #facet_grid(celltype ~ tissue, scales = "fixed") +
          geom_hline(data = merscope, mapping = aes(yintercept = pct), color = "red")

        subtitle = "Estimated proportions across all subjects (each bar is a single parameter set)"
        print(plt + plot_annotation(title = plot_id_tissue, subtitle = subtitle))
      }

      # Significance graphs per algorithm
      for (alg in unique(best_ests_params$algorithm)) {
        plots <- lapply(unique(best_ests_params$tissue), function(tiss) {
          ests_alg <- subset(best_ests_params, algorithm == alg &
                               tissue == tiss &
                               diagnosis %in% c("CT", "AD"))
          significance_sub <- subset(significance,
                                     param_id %in% ests_alg$param_id) %>%
            select(tissue, celltype, param_id, significant_01)

          ests_box <- merge(ests_alg, significance_sub,
                            by = c("tissue", "celltype", "param_id")) %>%
            mutate(value = percent) %>%
            Create_BoxStats(c("tissue", "celltype", "param_id", "title",
                              "algorithm", "diagnosis", "significant_01"))

          # TODO fill colors for diagnosis
          plt <- Plot_FacetBoxPlot(ests_box,
                                   x_axis = "celltype",
                                   facet_vars = c("tissue", "title"),
                                   fill_colors = NULL,
                                   fill = "diagnosis",
                                   color = "significant_01",
                                   width = 0.5) +
            scale_color_manual(values = c("#dddddd", "#000000"))
          return(plt)
        })

        plt <- plots[[1]]
        for (N in 2:length(plots)) {
          plt <- plt / plots[[N]]
        }
        print(plt + plot_annotation(title = plot_id,
                                    subtitle = "Estimates, AD vs Control"))
      }
    }

    dev.off()
  }
}

significance <- readRDS(file.path(dir_analysis,
                                  str_glue("significance_lists_{granularity}.rds")))

sig_toplevel <- do.call(rbind, significance$significance_props_toplevel)
#sig_avg <- do.call(rbind, significance$significance_props_all)

pdf(file.path(dir_figures,
              str_glue("significance_plots_{granularity}_summary.pdf")),
    width=10, height = 12)

sig_toplevel$normalization <- sapply(sig_toplevel$avg_id, function(A) {
  for (norm in c("counts|cpm", "tmm", "tpm")) {
    if (grepl(norm, A)) {
      if (norm == "counts|cpm") { return("CPM") }
      else { return(toupper(norm)) }
    }
  }
})

sig_toplevel$regression_method <- sapply(sig_toplevel$avg_id, function(A) {
  for (reg in c("none", "edger", "lme", "dream")) {
    if (grepl(reg, A)) {
      if (reg == "edger") { return("edgeR") }
      if (reg == "lme") { return("LME") }
      if (reg == "none") { return("No regression") }
      return(reg)
    }
  }
})

# Fill in any missing values for combinations of parameters
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

best_dt <- CountDf(c("tissue", "data_transform")) %>%
  group_by(tissue) %>%
  slice_max(order_by = count, with_ties = FALSE) %>%
  mutate(normalization = str_replace(data_transform, " \\+.*", ""),
         regression_method = str_replace(data_transform, ".*\\+ ", "")) %>%
  as.data.frame()

best_dt <- tidyr::expand_grid(best_dt, algorithm = unique(sig_toplevel$algorithm))
best_dt$normalization[best_dt$algorithm == "Music"] <- "counts/cpm"
best_dt$normalization[best_dt$normalization == "tmm" & best_dt$algorithm == "CibersortX"] <- "counts/cpm"

if (any(best_dt$normalization == "tpm")) {
  print("WARNING: at least one tissue has 'tpm' as the best normalization")
}

best_dt <- best_dt %>%
  mutate(data_transform = paste(normalization, "+", regression_method)) %>%
  select(tissue, algorithm, data_transform) %>%
  distinct()

sig_final <- merge(sig_toplevel, best_dt,
                   by = c("tissue", "algorithm", "data_transform"),
                   all.x = FALSE)

# Cap log2_fc values to +/- 1 so color scaling is better. Need to account for
# Inf and -Inf values. Set non-significant log2 values to NA.
sig_final$log2_fc[sig_final$fc == 0] <- NA
sig_final$log2_fc[is.infinite(sig_final$log2_fc)] <- 1
sig_final$log2_fc[sig_final$log2_fc > 1] <- 1
sig_final$log2_fc[sig_final$log2_fc < -1] <- -1

sig_final$log2_fc[sig_final$p_adj_thresh >= 0.05] <- NA

plot_limit <- c(-1, 1)
n_cols <- if (granularity == "broad_class") 4 else 6

plt <- ggplot(sig_final, aes(x = tissue, y = algorithm, fill = log2_fc)) +
  geom_tile(color = "black") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~celltype, ncol = n_cols) +
  #facet_grid(rows = vars(algorithm), cols = vars(celltype)) +
  coord_fixed() +
  scale_fill_distiller(palette = "RdBu", limit = plot_limit, na.value = "white") +
  ggtitle("Significant changes at p < 0.05 (using best data for each tissue)")
print(plt)

# Limit to cell types that have something significant
ok <- subset(sig_final, is.finite(log2_fc))
sig_final_ok <- subset(sig_final, celltype %in% unique(ok$celltype))
plt <- ggplot(sig_final_ok, aes(x = tissue, y = algorithm, fill = log2_fc)) + #, size = -log_p)) +
  geom_tile(color = "black") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~celltype, ncol = n_cols) +
  #facet_grid(rows = vars(algorithm), cols = vars(celltype)) +
  coord_fixed() +
  scale_fill_distiller(palette = "RdBu", limit = plot_limit, na.value = "white") +
  ggtitle("Significant changes at p < 0.05 (using best data for each tissue)")
print(plt)

# At p <= 0.01

sig_final_01 <- sig_final
sig_final_01$log2_fc[sig_final_01$p_adj_thresh >= 0.01] <- NA

plt <- ggplot(sig_final_01, aes(x = tissue, y = algorithm, fill = log2_fc)) +
  geom_tile(color = "black") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~celltype, ncol = n_cols) +
  #facet_grid(rows = vars(algorithm), cols = vars(celltype)) +
  coord_fixed() +
  scale_fill_distiller(palette = "RdBu", limit = plot_limit, na.value = "white") +
  ggtitle("Significant changes at p < 0.01 (using best data for each tissue)")
print(plt)

# Limit to cell types that have something significant
ok <- subset(sig_final_01, is.finite(log2_fc))
sig_final_ok <- subset(sig_final_01, celltype %in% unique(ok$celltype))
plt <- ggplot(sig_final_ok, aes(x = tissue, y = algorithm, fill = log2_fc)) + #, size = -log_p)) +
  geom_tile(color = "black") + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~celltype, ncol = n_cols) +
  #facet_grid(rows = vars(algorithm), cols = vars(celltype)) +
  coord_fixed() +
  scale_fill_distiller(palette = "RdBu", limit = plot_limit, na.value = "white") +
  ggtitle("Significant changes at p < 0.01 (using best data for each tissue)")
print(plt)



for (dt in unique(sig_toplevel$data_transform)) {
  sig_filt <- subset(sig_toplevel, data_transform == dt)

  if (nrow(sig_filt) == 0) {
    next
  }
  print(dt)

  # Cap log2_fc values to +/- 2 so color scaling is better. Need to account for
  # Inf and -Inf values. Set non-significant log2 values to NA.
  sig_filt$log2_fc[sig_filt$fc == 0] <- NA
  sig_filt$log2_fc[is.infinite(sig_filt$log2_fc)] <- 1
  sig_filt$log2_fc[sig_filt$log2_fc > 1] <- 1
  sig_filt$log2_fc[sig_filt$log2_fc < -1] <- -1

  sig_filt$log2_fc[sig_filt$p_adj_thresh >= 0.05] <- NA

  # Cap log_p values to -5 so color scaling is better
  #sig_filt$log_p[sig_filt$log_p < -5] <- -5

  #y_axis <- sig_filt %>%
  #  select(avg_id, tissue, algorithm) %>%
  #  distinct() %>%
  #  group_by(tissue, algorithm) %>%
  #  mutate(name = paste("Trial", 1:n(), sep = " "))

  #sig_filt <- merge(sig_filt, y_axis, by = c("avg_id", "tissue", "algorithm"))

  plot_limit <- c(-1, 1)
  n_cols <- if (granularity == "broad_class") 4 else 6

  plt <- ggplot(sig_filt, aes(x = tissue, y = algorithm, fill = log2_fc)) +
    geom_tile(color = "black") + theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~celltype, ncol = n_cols) +
    #facet_grid(rows = vars(algorithm), cols = vars(celltype)) +
    coord_fixed() +
    scale_fill_distiller(palette = "RdBu", limit = plot_limit, na.value = "white") +
    ggtitle(dt)
  print(plt)

  # Limit to cell types that have something significant
  ok <- subset(sig_filt, is.finite(log2_fc))
  sig_filt_ok <- subset(sig_filt, celltype %in% unique(ok$celltype))
  plt <- ggplot(sig_filt_ok, aes(x = tissue, y = algorithm, fill = log2_fc)) + #, size = -log_p)) +
    geom_tile(color = "black") + theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~celltype, ncol = n_cols) +
    #facet_grid(rows = vars(algorithm), cols = vars(celltype)) +
    coord_fixed() +
    scale_fill_distiller(palette = "RdBu", limit = plot_limit, na.value = "white") +
    ggtitle(dt)
  print(plt)
}

dev.off()
