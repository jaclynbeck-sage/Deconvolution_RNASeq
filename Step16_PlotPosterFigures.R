library(Matrix)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(dplyr)
library(stringr)
library(patchwork)
library(sagethemes)

source(file.path("functions", "FileIO_HelperFunctions.R"))
source(file.path("functions", "Step15_Analysis_HelperFunctions.R"))
source(file.path("functions", "Step15_Plotting_HelperFunctions.R"))

options(scipen = 999)

bulk_datasets <- c("Mayo", "MSBB", "ROSMAP")

best_errors_list <- Load_BestErrors()

quality_stats <- Load_QualityStats()

# Unpack error variables into environment for readability
list2env(best_errors_list, globalenv())

tissues_use <- c("TCX", "PHG", "ACC") # When subsetting
tissues_use_full <- c("Mayo TCX", "MSBB PHG", "ROSMAP ACC") # When subsetting


# Color setup ------------------------------------------------------------------

algs <- unique(best_errors$algorithm)
regs <- unique(best_errors$regression_method)

# Algorithms
algorithm_colors <- RColorBrewer::brewer.pal(length(algs), "Set2")
names(algorithm_colors) <- sort(algs)

# Regression methods
regression_colors <- c(
  sage_colors$royal[["600"]],
  sage_colors$turquoise[["500"]],
  sage_colors$butterscotch[["500"]],
  sage_colors$apple[["500"]]
) |> paste0("CC") # Add 80% opacity

granularity_colors <- c(
  sage_colors$slate[["600"]],
  sage_colors$blueberry[["400"]]
)

# Tissues (modified viridis turbo color scheme)
tiss <- colSums(table(best_errors$tissue, best_errors$test_data_name) > 1)
tissue_colors <- c(viridis::turbo(tiss[["Mayo"]], begin = 0.1, end = 0.2),
                   viridis::turbo(tiss[["MSBB"]], begin = 0.35, end = 0.6),
                   viridis::turbo(tiss[["ROSMAP"]], begin = 0.7, end = 0.9))
names(tissue_colors) <- sort(unique(best_errors$tissue_full))

# MSBB colors need to be darker and a little more differentiated -- original
# colors are 20EAABFF, 67FD68FF, AEFA37FF, E1DD37FF
tissue_colors[grepl("MSBB", names(tissue_colors))] <- c("#00CA8BFF", "#37CD38FF",
                                                        "#7ECA07FF", "#C1BD17FF")

# Slightly more pastel than default looks better
tissue_fill_colors <- str_replace(tissue_colors, "FF$", "88")
names(tissue_fill_colors) <- names(tissue_colors)


# ggplot2 settings that generally apply to all plots
poster_theme <- theme(
  # All text
  text = element_text(family = "DM Sans"),
  # Add space between axis titles and axis text
  axis.title.x = element_text(margin = margin(t = 10, unit = "pt")),
  axis.title.y = element_text(margin = margin(r = 10, unit = "pt")),
  # Plot title
  plot.title = element_text(face = "bold", size = 10),
  # Facet labels
  strip.background = element_blank(),
  strip.text = element_text(size = 10, face = "bold"),
  strip.placement.y = "outside"
)

numerical_legend <- theme(
  legend.title = element_text(size = 10, margin = margin(b = 10)),
  legend.text = element_text(size = 8),
  legend.key.width = unit(10, "pt")
)

angled_x_axis <- theme(axis.text.x = element_text(angle = 45, hjust = 1))

no_axis_lines <- theme(axis.line = element_blank(),
                       axis.ticks = element_blank())


# Figure 1B --------------------------------------------------------------------

baselines_plot_tissue <- baselines |>
  subset(granularity == "Broad class" & reference_data_name != "All zeros") |>
  Create_BoxStats(c("tissue_full", "error_metric")) |>
  mutate(best_val = case_when(error_metric == "Correlation" ~ max_val,
                              TRUE ~ min_val))

errs_box_noalg <- best_errors |>
  subset(granularity == "Broad class") |>
  Create_BoxStats(c("tissue_full", "error_metric",
                    "normalization", "regression_method"))

# One tissue where regression helps, one where it doesn't
plt1B <- ggplot(subset(errs_box_noalg,
                       tissue_full %in% c("MSBB PHG", "ROSMAP ACC")),
                aes(x = normalization, fill = regression_method,
                    y = median_val, ymin = min_val, ymax = max_val)) +
  geom_hline(aes(yintercept = best_val, color = "Baseline"),
             data = subset(baselines_plot_tissue,
                           tissue_full %in% c("MSBB PHG", "ROSMAP ACC")),
             linetype = "twodash") +
  geom_crossbar(position = position_dodge2(padding = 0),
                fatten = 1.5, width = 0.75) +
  theme_bw() +
  scale_fill_manual(values = regression_colors) +
  scale_color_manual(values = "slategray") +
  facet_grid(error_metric ~ tissue_full, scales = "free", switch = "y") +
  labs(fill = "Regression Method", color = "") +
  xlab(NULL) +
  ylab(NULL) +
  poster_theme +
  theme(panel.spacing.x = unit(10, "pt"),
        panel.spacing.y = unit(20, "pt"))

print(plt1B)

ggsave("plt1B.svg", plt1B, path = file.path("figures", "poster"),
       width = 750, height = 600, units = "px", dpi = 72)


# Figure 1C --------------------------------------------------------------------

plt1C <- ggplot(Count_Grouped(top3_by_tissue,
                              c("tissue", "normalization", "regression_method")),
                aes(x = normalization, y = regression_method,
                    color = count, size = count)) +
  geom_count() + theme_bw() +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5,
                      limits = c(0, 12),
                      breaks = c(0, 3, 6, 9, 12)) +
  facet_wrap(~tissue, nrow = 1) +
  labs(color = "Count") +
  xlab(NULL) +
  ylab(NULL) +
  guides(size = "none") +
  poster_theme + angled_x_axis + numerical_legend

print(plt1C)

ggsave("plt1C.svg", plt1C, path = file.path("figures", "poster"),
       width = 850, height = 170, units = "px", dpi = 72)


# Figure 1D --------------------------------------------------------------------

top3_by_tissue_broad <- subset(top3_by_tissue, granularity == "Broad class")
top3_by_tissue_sub <- subset(top3_by_tissue, granularity == "Sub class")

plt1Da <- ggplot(Count_Grouped(top3_by_tissue_broad,
                               c("tissue", "reference_data_name")),
                 aes(x = tissue, y = reference_data_name, color = count, size = count)) +
  geom_count() +
  theme_classic() +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5,
                      limits = c(1, 9),
                      breaks = c(1, 3, 5, 7, 9)) +
  labs(color = "Count") +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Broad class") +
  guides(size = "none") +
  poster_theme + angled_x_axis + numerical_legend

plt1Db <- ggplot(Count_Grouped(top3_by_tissue_sub,
                               c("tissue", "reference_data_name")),
                 aes(x = tissue, y = reference_data_name, color = count, size = count)) +
  geom_count() + theme_classic() +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5,
                      limits = c(1, 9),
                      breaks = c(1, 3, 5, 7, 9)) +
  labs(color = "Count") +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Sub class") +
  guides(size = "none") +
  poster_theme + angled_x_axis + numerical_legend

print(plt1Da / plt1Db)

ggsave("plt1D.svg", (plt1Da / plt1Db), path = file.path("figures", "poster"),
       width = 400, height = 380, units = "px", dpi = 72)


# Figure 1E --------------------------------------------------------------------

total_valid <- quality_stats$n_valid_by_algorithm |>
  # So something shows up for DeconRNASeq
  mutate(pct_valid = ifelse(pct_valid < 0.01, 0.01, pct_valid))

plt1E <- ggplot(total_valid,
                aes(x = algorithm, y = pct_valid, fill = granularity)) +
  geom_col(position = "dodge2", width = 0.6) +
  theme_bw() +
  scale_fill_manual(values = granularity_colors) +
  labs(fill = NULL) +
  xlab(NULL) +
  ylab("Percent") +
  poster_theme + angled_x_axis

print(plt1E)

ggsave("plt1E.svg", plt1E, path = file.path("figures", "poster"),
       width = 540, height = 250, units = "px", dpi = 72)


# Figure 2A --------------------------------------------------------------------

plt2Aa <- ggplot(Count_Grouped(top3_by_tissue_broad, c("tissue", "algorithm")),
                 aes(x = tissue, y = algorithm, color = count, size = count)) +
  geom_count() +
  theme_classic() +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5) +
  xlab(NULL) +
  ylab(NULL) +
  guides(color = "none", size = "none") +
  ggtitle("Broad class") +
  poster_theme + angled_x_axis + numerical_legend

# Sub class is missing several algorithms that dropped out, fill these in
alg_info <- best_errors |>
  select(tissue, algorithm) |>
  distinct()

alg_sub <- Count_Grouped(top3_by_tissue_sub, c("tissue", "algorithm")) |>
  merge(alg_info, all = TRUE)

plt2Ab <- ggplot(alg_sub,
                 aes(x = tissue, y = algorithm, color = count, size = count)) +
  geom_count() +
  theme_classic() +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5) +
  labs(color = "Count") +
  xlab(NULL) +
  ylab(NULL) +
  guides(size = "none") +
  ggtitle("Sub class") +
  poster_theme + angled_x_axis + numerical_legend

plt2A <- plt2Aa + plt2Ab

print(plt2A)

ggsave("plt2A.svg", plt2A, path = file.path("figures", "poster"),
       width = 700, height = 240, units = "px", dpi = 72)


# Figure 2B --------------------------------------------------------------------

better_stats_alg <- quality_stats$better_than_baseline_by_algorithm |>
  # So something shows up for 0 values
  mutate(pct_better_than_baseline = ifelse(pct_better_than_baseline < 0.01, 0.01,
                                           pct_better_than_baseline))

plt2B <- ggplot(better_stats_alg,
                aes(x = algorithm, y = pct_better_than_baseline, fill = granularity)) +
  geom_col(position = "dodge2") +
  scale_fill_manual(values = granularity_colors) +
  theme_bw() +
  labs(fill = NULL) +
  ylim(0, 1.05) +
  xlab(NULL) +
  ylab("Percent") +
  poster_theme + angled_x_axis

print(plt2B)

ggsave("plt2B.svg", plt2B, path = file.path("figures", "poster"),
       width = 600, height = 350, units = "px", dpi = 72)


# Figure 2C --------------------------------------------------------------------

inh_ratio <- quality_stats$exc_inh_ratio |>
  subset(granularity == "Broad class") |>
  mutate(median_inh_exc_ratio = 1/median_exc_inh_ratio)

plt2C <- ggplot(inh_ratio,
                aes(x = algorithm, y = median_inh_exc_ratio, fill = algorithm)) +
  geom_col() +
  geom_hline(yintercept = 0.5, linetype = "twodash", color = "slategray") +
  theme_bw() +
  facet_wrap(~tissue, nrow = 3) +
  scale_fill_manual(values = algorithm_colors) +
  ggtitle("Median ratio of Inhibitory:Excitatory neurons") +
  xlab(NULL) +
  ylab("Ratio") +
  guides(fill = "none") +
  poster_theme + angled_x_axis

print(plt2C)


# Figure 2? --------------------------------------------------------------------

if (FALSE) {
  bulk_dataset <- "Mayo"
  algorithm <- "Music"
  err_files <- Get_ErrorFiles(bulk_dataset, algorithm, granularity)

  marker_summary <- lapply(err_files, function(EF) {
    err_list <- readRDS(EF)

    if (length(err_list) == 0) {
      next
    }

    errs_all <- err_list$means |>
      subset(tissue != "All") |>
      merge(err_list$params, by.x = "param_id", by.y = "row.names") |>
      mutate(marker_combo = paste(marker_type, marker_subtype, marker_input_type),
             marker_combo = str_replace_all(marker_combo, " None", ""),
             algorithm = algorithm) |>
      Paper_Renames() |>
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

    errs_summary <- errs_all |>
      group_by(tissue, tissue_full, data_transform, marker_combo, marker_order,
               n_marker_type, total_markers_used) |>
      dplyr::summarize(Correlation = max(Correlation),
                       RMSE = min(RMSE),
                       MAPE = min(MAPE),
                       .groups = "drop")
    return(errs_summary)
  })

  marker_summary <- do.call(rbind, marker_summary)

  markers_plot <- best_errors_broad

  # TODO
}



# Figure 3A --------------------------------------------------------------------

significance_broad <- quality_stats$significance |>
  subset(granularity == "Broad class" & tissue != "CBE") |>
  Load_Significance(best_dt, p_sig = 0.01, log2_cap = 1)

plot_limit <- c(-1, 1)
n_cols <- 4 #if (granularity == "Broad class") 4 else 6

ok <- significance_broad |>
  group_by(celltype) |>
  summarize(n_significant = sum(is.finite(log2_fc))) |>
  subset(n_significant >= 3)
sig_final_ok <- subset(significance_broad, celltype %in% unique(ok$celltype))

plt3A <- ggplot(sig_final_ok, aes(x = algorithm, y = tissue, fill = log2_fc)) +
  geom_tile(color = "black") +
  theme_classic() +
  facet_wrap(~celltype, ncol = 4) +
  coord_fixed() +
  scale_fill_distiller(palette = "RdBu", limit = plot_limit, na.value = "white") +
  scale_y_discrete(limits = rev) +
  xlab(NULL) + ylab(NULL) +
  labs(fill = "log2-FC") +
  poster_theme + angled_x_axis + numerical_legend + no_axis_lines +
  theme(panel.spacing = unit(10, "pt"))

print(plt3A)

# Sizing is set up to make each grid ~130x150 px
ggsave("plt3A.svg", plt3A, path = file.path("figures", "poster"),
       width = 720, height = 260, units = "px", dpi = 72)


# Figure 3B --------------------------------------------------------------------

significance_sub <- quality_stats$significance |>
  subset(granularity == "Sub class" & tissue != "CBE") |>
  Load_Significance(best_dt, p_sig = 0.01, log2_cap = 1)

plot_limit <- c(-1, 1)
n_cols <- 6

ok <- significance_sub |>
  group_by(celltype) |>
  summarize(n_significant = sum(is.finite(log2_fc))) |>
  subset(n_significant >= 3)
sig_final_ok <- subset(significance_sub, celltype %in% unique(ok$celltype))

plt3B <- ggplot(sig_final_ok, aes(x = algorithm, y = tissue, fill = log2_fc)) +
  geom_tile(color = "black") + theme_classic() +
  facet_wrap(~celltype, ncol = 5) +
  coord_fixed() +
  scale_fill_distiller(palette = "RdBu", limit = plot_limit, na.value = "white") +
  xlab(NULL) + ylab(NULL) +
  labs(fill = "log2-FC") +
  poster_theme + angled_x_axis + numerical_legend + no_axis_lines +
  theme(panel.spacing = unit(10, "pt"))

print(plt3B)

# Sizing is set up to make each grid ~130x150 px
ggsave("plt3B.svg", plt3B, path = file.path("figures", "poster"),
       width = 850, height = 450, units = "px", dpi = 72)


# Plot 3C ----------------------------------------------------------------------

# From differential expression of each cell type vs all other excitatory
# neurons, plus markers from the Cain paper
markers <- list(
  "Exc.1" = c("CUX2", "GLIS3", "CARTPT"),
  "Exc.2" = c("TSHZ2", "GABRG1", "TDRD1", "RORB"),
  "Exc.4" = c("RSPO3", "TMEM212", "TMSB10", "RORB"),
  "Exc.7" = c("MCUB", "TMEM233", "PRRX1", "THEMIS"),
  "Exc.10" = c("POSTN", "SMYD1", "RGS12", "THEMIS")
)

genes <- unlist(markers) |> unique()

#pb <- Load_PseudobulkPureSamples("cain", "sub_class", output_type = "log_cpm")
#pb_mat <- assay(pb, "counts")[genes, pb$celltype %in% names(markers)]
#pb_mat <- scale(t(pb_mat)) |> t()

sig <- Load_SignatureMatrix("cain", "sub_class", "log_cpm")
sig <- sig[genes, names(markers)]

pb_plot <- pb_mat |>
  as.data.frame() |>
  tibble::rownames_to_column("gene") |>
  tidyr::pivot_longer(-gene, names_to = "sample", values_to = "expression") |>
  merge(as.data.frame(colData(pb)))

ggplot(pb_plot, aes(x = sample, y = gene, fill = expression)) +
  geom_tile() +
  theme_classic() +
  facet_wrap(~celltype, nrow = 1) +
  scale_fill_viridis() +
  poster_theme + no_axis_lines +
  theme(axis.text.x = element_blank())

anno <- subset(colData(pb), sample %in% colnames(pb_mat)) |>
  as.data.frame() |>
  select(celltype)

pheatmap::pheatmap(t(sig),
                   color = viridis(100),
                   show_rownames = TRUE,
                   #show_colnames = FALSE,
                   treeheight_row = 0,
                   treeheight_col = 0,
                   scale = "column")
