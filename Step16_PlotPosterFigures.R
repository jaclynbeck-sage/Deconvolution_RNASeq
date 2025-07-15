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


# Color setup ------------------------------------------------------------------

# Algorithms
algs <- unique(best_errors$algorithm)
algorithm_colors <- sage_hue_pal(level = "400")(length(algs))
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

# ggplot2 settings that generally apply to all plots
poster_theme <- theme(
  # All text
  text = element_text(family = "DM Sans"),
  # Add space between axis titles and axis text
  axis.title.x = element_text(size = 20, margin = margin(t = 10, unit = "pt")),
  axis.title.y = element_text(size = 20, margin = margin(r = 10, unit = "pt")),
  # Larger axis labels
  axis.text = element_text(size = 16),
  # Plot title
  plot.title = element_text(face = "bold", size = 24),
  # Facet labels
  strip.background = element_blank(),
  strip.text = element_text(size = 24, face = "bold"),
  strip.placement.y = "outside",
  # Legend text
  legend.title = element_text(size = 16, margin = margin(b = 10)),
  legend.text = element_text(size = 16),
)

dotplot_theme <- theme(
  axis.text = element_text(size = 14),
  title = element_text(size = 18),
  legend.text = element_text(size = 14)
)

x_axis_45 <- theme(axis.text.x = element_text(angle = 45, hjust = 1))
x_axis_90 <- theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

no_axis_lines <- theme(axis.line = element_blank(),
                       axis.ticks = element_blank())


# Figure 1B --------------------------------------------------------------------

baselines_plot_tissue <- baselines |>
  subset(reference_data_name != "All zeros" & error_metric == "Correlation") |>
  Create_BoxStats(c("tissue_full", "error_metric")) |>
  mutate(best_val = case_when(error_metric == "Correlation" ~ max_val,
                              TRUE ~ min_val))

errs_box_noalg <- best_errors |>
  subset(error_metric == "Correlation") |>
  Create_BoxStats(c("tissue_full", "error_metric",
                    "normalization", "regression_method"))

# One tissue where regression helps, one where it doesn't
plt1B <- ggplot(subset(errs_box_noalg,
                       tissue_full %in% c("Mayo TCX", "MSBB PHG", "ROSMAP ACC")),
                aes(x = normalization, fill = regression_method,
                    y = median_val, ymin = min_val, ymax = max_val)) +
  geom_hline(aes(yintercept = best_val, color = "Baseline"),
             data = subset(baselines_plot_tissue,
                           tissue_full %in% c("Mayo TCX", "MSBB PHG", "ROSMAP ACC")),
             linetype = "twodash") +
  geom_crossbar(position = position_dodge2(padding = 0),
                fatten = 1.5, width = 0.75) +
  theme_bw() +
  scale_fill_manual(values = regression_colors) +
  scale_color_manual(values = sage_colors$blueberry[["700"]]) +
  facet_wrap(~tissue_full) +
  labs(fill = "Regression Method", color = "") +
  xlab(NULL) +
  ylab("Correlation") +
  poster_theme +
  theme(panel.spacing.x = unit(20, "pt")) #,
  #      panel.spacing.y = unit(20, "pt"))

print(plt1B)

ggsave("plt1B.svg", plt1B, path = file.path("figures", "poster"),
       width = 800, height = 240, units = "px", dpi = 72)


# Figure 1C --------------------------------------------------------------------
# Not used

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
  poster_theme + x_axis_45

print(plt1C)

#ggsave("plt1C.svg", plt1C, path = file.path("figures", "poster"),
#       width = 850, height = 170, units = "px", dpi = 72)


# Figure 1D --------------------------------------------------------------------
# Not used

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
  poster_theme + x_axis_45

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
  ggtitle("Subtypes") +
  guides(size = "none") +
  poster_theme + x_axis_45

print(plt1Da / plt1Db)

#ggsave("plt1D.svg", (plt1Da / plt1Db), path = file.path("figures", "poster"),
#       width = 400, height = 380, units = "px", dpi = 72)


# Figure 1E --------------------------------------------------------------------

total_valid <- quality_stats$n_valid_by_algorithm |>
  # So something shows up for DeconRNASeq
  mutate(pct_valid = ifelse(pct_valid < 0.01, 0.01, pct_valid),
         granularity = ifelse(granularity == "Sub class", "Subtypes", granularity))

plt1E <- ggplot(total_valid,
                aes(x = algorithm, y = pct_valid, fill = granularity)) +
  geom_col(position = "dodge2", width = 0.7) +
  theme_bw() +
  scale_fill_manual(values = granularity_colors) +
  labs(fill = NULL) +
  xlab(NULL) +
  ylab("Percent") +
  poster_theme + x_axis_90 +
  theme(legend.position = "bottom")

print(plt1E)

ggsave("plt1E.svg", plt1E, path = file.path("figures", "poster"),
       width = 315, height = 285, units = "px", dpi = 72)


# Figure 2A --------------------------------------------------------------------

best_algs <- Count_Grouped(top3_by_tissue_broad, c("tissue", "algorithm")) |>
  subset(tissue != "CBE")

plt2Aa <- ggplot(best_algs,
                 aes(x = tissue, y = algorithm, color = count, size = count)) +
  geom_count() +
  theme_classic() +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5) +
  scale_y_discrete(limits = rev) +
  xlab(NULL) +
  ylab(NULL) +
  guides(color = "none", size = "none") +
  ggtitle("Broad class") +
  poster_theme + x_axis_90 + dotplot_theme

# Sub class is missing several algorithms that dropped out, fill these in
alg_info <- best_errors |>
  select(tissue, algorithm) |>
  distinct()

alg_sub <- Count_Grouped(top3_by_tissue_sub, c("tissue", "algorithm")) |>
  merge(alg_info, all = TRUE) |>
  subset(tissue != "CBE")

plt2Ab <- ggplot(alg_sub,
                 aes(x = tissue, y = algorithm, color = count, size = count)) +
  geom_count() +
  theme_classic() +
  coord_fixed() +
  scale_color_viridis(option = "plasma", direction = -1, begin = 0.5) +
  scale_y_discrete(limits = rev) +
  labs(color = "Count") +
  xlab(NULL) +
  ylab(NULL) +
  guides(size = "none") +
  ggtitle("Subtypes") +
  poster_theme + x_axis_90 + dotplot_theme

plt2A <- plt2Aa + plt2Ab

print(plt2A)

ggsave("plt2A.svg", plt2A, path = file.path("figures", "poster"),
       width = 620, height = 230, units = "px", dpi = 72)


# Figure 2B --------------------------------------------------------------------

better_stats_alg <- quality_stats$better_than_baseline_by_algorithm |>
  # So something shows up for 0 values
  mutate(pct_better_than_baseline = ifelse(pct_better_than_baseline < 0.01, 0.01,
                                           pct_better_than_baseline),
         granularity = ifelse(granularity == "Sub class", "Subtypes", granularity))

plt2B <- ggplot(better_stats_alg,
                aes(x = algorithm, y = pct_better_than_baseline, fill = granularity)) +
  geom_col(position = "dodge2", width = 0.7) +
  scale_fill_manual(values = granularity_colors) +
  theme_bw() +
  labs(fill = NULL) +
  ylim(0, 1.05) +
  xlab(NULL) +
  ylab("Percent") +
  poster_theme + x_axis_90 +
  theme(legend.position = "bottom")

print(plt2B)

ggsave("plt2B.svg", plt2B, path = file.path("figures", "poster"),
       width = 315, height = 285, units = "px", dpi = 72)


# Figure 2C --------------------------------------------------------------------

inh_ratio <- quality_stats$exc_inh_ratio |>
  subset(tissue %in% c("TCX", "PHG", "ACC")) |>
  mutate(median_inh_exc_ratio = 1/median_exc_inh_ratio,
         granularity = str_replace(granularity, "Sub class", "Subtypes"),
         tissue_full = case_match(tissue,
                                  "TCX" ~ "Mayo TCX",
                                  "PHG" ~ "MSBB PHG",
                                  "ACC" ~ "ROSMAP ACC"))

plt2C <- ggplot(inh_ratio,
                aes(x = algorithm, y = median_inh_exc_ratio, fill = algorithm)) +
  geom_col() +
  geom_hline(yintercept = 0.5, linetype = "twodash", color = sage_colors$blueberry[["700"]]) +
  theme_bw() +
  #facet_wrap(~tissue_full, nrow = 1) +
  facet_grid(tissue_full ~ granularity) +
  scale_fill_manual(values = algorithm_colors) +
  scale_y_continuous(limits = c(0, 1)) +
  xlab(NULL) +
  ylab("Ratio") +
  guides(fill = "none") +
  poster_theme + x_axis_90 +
  theme(panel.spacing = unit(20, "pt"),
        strip.text = element_text(size = 23))

print(plt2C)

ggsave("plt2C.svg", plt2C, path = file.path("figures", "poster"),
       width = 450, height = 650, units = "px", dpi = 72)


# Figure 3A --------------------------------------------------------------------

significance_broad <- quality_stats$significance |>
  subset(granularity == "Broad class" & tissue != "CBE" & algorithm != "Baseline") |>
  Load_Significance(best_dt, p_sig = 0.01, log2_cap = 1)

plot_limit <- c(-1, 1)
n_cols <- 4

ok <- significance_broad |>
  group_by(celltype) |>
  summarize(n_significant = sum(is.finite(log2_fc))) |>
  subset(n_significant >= 3)
sig_final_ok <- subset(significance_broad, celltype %in% unique(ok$celltype)) |>
  mutate(celltype = ifelse(celltype == "Oligodendrocyte", "Oligo.", celltype))

plt3A <- ggplot(sig_final_ok, aes(x = algorithm, y = tissue, fill = log2_fc)) +
  geom_tile(color = "black") +
  theme_classic() +
  facet_wrap(~celltype, ncol = 4) +
  coord_fixed() +
  scale_fill_distiller(palette = "RdBu", limit = plot_limit, na.value = "white") +
  scale_y_discrete(limits = rev) +
  xlab(NULL) + ylab(NULL) +
  labs(fill = "log2-FC") +
  poster_theme + x_axis_90 + no_axis_lines +
  theme(panel.spacing = unit(10, "pt"),
        strip.text = element_text(size = 18, face = "bold"))

print(plt3A)

# Sizing is set up to make each grid ~130x170 px
ggsave("plt3A.svg", plt3A, path = file.path("figures", "poster"),
       width = 740, height = 335, units = "px", dpi = 72)


# Figure 3B --------------------------------------------------------------------

significance_sub <- quality_stats$significance |>
  subset(granularity == "Sub class" & tissue != "CBE" & algorithm != "Baseline") |>
  Load_Significance(best_dt, p_sig = 0.01, log2_cap = 1)

plot_limit <- c(-1, 1)
n_cols <- 6

lev <- levels(significance_sub$celltype) |> str_replace("Oligodendrocyte", "Oligo.")

ok <- significance_sub |>
  group_by(celltype) |>
  summarize(n_significant = sum(is.finite(log2_fc))) |>
  subset(n_significant >= 3)
sig_final_ok <- subset(significance_sub, celltype %in% unique(ok$celltype)) |>
  mutate(celltype = ifelse(celltype == "Oligodendrocyte", "Oligo.", as.character(celltype)),
         celltype = factor(celltype, levels = lev))

plt3B <- ggplot(sig_final_ok, aes(x = algorithm, y = tissue, fill = log2_fc)) +
  geom_tile(color = "black") + theme_classic() +
  facet_wrap(~celltype, ncol = 5) +
  coord_fixed() +
  scale_fill_distiller(palette = "RdBu", limit = plot_limit, na.value = "white") +
  scale_y_discrete(limits = rev) +
  xlab(NULL) + ylab(NULL) +
  labs(fill = "log2-FC") +
  poster_theme + x_axis_90 + no_axis_lines +
  theme(panel.spacing = unit(10, "pt"),
        strip.text = element_text(size = 18, face = "bold"))

print(plt3B)

# Sizing is set up to make each grid ~130x150 px
ggsave("plt3B.svg", plt3B, path = file.path("figures", "poster"),
       width = 900, height = 530, units = "px", dpi = 72)


# Plot 3C ----------------------------------------------------------------------

# From differential expression of each cell type vs all other excitatory
# neurons, plus markers from the Cain paper
markers <- list(
  "Exc.1" = c("CUX2", "CARTPT", "GLIS3"),
  "Exc.2" = c("TSHZ2", "GABRG1", "TDRD1", "RORB"),
  "Exc.4" = c("TMEM212", "TMSB10", "RSPO3", "RORB"),
  "Exc.7" = c("MCUB", "TMEM233", "PRRX1", "THEMIS"),
  "Exc.10" = c("SMYD1", "POSTN", "RGS12", "THEMIS")
)

genes <- unlist(markers) |> unique()

sigs <- lapply(c("cain", "lau", "leng", "mathys", "seaRef"), function(ds) {
  sig <- Load_SignatureMatrix(ds, "sub_class", "log_cpm")
  sig[genes, names(markers)] |>
    as.data.frame() |>
    tibble::rownames_to_column("gene") |>
    tidyr::pivot_longer(-gene, names_to = "celltype", values_to = "expression")
})

# Scale mean expression for each gene
sigs <- do.call(rbind, sigs) |>
  group_by(gene, celltype) |>
  summarize(expression = mean(expression), .groups = "drop") |>
  group_by(gene) |>
  mutate(expression = as.numeric(scale(expression))) |>
  ungroup()

gene_order <- genes
cell_order <- c("Exc.1", "Exc.2", "Exc.4", "Exc.7", "Exc.10")

sig_plot <- sigs |>
  mutate(gene = factor(gene, levels = gene_order),
         celltype = factor(celltype, levels = cell_order))

plt3C <- ggplot(sig_plot, aes(x = celltype, y = gene, fill = expression)) +
  geom_tile(color = "lightgray") +
  theme_classic() +
  scale_fill_viridis() +
  scale_y_discrete(limits = rev) +
  labs(x = NULL, y = NULL, fill = str_wrap("Scaled log2-CPM", width = 10)) +
  poster_theme + no_axis_lines + x_axis_90 +
  #theme(legend.title = element_text(size = 12)) +
  coord_fixed()

print(plt3C)

ggsave("plt3C.svg", plt3C, path = file.path("figures", "poster"),
       width = 320, height = 430, units = "px", dpi = 72)

# TODO is inh:exc ratio percent RNA comparison or percent cells comparison? this might matter.
