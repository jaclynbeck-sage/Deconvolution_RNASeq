library(Matrix)
library(ggplot2)
library(viridis)
library(dplyr)
library(stringr)
library(reshape2)
library(patchwork)

source(file.path("functions", "Step11_Error_HelperFunctions.R"))
source(file.path("functions", "Plotting_HelperFunctions.R"))

options(scipen = 999)

datasets <- c("cain", "lau", "leng", "mathys", "seaRef")

est_fields = list("CibersortX" = "estimates",
                  "DeconRNASeq" = "estimates",
                  "Dtangle" = "estimates",
                  "DWLS" = "estimates",
                  "HSPE" = "estimates",
                  "Music" = "Est.pctRNA.weighted",
                  "Scaden" = "estimates",
                  "Baseline" = "estimates")

algorithms <- names(est_fields)

granularity <- c("broad_class")

bulk_datasets <- c("Mayo", "MSBB", "ROSMAP")

best_errors_list <- Get_AllBestErrorsAsDf(bulk_datasets, granularity)
best_errors <- best_errors_list$errors
quality_stats <- best_errors_list$quality_stats

best_errs_filt <- best_errors %>%
  mutate(reference_data_name = str_replace(reference_data_name, "_.*", ""),
         reference_data_name = if_else(algorithm == "Baseline",
                                       reference_data_name, "single cell"),
         tissue = paste(test_data_name, tissue),
         regression_method = str_replace(regression_method, "none", "no regression"),
         normalization = str_replace(normalization, "counts", "cpm"),
         normalization = str_replace(normalization, "cpm", "counts/cpm"),
         normalization = str_replace(normalization, "log_", "")) %>%
  subset(!(grepl("All", tissue)))

baselines <- best_errs_filt$algorithm == "Baseline"
best_errs_filt$algorithm[baselines] <- paste0(best_errs_filt$algorithm[baselines],
                                              " (", best_errs_filt$reference_data_name[baselines],
                                              ")")

errs_melt <- best_errs_filt  %>%
  group_by(tissue, solve_type, reference_data_name,
           test_data_name, algorithm, normalization,
           regression_method) %>%
  summarize(cor = max(cor),
            rMSE = min(rMSE),
            mAPE = min(mAPE),
            .groups = "drop") %>%
  melt(variable.name = "error_type")

errs_melt$data_transform <- paste(errs_melt$normalization, "/", errs_melt$regression_method)
errs_melt$data_transform <- str_replace(errs_melt$data_transform, "log_", "")

# TODO temporary
#errs_melt <- subset(errs_melt, !grepl("Baseline", algorithm))

##### Color setup #####

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
normalization_colors <- RColorBrewer::brewer.pal(length(regs), "Set1")

# Regression methods
regression_colors <- viridis::turbo(length(regs), begin = 0.1, end = 0.9)#RColorBrewer::brewer.pal(length(regs), "Accent")

# Tissues (modified viridis turbo color scheme)
tiss <- colSums(table(errs_melt$tissue, errs_melt$test_data_name) > 1)
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

params <- list(solve_type = "signature",
               normalization = unique(errs_melt$normalization),
               regression_method = unique(errs_melt$regression_method))

plt1 <- Plot_ErrsByAlgorithm(errs_melt, params, "cor", fill = "tissue",
                             fill_colors = tissue_colors,
                             facet_var = c("test_data_name", "algorithm"))

plt2 <- ggplot(subset(errs_melt, error_type == "cor"),
               aes(x = algorithm, y = value, fill = tissue)) +
  geom_boxplot() + theme_bw() + scale_fill_manual(values = tissue_colors) +
  facet_wrap(~test_data_name) #+ geom_errorbar(aes(ymin = min_val, ymax = max_val, group = tissue))

############ Baseline vs best of any param set for each algorithm ##############

box_stats <- errs_melt %>%
  group_by(tissue, test_data_name, algorithm, error_type) %>%
  summarize(max_val = max(value),
            min_val = min(value),
            median_val = median(value),
            upper_quartile = quantile(value, probs = 0.75, na.rm = TRUE),
            lower_quartile = quantile(value, probs = 0.25, na.rm = TRUE),
            .groups = "drop")

plt3 <- ggplot(subset(box_stats, error_type == "cor" & algorithm != "Baseline (zeros)"),
               aes(x = algorithm, ymin = min_val, ymax = max_val,
                   lower = lower_quartile, middle = median_val,
                   upper = upper_quartile, fill = tissue)) +
  geom_boxplot(stat = "identity") +
  theme_bw() + scale_fill_manual(values = tissue_colors) +
  facet_wrap(~test_data_name) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

plt4 <- ggplot(subset(box_stats, error_type == "rMSE"),
               aes(x = algorithm, ymin = min_val, ymax = max_val,
                   lower = lower_quartile, middle = median_val,
                   upper = upper_quartile, fill = tissue)) +
  geom_boxplot(stat = "identity", ) +
  theme_bw() + scale_fill_manual(values = tissue_colors) +
  facet_wrap(~test_data_name) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

plt5 <- ggplot(subset(box_stats, error_type == "mAPE" & algorithm != "Baseline (zeros)"),
               aes(x = algorithm, ymin = min_val, ymax = max_val,
                   lower = lower_quartile, middle = median_val,
                   upper = upper_quartile, fill = tissue)) +
  geom_boxplot(stat = "identity", ) +
  theme_bw() + scale_fill_manual(values = tissue_colors) +
  facet_wrap(~test_data_name) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))


########### Normalization vs regression ##################

errs_melt2 <- melt(best_errs_filt, variable.name = "error_type")  %>%
  group_by(tissue, solve_type, test_data_name, algorithm, error_type,
           normalization, regression_method) %>%
  summarize(max_val = max(value),
            min_val = min(value),
            median_val = median(value),
            upper_quartile = quantile(value, probs = 0.75, na.rm = TRUE),
            lower_quartile = quantile(value, probs = 0.25, na.rm = TRUE),
            .groups = "drop")

# TODO we need to split by tissue, but duing a tissue x algorithm grid is pretty large.
# How to collapse it a bit?
plt6 <- ggplot(subset(errs_melt2, error_type == "cor" & algorithm != "Baseline (zeros)"),
               aes(x = normalization, ymin = min_val, ymax = max_val,
                   lower = lower_quartile, middle = median_val,
                   upper = upper_quartile, fill = regression_method)) +
  geom_boxplot(stat = "identity") +
  theme_bw() + scale_fill_manual(values = regression_colors) +
  facet_grid(rows = vars(tissue), cols = vars(algorithm)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

plt7 <- ggplot(subset(errs_melt2, error_type == "rMSE"),
               aes(x = normalization, ymin = min_val, ymax = max_val,
                   lower = lower_quartile, middle = median_val,
                   upper = upper_quartile, fill = regression_method)) +
  geom_boxplot(stat = "identity") +
  theme_bw() + scale_fill_manual(values = regression_colors) +
  facet_grid(rows = vars(tissue), cols = vars(algorithm)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

plt8 <- ggplot(subset(errs_melt2, error_type == "mAPE" & algorithm != "Baseline (zeros)"),
               aes(x = normalization, ymin = min_val, ymax = max_val,
                   lower = lower_quartile, middle = median_val,
                   upper = upper_quartile, fill = regression_method)) +
  geom_boxplot(stat = "identity") +
  theme_bw() + scale_fill_manual(values = regression_colors) +
  facet_grid(rows = vars(tissue), cols = vars(algorithm)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))



for (err_metric in c("cor", "rMSE", "mAPE")) {
  errs_sub <- subset(melt(best_errors, variable.name = "error_type"),
                     solve_type == "signature" & error_type == err_metric &
                       tissue == "All" & grepl(err_metric, metrics))
  errs_sub$normalization <- str_replace(errs_sub$normalization, "log_", "")
  errs_sub$normalization <- str_replace(errs_sub$normalization, "counts", "cpm")

  ggplot(errs_sub, aes(x = normalization, y = value, fill = regression_method)) +
    geom_boxplot() + theme_bw() + facet_wrap(~test_data_name, ncol = 3) +
    #facet_grid(rows = vars(test_data_name), cols = vars(tissue)) +
    scale_fill_manual(values = regression_colors)

  errs_sub2 <- subset(melt(best_errors, variable.name = "error_type"),
                      test_data_name == "Mayo" &
                        solve_type == "signature" & error_type == err_metric &
                        tissue != "All" & grepl(err_metric, metrics))
  errs_sub2$normalization <- str_replace(errs_sub2$normalization, "log_", "")
  errs_sub2$normalization <- str_replace(errs_sub2$normalization, "counts", "cpm")

  ggplot(errs_sub2, aes(x = normalization, y = value, fill = regression_method)) +
    geom_boxplot() + theme_bw() + #facet_wrap(~tissue)
    facet_grid(rows = vars(tissue), cols = vars(algorithm)) +
    scale_fill_manual(values = regression_colors)

  errs_sub3 <- subset(melt(best_errors, variable.name = "error_type"),
                      test_data_name == "MSBB" &
                        solve_type == "signature" & error_type == err_metric &
                        tissue != "All" & grepl(err_metric, metrics))
  errs_sub3$normalization <- str_replace(errs_sub3$normalization, "log_", "")
  errs_sub3$normalization <- str_replace(errs_sub3$normalization, "counts", "cpm")

  ggplot(errs_sub3, aes(x = normalization, y = value, fill = regression_method)) +
    geom_boxplot() + theme_bw() + #facet_wrap(~tissue)
    facet_grid(rows = vars(tissue), cols = vars(algorithm)) +
    scale_fill_manual(values = regression_colors)

  errs_sub4 <- subset(melt(best_errors, variable.name = "error_type"),
                      test_data_name == "ROSMAP" &
                        solve_type == "signature" & error_type == err_metric &
                        tissue != "All" & grepl(err_metric, metrics))
  errs_sub4$normalization <- str_replace(errs_sub4$normalization, "log_", "")
  errs_sub4$normalization <- str_replace(errs_sub4$normalization, "counts", "cpm")

  ggplot(errs_sub4, aes(x = normalization, y = value, fill = regression_method)) +
    geom_boxplot() + theme_bw() + #facet_wrap(~tissue)
    facet_grid(rows = vars(tissue), cols = vars(algorithm)) +
    scale_fill_manual(values = regression_colors)
}


quality_stats$normalization <- str_replace(quality_stats$normalization, "log_", "")
quality_stats$normalization <- str_replace(quality_stats$normalization, "counts", "cpm")
quality_stats <- subset(quality_stats, algorithm != "Baseline")

############ Bad inhibitory ratio, norm vs regression ##############

qbox_stats <- quality_stats %>% group_by(test_data_name, algorithm) %>%
  summarize(max_val = max(mean_pct_bad_inhibitory_ratio),
            min_val = min(mean_pct_bad_inhibitory_ratio),
            median_val = median(mean_pct_bad_inhibitory_ratio),
            upper_quartile = quantile(mean_pct_bad_inhibitory_ratio, probs = 0.75, na.rm = TRUE),
            lower_quartile = quantile(mean_pct_bad_inhibitory_ratio, probs = 0.25, na.rm = TRUE),
            .groups = "drop")

ggplot(quality_stats, aes(x = normalization, y = mean_pct_bad_inhibitory_ratio,
                          fill = regression_method)) +
  geom_boxplot() + theme_bw() + #facet_wrap(~tissue)
  facet_grid(rows = vars(test_data_name), cols = vars(algorithm)) +
  scale_fill_manual(values = regression_colors)

ggplot(qbox_stats, aes(x = algorithm, ymin = min_val, ymax = max_val,
                          lower = lower_quartile, middle = median_val,
                          upper = upper_quartile, fill = algorithm)) +
  geom_boxplot(stat = "identity") +
  theme_bw() + facet_wrap(~test_data_name) + #, ncol = 1) +
  #facet_grid(rows = vars(test_data_name), cols = vars(algorithm)) +
  scale_fill_manual(values = algorithm_colors)

########### params passing QC, norm vs regression #######################

# Fill in stats for reference data sets that didn't have any valid param sets
# and therefore aren't in this data frame
file_params <- quality_stats %>% subset(algorithm == "CibersortX") %>%
  select(test_data_name, granularity, normalization,
         regression_method, algorithm) %>%
  distinct()

for (row in 1:nrow(file_params)) {
  tmp <- subset(quality_stats,
                test_data_name == file_params$test_data_name[row] &
                  granularity == file_params$granularity[row] &
                  normalization == file_params$normalization[row] &
                  regression_method == file_params$regression_method[row] &
                  algorithm == "CibersortX")
  missing <- setdiff(datasets, tmp$reference_data_name)
  for (M in missing) {
    new_df <- cbind(reference_data_name = M, file_params[row,])
    new_df$mean_pct_bad_inhibitory_ratio <- 0
    new_df$file_id <- "NA"
    new_df$pct_valid <- 0

    quality_stats <- rbind(quality_stats, new_df[,colnames(quality_stats)])
  }
}

ggplot(quality_stats, aes(x = normalization, y = pct_valid,
                          fill = regression_method)) +
  geom_boxplot() + theme_bw() + #facet_wrap(~tissue)
  facet_grid(rows = vars(test_data_name), cols = vars(algorithm)) +
  scale_fill_manual(values = regression_colors)

qbox_stats2 <- quality_stats %>% group_by(test_data_name, algorithm) %>%
  summarize(max_val = max(pct_valid),
            min_val = min(pct_valid),
            median_val = median(pct_valid),
            upper_quartile = quantile(pct_valid, probs = 0.75, na.rm = TRUE),
            lower_quartile = quantile(pct_valid, probs = 0.25, na.rm = TRUE),
            .groups = "drop")

ggplot(qbox_stats2, aes(x = algorithm, ymin = min_val, ymax = max_val,
                       lower = lower_quartile, middle = median_val,
                       upper = upper_quartile, fill = algorithm)) +
  geom_boxplot(stat = "identity") +
  theme_bw() + facet_wrap(~test_data_name) + #, ncol = 1) +
  #facet_grid(rows = vars(test_data_name), cols = vars(algorithm)) +
  scale_fill_manual(values = algorithm_colors)

file_params <- best_errors %>%
  select(reference_data_name, test_data_name, algorithm,
         reference_input_type, normalization,
         regression_method) %>%
  distinct()

# TODO hspe is wrong
for (dataset in datasets) {
  pdf(file.path(dir_figures, paste0("error_plots_", dataset, "_broad_detailed.pdf")), width=10, height = 12)

  for (bulk_dataset in bulk_datasets) {
    bulk_se <- Load_BulkData(bulk_dataset)
    metadata <- as.data.frame(colData(bulk_se))

    # TODO other tissues and lm?
    best_errors_sub <- subset(best_errors, tissue == "DLPFC" &
                                solve_type == "signature" &
                                regression_method == "edger" &
                                normalization %in% c("counts", "cpm", "log_cpm") &
                                test_data_name == bulk_dataset)

    best_ests <- Get_AllEstimatesAsDf(dataset, bulk_dataset, algorithms,
                                      granularity, best_errors_sub, est_fields)

    # TODO fix
    #best_params <- subset(best_errors, reference_data_name == dataset &
    #                        algorithm %in% unique(ests_alg$algorithm) &
    #                        grepl(bulk_dataset, test_data_name))
    #best_params <- best_params %>% group_by(algorithm) %>%
    #  mutate(#title = paste("Best", metrics, "/", test_datasets), # TODO fix
    #         title_short = paste(algorithm, "params", as.numeric(factor(title)))) %>%
    #  ungroup()
    #best_params$param_id <- sapply(best_params$params, paste, collapse = " ")

    #best_ests <- ests_alg %>% merge(best_params, by = c("algorithm", "param_id"))
    best_ests <- best_ests %>% group_by(algorithm) %>%
      mutate(title = paste(algorithm, "params",
                           as.numeric(factor(param_id))),
             title_short = title) %>%
      ungroup()

    best_ests_dlpfc <- subset(best_ests, sample %in% metadata$sample[metadata$tissue == "DLPFC"])

    plot_id <- paste("Reference dataset:", dataset, "/ Bulk dataset:", bulk_dataset,
                     "/ Tissue: DLPFC / Normalization: CPM / Regression: edgeR")

    merscope <- data.frame(celltype = c("Astrocyte", "Excitatory", "Inhibitory",
                                        "Microglia", "Oligodendrocyte", "OPC",
                                        "Vascular"),
                           pct = c(10, 34, 7, 6, 32.5, 3.5, 11.5)/100)
    merscope = subset(merscope, celltype %in% c("Astrocyte", "Oligodendrocyte"))
    merscope = rbind(merscope, data.frame(celltype = c("Endothelial", "Pericyte"),
                                          pct = c(8, 3.5)/100))

    plt <- ggplot(best_ests_dlpfc, aes(x = algorithm, y = pct_est, color = title)) +
      geom_boxplot(width = 0.5) + theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      theme(legend.position = "bottom") +
      facet_wrap(~celltype, scales = "fixed", nrow = 3) +
      geom_hline(data = merscope, mapping = aes(yintercept = pct), color = "red")
    #ncol = ceiling(length(unique(best_ests$celltype)) / 2))

    print(plt + plot_annotation(title = plot_id, subtitle = "Estimated proportions across all subjects (each bar is a single parameter set)"))

    # TODO make functions to do this and move these files to a better place
    ests_ad <- merge(best_ests, metadata, by = "sample") #%>%
    #subset(diagnosis %in% c("CT", "AD"))

    ests_ad$celltype <- factor(ests_ad$celltype)

    significant <- list()
    for (param_set in unique(ests_ad$title_short)) {
      ests_param <- subset(ests_ad, title_short == param_set)
      anov <- aov(pct_est ~ diagnosis*celltype*tissue, data = ests_param)
      tuk <- TukeyHSD(anov, "diagnosis:celltype:tissue")

      ct_v_tissue <- expand.grid(levels(ests_ad$celltype), levels(ests_ad$tissue))
      ct_v_tissue <- paste(ct_v_tissue$Var1, ct_v_tissue$Var2, sep = ":")

      comparisons <- paste0("CT:", ct_v_tissue, "-AD:", ct_v_tissue)
      tuk <- as.data.frame(tuk[[1]][comparisons,])
      tuk$p_adj <- tuk$`p adj`
      tuk$celltype <- str_split(rownames(tuk), pattern = ":", simplify = TRUE)[,2]
      tuk$tissue <- str_split(rownames(tuk), pattern = ":", simplify = TRUE)[,5]
      tuk$significant <- tuk$p_adj <= 0.05
      tuk$title_short <- param_set
      significant[[param_set]] <- tuk
    }

    significant <- do.call(rbind, significant) %>% select(p_adj, celltype, tissue, significant, title_short)

    ests_ad <- merge(ests_ad, significant, by = c("celltype", "tissue", "title_short"))

    for (alg in unique(ests_ad$algorithm)) {
      ests_params <- subset(ests_ad, algorithm == alg & diagnosis %in% c("CT", "AD"))
      plt <- ggplot(ests_params, aes(x = celltype, y = pct_est, fill = diagnosis, color = significant), group = tissue) +
        geom_boxplot(width = 0.5) + theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_color_manual(values = c("#dddddd", "#000000")) +
        facet_wrap(~title_short + tissue, scales = "free",
                   ncol = length(unique(ests_ad$tissue)))
      print(plt + plot_annotation(title = plot_id, subtitle = paste(alg, "estimates, AD vs Control")))
    }
    # TODO can we break down goodness-of-fit errors by cell type and subject?
  }

  dev.off()
}
