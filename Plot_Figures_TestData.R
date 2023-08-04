library(Matrix)
library(ggplot2)
library(viridis)
library(dplyr)
library(stringr)
library(reshape2)
library(patchwork)

source(file.path("functions", "Plotting_HelperFunctions.R"))
source(file.path("functions", "General_HelperFunctions.R"))
source(file.path("functions", "Error_HelperFunctions.R"))
source(file.path("functions", "FileIO_HelperFunctions.R"))

options(scipen = 999)

datasets <- c("cain", "lau", "leng", "mathys", "morabito")#, "seaRef") #, "seaAD")

est_fields = list("deconRNASeq" = "out.all",
                  "dtangle" = "estimates",
                  "hspe" = "estimates",
                  "music" = "Est.pctRNA.weighted"
)

algorithms <- c(names(est_fields), "random")

granularity <- c("broad")

bulk_datasets <- c("Mayo", "MSBB", "ROSMAP")

best_params_all <- Get_AllBestParamsAsDf(reference_datasets, granularity)

errs_all <- Get_AllErrorsAsDf(bulk_datasets, reference_datasets, algorithms,
                              granularity, best_params_all$params)

errs_all <- melt(errs_all) %>% dplyr::rename(error_type = "variable")

errs_total <- subset(errs_all, tissue == "All")

errs_tissue <- subset(errs_all, tissue != "All")
errs_tissue$tissue <- paste(errs_tissue$test_data_name, errs_tissue$tissue)


##### Color setup #####

refs <- unique(errs_all$reference_data_name)
algs <- unique(errs_all$method)
tiss <- colSums(table(errs_tissue$tissue, errs_tissue$test_data_name) > 1)

# Reference data sets
reference_colors <- RColorBrewer::brewer.pal(length(refs), "Set3")
names(reference_colors) <- sort(refs)

# Algorithms
algorithm_colors <- RColorBrewer::brewer.pal(length(algs), "Set2")
names(algorithm_colors) <- sort(algs)

# Tissues (modified viridis turbo color scheme)
tiss <- colSums(table(errs_tissue$tissue, errs_tissue$test_data_name) > 1)
tissue_colors <- c(viridis::turbo(tiss[["Mayo"]], begin = 0.1, end = 0.2),
                   viridis::turbo(tiss[["MSBB"]], begin = 0.4, end = 0.55),
                   viridis::turbo(tiss[["ROSMAP"]], begin = 0.7, end = 0.9))
names(tissue_colors) <- sort(unique(errs_tissue$tissue))

# Slightly more pastel than default looks better
tissue_fill_colors <- str_replace(tissue_colors, "FF$", "88")
names(tissue_fill_colors) <- names(tissue_colors)

# Bulk data sets -- match some of the tissue colors
bulk_colors <- tissue_fill_colors[c(1, 4, 9)]
names(bulk_colors) <- sort(bulk_datasets)

##### Comparison across data sets

# TODO run an ANOVA to check if errors are significantly different between
# data sets and against "random"

# Correlation errors can all be plotted with the same axis constraints
ros_names <- unique(grep("_[props|pct]", errs_total$error_type, value = TRUE))
cor_names <- setdiff(unique(grep("cor", errs_total$error_type, value = TRUE)),
                     ros_names)
other_errs <- setdiff(errs_total$error_type, c(cor_names, ros_names))

pdf(file.path(dir_figures, "error_plots_ROSMAP_broad.pdf"), width=10, height = 12)

plt1 <- Plot_ErrsByMethod(errs_total, cor_names, fill = "test_data_name", fill_colors = bulk_colors)
plt2 <- Plot_ErrsByDataset(errs_total, cor_names, fill = "test_data_name", fill_colors = bulk_colors)
plt3 <- Plot_ErrsByMethod(errs_total, cor_names, x_axis = "test_data_name",
                          fill = "reference_data_name", fill_colors = reference_colors)
plt4 <- Plot_ErrsByDataset(errs_total, cor_names, x_axis = "test_data_name",
                           fill = "method", fill_colors = algorithm_colors)

print(plt1 / plt2 / plt3 / plt4 + plot_annotation(title = "Correlation by algorithm and dataset"))

# By tissue

plt5 <- Plot_ErrsByMethod(errs_tissue, cor_names, color = "tissue",
                          fill = "tissue", fill_colors = tissue_fill_colors) +
          scale_color_manual(values = tissue_colors)
plt6 <- Plot_ErrsByDataset(errs_tissue, cor_names, color = "tissue",
                           fill = "tissue", fill_colors = tissue_fill_colors) +
          scale_color_manual(values = tissue_colors)

print(plt5 / plt6 + plot_annotation(title = "Correlation by algorithm and dataset (by tissue, combined view)"))

plt7 <- Plot_ErrsByMethod(errs_tissue, cor_names,
                          facet_var = c("method", "test_data_name"),
                          color = "tissue", fill = "reference_data_name",
                          fill_colors = reference_colors) +
          scale_color_manual(values = tissue_fill_colors)

print(plt7 + plot_annotation(title = "Correlation by algorithm (by tissue, split view)"))

plt8 <- Plot_ErrsByDataset(errs_tissue, cor_names,
                           facet_var = c("reference_data_name", "test_data_name"),
                           color = "tissue", fill = "method",
                           fill_colors = algorithm_colors) +
          scale_color_manual(values = tissue_fill_colors)

print(plt8 + plot_annotation(title = "Correlation by dataset (by tissue, split view)"))


# Each error type needs a different fixed axis, but we plot all on the same
# page as separate rows
plts <- list()
for (B in other_errs) {
  plt <- Plot_ErrsByMethod(errs_total, B, facet_var = c("error_type", "method"),
                           fill = "test_data_name", fill_colors = bulk_colors)
  plts[[B]] <- plt
}

print(wrap_plots(plts, nrow = length(other_errs)) + plot_annotation(title = 'Errors by algorithm'))

plts <- list()
for (B in other_errs) {
  plt <- Plot_ErrsByDataset(errs_total, B, facet_var = c("error_type", "reference_data_name"),
                            fill = "test_data_name", fill_colors = bulk_colors)
  plts[[B]] <- plt
}

print(wrap_plots(plts, nrow = length(other_errs)) + plot_annotation(title = 'Errors by dataset'))

# By tissue

for (B in other_errs) {
  plt <- Plot_ErrsByMethod(errs_tissue, B,
                           facet_var = c("method", "test_data_name"),
                           color = "tissue", fill = "reference_data_name",
                           fill_colors = reference_colors) +
            scale_color_manual(values = tissue_fill_colors)
  print(plt + plot_annotation(title = paste0(B, ' by algorithm (by tissue)')))
}

for (B in other_errs) {
  plt <- Plot_ErrsByDataset(errs_tissue, B,
                            facet_var = c("reference_data_name", "test_data_name"),
                            color = "tissue", fill = "method",
                            fill_colors = algorithm_colors) +
            scale_color_manual(values = tissue_fill_colors)
  print(plt + plot_annotation(title = paste0(B, ' by dataset (by tissue)')))
}

# ROSMAP errors

errs_ros <- subset(errs_total, test_data_name == "ROSMAP")

# Difference between percent RNA and proportions is negligible right now, only
# graph one of them
ros_names2 <- grep("_pct", ros_names, value = TRUE)

plts <- list()
for (B in sort(ros_names2)) {
  plt <- Plot_ErrsByMethod(errs_ros, B, facet_var = c("error_type", "method"),
                           fill = "reference_data_name",
                           fill_colors = reference_colors)
  plts[[B]] <- plt
}

print(wrap_plots(plts, nrow = length(ros_names2)) + plot_annotation(title = 'ROSMAP Errors vs ROSMAP IHC'))

errs_ros_tissue <- subset(errs_tissue, test_data_name == "ROSMAP")
errs_ros_tissue$tissue <- str_replace(errs_ros_tissue$tissue, "ROSMAP ", "")

plts <- list()
for (B in sort(ros_names2)) {
  plt <- Plot_ErrsByMethod(errs_ros_tissue, B, x_axis = "tissue",
                           facet_var = c("error_type", "method"),
                           fill = "reference_data_name", fill_colors = reference_colors)
  plts[[B]] <- plt
}

print(wrap_plots(plts, nrow = length(ros_names2)) + plot_annotation(title = 'ROSMAP Errors vs ROSMAP IHC (by tissue)'))

dev.off()

#### Individual dataset examination #####

# TODO hspe is wrong
for (dataset in datasets) {
  pdf(file.path(dir_figures, paste0("error_plots_", dataset, "_broad_detailed.pdf")), width=10, height = 12)

  for (bulk_dataset in bulk_datasets) {
    bulk_se <- Load_BulkData(bulk_dataset)
    metadata <- as.data.frame(colData(bulk_se)) %>% dplyr::rename(subject = specimenID)

    ests_alg <- list()

    for (algorithm in algorithms) {
      if (algorithm == "random") {
        next
      }
      est_field <- est_fields[[algorithm]]

      ests <- Load_AlgorithmOutputList(algorithm, dataset, bulk_dataset, granularity)
      ests <- ests[names(ests) %in% errs_all$name]

      if (length(ests) == 0) {
        next
      }

      params <- lapply(ests, function(X) {
        as.list(X[["params"]] %>% select(-test_data_name))
      })

      param_ids <- sapply(params, paste, collapse = " ")

      ests_melt <- lapply(ests, "[[", est_field)

      #if (bulk_dataset == "ROSMAP") {
      #  ests_melt <- lapply(ests_melt, ConvertToROSMAPCelltypes, remove_unused = FALSE)
      #}

      ests_melt <- lapply(names(ests_melt), FUN = function(X) {
        ests_melt[[X]] <- as.data.frame(ests_melt[[X]])
        ests_melt[[X]]$name <- X
        ests_melt[[X]]$subject <- rownames(ests_melt[[X]])
        ests_melt[[X]]$param_id <- param_ids[X]
        ests_melt[[X]]
      })
      ests_melt <- do.call(rbind, ests_melt)
      ests_melt <- melt(ests_melt) %>% dplyr::rename(celltype = variable,
                                                     prop_est = value)

      #ests_melt$name <- str_replace(ests_melt$name, "_ROSMAP", "")
      ests_melt$algorithm <- algorithm

      ests_alg[[algorithm]] <- ests_melt
    }

    ests_alg <- do.call(rbind, ests_alg)

    best_params <- subset(best_params_all, reference_data_name == dataset &
                            algorithm %in% unique(ests_alg$algorithm) &
                            grepl(bulk_dataset, test_datasets))
    best_params <- best_params %>% group_by(algorithm) %>%
      mutate(title = paste("Best", metrics, "/", test_datasets),
             title_short = paste(algorithm, "params", as.numeric(factor(title)))) %>%
      ungroup()
    best_params$param_id <- sapply(best_params$params, paste, collapse = " ")

    best_ests <- ests_alg %>% merge(best_params, by = c("algorithm", "param_id"))

    # TODO this needs some fixing
    if (bulk_dataset == "ROSMAP") {
      ihc_props <- as.matrix(read.csv(file_rosmap_ihc_proportions, row.names = 1))
      A <- Load_AvgLibSize(dataset, granularity)
      A2 <- A

      # Some datasets are missing some vascular types
      for (col in c("Endo", "Peri", "VLMC")) {
        if (!(col %in% names(A2))) {
          A2 <- c(A2, 0)
          names(A2)[length(A2)] <- col
        }
      }

      # TODO this isn't quite right and the A matrix needs to be re-processed
      A2["Neuro"] <- mean(A2[c("Exc", "Inh")])
      A2["Oligo"] <- mean(A2[c("Oligo","OPC")])
      A2["Endo"] <- max(A2[c("Endo", "Peri", "VLMC")])

      A2 <- A2[colnames(ihc_props)]
      A2 <- A2 / sum(A2)
      ihc_pct <- ConvertPropCellsToPctRNA(ihc_props, A2)

      ihc_props <- melt(ihc_props) %>% dplyr::rename(subject = Var1, celltype = Var2, prop_est = value)
      ihc_props$title <- "IHC cell proportions"

      ihc_pct <- melt(ihc_pct) %>% dplyr::rename(subject = Var1, celltype = Var2, prop_est = value)
      ihc_pct$title <- "IHC RNA percent"

      tmp <- rbind(ihc_props, ihc_pct)
      tmp$algorithm <- "IHC"
      tmp$param_id <- 0
      tmp$name <- 0
      tmp$params <- 0
      tmp$metrics <- 0
      tmp$test_datasets <- bulk_dataset
      tmp$total <- 1
      tmp$reference_data_name <- dataset
      tmp$title_short <- tmp$title
      tmp <- tmp[,colnames(best_ests)]

      best_ests <- rbind(best_ests, tmp)
    }

    tmp <- best_ests %>% select(name, title) %>% distinct()
    titles <- tmp$title
    names(titles) <- tmp$name

    plot_id <- paste("Reference dataset:", dataset, " / Bulk dataset:", bulk_dataset)

    plt <- ggplot(best_ests, aes(x = algorithm, y = prop_est, color = title)) +
      geom_boxplot(width = 0.5) + theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      theme(legend.position = "bottom") +
      facet_wrap(~celltype, scales = "fixed",
                 nrow = 3)
    #ncol = ceiling(length(unique(best_ests$celltype)) / 2))

    print(plt + plot_annotation(title = plot_id, subtitle = "Estimated proportions across all subjects (each bar is a single parameter set)"))

    # TODO make functions to do this and move these files to a better place
    ests_ad <- merge(best_ests, metadata, by = "subject") %>%
                subset(diagnosis %in% c("CT", "AD")) %>%
                subset(celltype %in% levels(ests_alg$celltype)) # Gets rid of added cell types from ROSMAP IHC
    ests_ad$celltype <- factor(ests_ad$celltype)

    significant <- list()
    for (param_set in unique(ests_ad$title_short)) {
      ests_param <- subset(ests_ad, title_short == param_set)
      anov <- aov(prop_est ~ diagnosis*celltype*tissue, data = ests_param)
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
      ests_params <- subset(ests_ad, algorithm == alg)
      plt <- ggplot(ests_params, aes(x = celltype, y = prop_est, fill = diagnosis, color = significant), group = tissue) +
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


