library(Matrix)
library(ggplot2)
library(viridis)
library(dplyr)
library(stringr)
library(SummarizedExperiment)
library(reshape2)
library(patchwork)

source("Filenames.R")
source(file.path("functions", "Plotting_HelperFunctions.R"))
source(file.path("functions", "General_HelperFunctions.R"))

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito",
              "seaRef") #, "seaAD")

datatypes <- c("donors", "training")

est_fields = list("dtangle" = "estimates",
                  "music_wt" = "Est.prop.weighted",
                  "music_nnls" = "Est.prop.allgene",
                  #"music2" = "Est.prop",
                  "hspe" = "estimates",
                  "deconRNASeq" = "Est.prop")

convert <- list("dtangle" = FALSE,
                "music_wt" = TRUE,
                "music_nnls" = TRUE,
                "music2" = TRUE,
                "hspe" = FALSE,
                "deconRNASeq" = FALSE)

algorithms <- names(est_fields)

best_params_all <- lapply(datasets, function(dataset) {
  readRDS(file.path(dir_output, paste0("best_params_", dataset,
                                       "_broad.rds")))
})
best_params_all <- do.call(rbind, best_params_all) %>% mutate(dataset = str_replace(name, "_.*", ""))

# Read in all errors for each dataset & data type,
# put in one big dataframe
errs_means <- Extract_MeanErrors(datasets, datatypes, best_params_all, dir_output)

# ignore goodness of fit for now
errs_means <- subset(errs_means, !grepl("gof", error.type))

##### Comparison across data sets

# Correlation errors can all be plotted with the same axis constraints
cor_names <- unique(grep("cor_", errs_means$error.type, value = TRUE))

pdf(file.path(dir_figures, "error_plots_ground_truth_summary.pdf"), width=8, height = 8)

plt <- Plot_ErrsByMethod(errs_means, cor_names, "donors", "dataset")
print(plt + plot_annotation(title = 'Correlation by algorithm (donors)'))

plt <- Plot_ErrsByMethod(errs_means, cor_names, "training", "dataset")
print(plt + plot_annotation(title = 'Correlation by algorithm (training set)'))

plt <- Plot_ErrsByDataset(errs_means, cor_names, "donors", "method")
print(plt + plot_annotation(title = 'Correlation by dataset (donors)'))

plt <- Plot_ErrsByDataset(errs_means, cor_names, "training", "method")
print(plt + plot_annotation(title = 'Correlation by dataset (training set)'))

# Each error type needs a different fixed axis, but we plot all 3 on the same
# page as separate rows
other_errs <- setdiff(errs_means$error.type, cor_names)
plts <- list()
for (B in other_errs) {
  plt <- Plot_ErrsByMethod(errs_means, B, "donors", "dataset")
  plts[[B]] <- plt
}

print(wrap_plots(plts, nrow = 3) + plot_annotation(title = 'Errors by algorithm (donors)'))

plts <- list()
for (B in other_errs) {
  plt <- Plot_ErrsByMethod(errs_means, B, "training", "dataset")
  plts[[B]] <- plt
}

print(wrap_plots(plts, nrow = 3) + plot_annotation(title = 'Errors by algorithm (training set)'))


plts <- list()
for (B in other_errs) {
  plt <- Plot_ErrsByDataset(errs_means, B, "donors", "method")
  plts[[B]] <- plt
}

print(wrap_plots(plts, nrow = 3) + plot_annotation(title = 'Errors by dataset (donors)'))

plts <- list()
for (B in unique(other_errs)) {
  plt <- Plot_ErrsByDataset(errs_means, B, "training", "method")
  plts[[B]] <- plt
}

print(wrap_plots(plts, nrow = 3) + plot_annotation(title = 'Errors by dataset (training set)'))


plt <- Plot_ErrsByMethod(errs_means, cor_names, datatypes, "datatype")
print(plt + plot_annotation(title = 'Correlation by algorithm (donors vs training)'))

plt <- Plot_ErrsByDataset(errs_means, cor_names, datatypes, "datatype")
print(plt + plot_annotation(title = 'Correlation by dataset (donors vs training)'))


plts <- list()
for (B in unique(other_errs)) {
  plt <- Plot_ErrsByMethod(errs_means, B, datatypes, "datatype")
  plts[[B]] <- plt
}

print(wrap_plots(plts, nrow = 3) + plot_annotation(title = 'Errors by algorithm (donors vs training)'))

plts <- list()
for (B in unique(other_errs)) {
  plt <- Plot_ErrsByDataset(errs_means, B, datatypes, "datatype")
  plts[[B]] <- plt
}

print(wrap_plots(plts, nrow = 3) + plot_annotation(title = 'Errors by dataset (donors vs training)'))

dev.off()

#### Individual data set examination #####

for (dataset in datasets) {
  pdf(file.path(dir_figures, paste0("error_plots_ground_truth_", dataset, "_detailed.pdf")), width=8, height = 8)
  dat <- dataset # Issues with the column name and condition being the same
  errs_dataset <- subset(errs_means, dataset == dat)

  errs_cor <- subset(errs_dataset, error.type %in% cor_names)
  errs_other <- subset(errs_dataset, error.type %in% other_errs)

  plt1 <- ggplot(errs_cor, aes(x = method, y = value, color = method)) +
            geom_boxplot() + theme_bw() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            facet_wrap(~error.type, scales = "fixed")

  plt2 <- ggplot(errs_other, aes(x = method, y = value, color = method)) +
    geom_boxplot() + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~error.type, scales = "free")

  print(plt1 / plt2 + plot_annotation(title = paste0("Errors for ", dataset, ": best parameters")))

  for (datatype in datatypes) {
    se <- readRDS(file.path(dir_pseudobulk,
                            paste0("pseudobulk_", dataset, "_", datatype,
                                   "_broadcelltypes.rds")))
    pctRNA <- as.matrix(metadata(se)[["pctRNA"]])

    for (algorithm in algorithms) {
      plot_id <- paste(dataset, datatype, algorithm, sep = " / ")

      est_field <- est_fields[[algorithm]]

      dtype <- datatype # Issues with the column name and condition being the same
      errs_algorithm <- subset(errs_means, dataset == dat & datatype == dtype &
                                 method == algorithm)

      alg_name <- algorithm
      if (alg_name == "music_wt" | alg_name == "music_nnls") {
        alg_name <- "music"
      }
      alg_file <- file.path(dir_params_lists,
                            paste0(alg_name, "_list_", dataset, "_",
                                   datatype, "_broad.rds"))

      if (!file.exists(alg_file)) {
        next
      }
      ests <- readRDS(file = file.path(dir_params_lists,
                                       paste0(alg_name, "_list_", dataset, "_",
                                              datatype, "_broad.rds")))
      if (convert[[algorithm]] == TRUE) {
        A <- readRDS(file.path(dir_input, paste0(dataset, "_A_matrix.rds")))
        for (i in 1:length(ests)) {
          ests[[i]][[est_field]] <- ConvertPropCellsToPctRNA(ests[[i]][[est_field]], A[["A_broad"]])
        }
      }

      ests_df <- MakePropsDataframe(pctRNA, ests, est_field)

      alg <- algorithm # Issues with the column name and condition being the same
      best_params <- subset(best_params_all, dataset == dat & algorithm == alg)

      best_ests <- subset(ests_df, name %in% best_params$name) %>%
        merge(best_params, by = "name") %>%
        mutate(title = paste("Best", metrics, "/", datatypes))
      tmp <- best_ests %>% select(name, title) %>% distinct()
      titles <- tmp$title
      names(titles) <- tmp$name

      plt <- ggplot(best_ests, aes(x = prop_truth, y = prop_est, color = celltype)) +
        geom_jitter(size = 1) + geom_abline(slope = 1) +
        theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(aspect.ratio = 1) + lims(x = c(0,1), y = c(0,1)) +
        facet_wrap(~name, scales = "free", nrow = 2, labeller = labeller(name = titles))
      print(plt + plot_annotation(title = plot_id, subtitle = "Best Estimates vs Ground Truth"))

      # How well is each cell type characterized?

      errs <- readRDS(file = file.path(dir_output, paste0("errors_", dataset, "_",
                                                          datatype, "_broad_shortsig.rds")))

      errs_by_celltype <- Extract_ErrorsByCelltype(errs, algorithm, best_params_all)

      errs1 <- subset(errs_by_celltype, !grepl("gof", error.type))
      #errs2 <- subset(errs_by_celltype, grepl("gof", error.type))

      plt <- ggplot(errs1, aes(x = celltype, y = value, fill = celltype)) +
        geom_violin() + geom_jitter(size = 0.5, width = 0.1) +
        theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        facet_wrap(~error.type, scales = "free")
      print(plt + plot_annotation(title = plot_id, subtitle = "Best errors for individual cell types"))

      ##### Donors only
      if (datatype == "donors") {
        # How well was each subject characterized? -- not as useful
        errs_by_subject <- Extract_ErrorsBySubject(errs, algorithm, best_params_all)

        errs1 <- subset(errs_by_subject, !grepl("gof", error.type))
        plt <- ggplot(errs1, aes(x = subject, y = value, fill = subject)) +
          geom_violin() + geom_jitter(size = 0.5, width = 0.1) +
          theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          theme(legend.position = "none") +
          facet_wrap(~error.type, scales = "free")
        print(plt + plot_annotation(title = plot_id, subtitle = "Best errors for individual donors"))
      }

      # How well is each cell type characterized by proportion?

      # TODO it would be useful to highlight actual biological range from donors
      # for each cell type
      plt1 <- Plot_TruthVsEstimates_Dots(best_ests, titles, n_col = ceiling(ncol(pctRNA)/2))
      print(plt1 + plot_annotation(title = plot_id, subtitle = "Best Estimates vs Ground Truth by Cell Type"))

      # Lines instead of points
      ests_avg <- best_ests %>% mutate(prop_truth = round(prop_truth, digits = 2)) %>%
        group_by(name, prop_truth, celltype) %>%
        mutate(avg_est = mean(prop_est)) %>%
        select(-prop_est) %>% distinct()

      plt2 <- Plot_TruthVsEstimates_Lines(ests_avg, titles, n_col = ceiling(ncol(pctRNA)/2))
      print(plt2 + plot_annotation(title = plot_id, subtitle = "Best Estimates vs Ground Truth by Cell Type (Average)"))

      # Zoom in on the crowded area < 30%
      ests2 <- subset(best_ests, prop_est <= 0.3 | prop_truth <= 0.3)
      plt <- Plot_TruthVsEstimates_Dots(ests2, titles, n_col = ceiling(ncol(pctRNA)/2), axis_limits = c(0, 0.35))
      print(plt + plot_annotation(title = plot_id, subtitle = "Best Estimates vs Ground Truth by Cell Type, zoomed to 0-30% RNA"))

      ests3 <- subset(ests_avg, prop_truth <= 0.3 | avg_est <= 0.3)
      plt <- Plot_TruthVsEstimates_Lines(ests3, titles, n_col = ceiling(ncol(pctRNA)/2), axis_limits = c(0, 0.35))
      print(plt + plot_annotation(title = plot_id, subtitle = "Best Estimates vs Ground Truth by Cell Type (Average), zoomed to 0-30% RNA"))
    }
  }

  dev.off()
}


