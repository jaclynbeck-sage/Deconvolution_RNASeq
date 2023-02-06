library(Matrix)
library(ggplot2)
#library(viridis)
library(dplyr)
library(stringr)
#library(SingleCellExperiment)
library(reshape2)
library(patchwork)

source("Filenames.R")
source(file.path("functions", "Plotting_HelperFunctions.R"))
source(file.path("functions", "General_HelperFunctions.R"))

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito",
              "seaRef") #, "seaAD")

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
errs_all <- lapply(datasets, function(dataset) {
  errs <- readRDS(file.path(dir_output, paste0("errors_", dataset, "_broad_ROSMAP.rds")))
  errs <- lapply(names(errs), FUN = function(X) {
    errs[[X]]$method <- X
    errs[[X]]$name <- rownames(errs[[X]])
    return(errs[[X]])
  })
  errs <- melt(do.call(rbind, errs))
  errs$dataset <- dataset
  return(errs)
})

errs_all <- do.call(rbind, errs_all) %>% rename(error.type = variable)

##### Comparison across data sets

# Correlation errors can all be plotted with the same axis constraints
cor_names <- unique(grep("cor", errs_all$error.type, value = TRUE))

pdf(file.path(dir_figures, "error_plots_ROSMAP.pdf"), width=8, height = 8)

plt1 <- Plot_ErrsByMethod(errs_all, cor_names, c(), "dataset")
plt2 <- Plot_ErrsByDataset(errs_all, cor_names, c(), "method")

print(plt1 / plt2 + plot_annotation(title = "Correlation by algorithm and dataset"))

# Each error type needs a different fixed axis, but we plot all 3 on the same
# page as separate rows
other_errs <- setdiff(errs_all$error.type, cor_names)
plts <- list()
for (B in other_errs) {
  plt <- Plot_ErrsByMethod(errs_all, B, c(), "dataset")
  plts[[B]] <- plt
}

print(wrap_plots(plts, nrow = 3) + plot_annotation(title = 'Errors by algorithm'))

plts <- list()
for (B in other_errs) {
  plt <- Plot_ErrsByDataset(errs_all, B, c(), "method")
  plts[[B]] <- plt
}

print(wrap_plots(plts, nrow = 3) + plot_annotation(title = 'Errors by dataset'))

#### Individual dataset examination #####

for (dataset in datasets) {
  #pdf(file.path(dir_figures, paste0("error_plots_ground_truth_", dataset, "_detailed.pdf")), width=8, height = 8)
  dat <- dataset # Issues with the column name and condition being the same
  errs_dataset <- subset(errs_all, dataset == dat)

  ests_alg <- list()
  for (algorithm in algorithms) {
    est_field <- est_fields[[algorithm]]
    alg_name <- algorithm
    if (alg_name == "music_wt" | alg_name == "music_nnls") {
      alg_name <- "music"
    }
    alg_file <- file.path(dir_rosmap,
                          paste0(alg_name, "_list_", dataset, "_broad_ROSMAP.rds"))

    if (!file.exists(alg_file)) {
      next
    }
    ests <- readRDS(file = alg_file)

    if (convert[[algorithm]] == TRUE) {
      A <- readRDS(file.path(dir_input, paste0(dataset, "_A_matrix.rds")))
      for (i in 1:length(ests)) {
        ests[[i]][[est_field]] <- ConvertPropCellsToPctRNA(ests[[i]][[est_field]], A[["A_broad"]])
      }
    }

    ests_melt <- lapply(ests, "[[", est_field)
    ests_melt <- lapply(names(ests_melt), FUN = function(X) {
      ests_melt[[X]] <- as.data.frame(ests_melt[[X]])
      ests_melt[[X]]$name <- X
      ests_melt[[X]]$subject <- rownames(ests_melt[[X]])
      ests_melt[[X]]
    })
    ests_melt <- do.call(rbind, ests_melt)
    ests_melt <- melt(ests_melt)
    colnames(ests_melt) <- c("name", "subject", "celltype", "prop_est")
    ests_melt$name <- str_replace(ests_melt$name, "_ROSMAP", "")
    ests_melt$algorithm <- algorithm

    ests_alg[[algorithm]] <- ests_melt
  }

  ests_alg <- do.call(rbind, ests_alg)

  best_params <- subset(best_params_all, dataset == dat & algorithm %in% unique(ests_alg$algorithm))
  best_params <- best_params %>% group_by(algorithm) %>%
                    mutate(title = paste("Best", metrics, "/", datatypes),
                           title_short = paste(algorithm, "params", as.numeric(factor(title)))) %>%
                    ungroup()

  best_ests <- ests_alg %>% merge(best_params, by = c("name", "algorithm"))

  tmp <- best_ests %>% select(name, title) %>% distinct()
  titles <- tmp$title
  names(titles) <- tmp$name

  plot_id <- paste("Reference dataset:", dataset)

  plt <- ggplot(best_ests, aes(x = algorithm, y = prop_est, color = algorithm, fill = title_short)) +
    geom_boxplot(width = 0.5) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    facet_wrap(~celltype + dataset, scales = "free",
               nrow = 3)
               #ncol = ceiling(length(unique(best_ests$celltype)) / 2))

  print(plt + plot_annotation(title = plot_id, subtitle = "Estimated proportions across all subjects (each bar is a single parameter set)"))

  metadata <- read.table(file.path(dir_input, "ageCensoredCovariates.tsv"), header = TRUE, sep = "\t") # TODO use the canonical files from ROSMAP
  metadata <- metadata %>% select(SampleID, c(Diagnosis, Diagnosis.msex)) %>%
                rename(subject = SampleID)
  metadata$subject <- make.names(metadata$subject)

  ests_ad <- merge(best_ests, metadata, by = "subject") %>% subset(Diagnosis %in% c("CONTROL", "AD"))

  significant <- list()
  for (param_set in unique(ests_ad$title_short)) {
    ests_param <- subset(ests_ad, title_short == param_set)
    anov <- aov(prop_est ~ Diagnosis*celltype, data = ests_param)
    tuk <- TukeyHSD(anov, "Diagnosis:celltype")

    comparisons <- paste0("CONTROL:", levels(ests_ad$celltype), "-AD:", levels(ests_ad$celltype))
    tuk <- as.data.frame(tuk[[1]][comparisons,])
    tuk$p_adj <- tuk$`p adj`
    tuk$celltype <- levels(ests_ad$celltype)
    tuk$significant <- tuk$p_adj <= 0.05
    tuk$title_short <- param_set
    significant[[param_set]] <- tuk
  }

  significant <- do.call(rbind, significant) %>% select(p_adj, celltype, significant, title_short)

  ests_ad <- merge(ests_ad, significant, by = c("celltype", "title_short"))

  for (alg in unique(ests_ad$algorithm)) {
    ests_params <- subset(ests_ad, algorithm == alg)
    plt <- ggplot(ests_params, aes(x = celltype, y = prop_est, fill = Diagnosis, color = significant)) +
      geom_boxplot(width = 0.5) + theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_color_manual(values = c("#dddddd", "#000000")) +
      facet_wrap(~title_short, scales = "free",
                 ncol = 1)
    print(plt + plot_annotation(title = plot_id, subtitle = paste(alg, "estimates, AD vs Control")))
  }
  # TODO can we break down goodness-of-fit errors by cell type and subject?
}

dev.off()
