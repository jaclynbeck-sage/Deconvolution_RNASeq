library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)
source("Filenames.R")
source(file.path("functions", "Error_HelperFunctions.R"))

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito",
              "seaRef") #, "seaAD")

datatypes <- c("donors", "training")

algorithm <- "deconRNASeq"

for (datatype in datatypes) {
  params_list <- list()

  for (dataset in datasets) {
    err_file <- file.path(dir_output,
                          paste0("errors_", dataset, "_", datatype, "_broad.rds"))
    if (!file.exists(err_file)) {
      next
    }

    err_list <- readRDS(err_file)
    params_df <- do.call(rbind, err_list[[algorithm]][["params"]])

    err_list <- err_list[[algorithm]][["means"]]

    params_df <- cbind(params_df, err_list)
    params_list[[dataset]] <- params_df
  }
  bests <- do.call(rbind, params_list)

  if (algorithm == "dtangle") {
    bests <- bests %>%
              mutate(marker = paste(marker_method, n_markers)) %>%
              select(-marker_method, -n_markers)

    tmp <- bests %>% select(dataset, input_type, cor_celltype:mAPE) %>%
              distinct() %>% group_by(dataset, input_type) %>%
              mutate(across(cor_celltype:mAPE, mean)) %>% distinct()

    tmp2 <- bests %>% select(dataset, gamma_name, cor_celltype:mAPE) %>%
      distinct() %>% group_by(dataset, gamma_name) %>%
      mutate(across(cor_celltype:mAPE, mean)) %>% distinct()

    tmp3 <- bests %>% select(dataset, sum_fn_type, cor_celltype:mAPE) %>%
      distinct() %>% group_by(dataset, sum_fn_type) %>%
      mutate(across(cor_celltype:mAPE, mean)) %>% distinct()

    tmp4 <- bests %>% select(dataset, marker, cor_celltype:mAPE) %>%
      distinct() %>% group_by(dataset, marker) %>%
      mutate(across(cor_celltype:mAPE, mean)) %>% distinct()
  }
  else if (algorithm == "deconRNASeq") {
    bests <- bests %>%
      mutate(marker = paste(marker_type, filter_level, n_markers)) %>%
      select(-marker_type, -filter_level, -n_markers)
  }

  tmp5 <- bests %>% group_by(dataset) %>% top_n(n = 10, wt = cor_celltype)
  tmp6 <- bests %>% group_by(dataset) %>% top_n(n = 10, wt = cor_subject)
  tmp7 <- bests %>% group_by(dataset) %>% top_n(n = 10, wt = -rMSE)
  tmp8 <- bests %>% group_by(dataset) %>% top_n(n = 10, wt = -mAPE)

  if (algorithm == "dtangle") {
    for (tmp in list(tmp5, tmp6, tmp7, tmp8)) {
      print(table(tmp$gamma_name) / table(bests$gamma_name)[names(table(tmp$gamma_name))])
      print(table(tmp$sum_fn_type) / table(bests$sum_fn_type)[names(table(tmp$sum_fn_type))])
      print(table(tmp$marker) / table(bests$marker)[names(table(tmp$marker))])
      print(table(tmp$input_type) / table(bests$input_type)[names(table(tmp$input_type))])
    }
  }

  else if (algorithm = "deconRNASeq") {
    for (tmp in list(tmp5, tmp6, tmp7, tmp8)) {
      print(table(tmp$use_scale) / table(bests$use_scale)[names(table(tmp$use_scale))])
      print(table(tmp$marker) / table(bests$marker)[names(table(tmp$marker))])
    }
  }
}
