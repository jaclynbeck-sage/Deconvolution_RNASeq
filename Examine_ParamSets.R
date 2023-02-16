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
      break
    }

    err_list <- readRDS(err_file)
    err_list <- err_list[[algorithm]][["means"]]

    # TODO temporary
    if (algorithm == "dtangle") {
      rownames(err_list) <- str_replace(rownames(err_list), "broad_method", "broad_input_singlecell_method")
    }

    if (algorithm == "dtangle") {
      params_df <- extract_dtangle_params(rownames(err_list))
    }
    else if (algorithm == "deconRNASeq") {
      params_df <- extract_deconRNASeq_params(rownames(err_list))
    }
    params_df <- cbind(params_df, err_list)
    params_df$dataset <- dataset
    params_list[[dataset]] <- params_df
  }
  params_df <- do.call(rbind, params_list)

  if ("gof.cor" %in% colnames(params_df)) {
    bests <- params_df %>% select(-gof.cor:-gof.mAPE)
  }
  else {
    bests <- params_df
  }

  if (algorithm == "dtangle") {
    bests <- bests %>%
              mutate(marker = paste(marker_meth, n_markers),
                     input_type = paste(input, normtype)) %>%
              select(-marker_meth, -n_markers, -input, -normtype)

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
      mutate(marker = paste(filterlvl, n_markers)) %>%
      select(-filterlvl, -n_markers)
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
      print(table(tmp$usescale) / table(bests$usescale)[names(table(tmp$usescale))])
      print(table(tmp$normtype) / table(bests$normtype)[names(table(tmp$normtype))])
      print(table(tmp$marker) / table(bests$marker)[names(table(tmp$marker))])
    }
  }

}
