library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)

source(file.path("functions", "FileIO_HelperFunctions.R"))

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito")#,
              #"seaRef") #, "seaAD")

datatypes <- c("donors")#, "training")

algorithm <- "dtangle"

for (datatype in datatypes) {
  params_list <- list()

  for (dataset in datasets) {

    err_list <- Load_ErrorList("dtangle", dataset, datatype, "broad")

    if (length(err_list) == 0) {
      next
    }

    params_df <- do.call(rbind, err_list[["params"]])

    err_list <- err_list[["means"]]

    params_df <- cbind(params_df, err_list)
    params_list[[dataset]] <- params_df
  }
  bests <- do.call(rbind, params_list)

  if (algorithm == "dtangle") {
    bests <- bests %>%
              mutate(marker = paste(marker_type, marker_subtype, n_markers)) %>%
              select(-marker_type, -marker_subtype, -n_markers)
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
      print(table(tmp$gamma_name))# / table(bests$gamma_name)[names(table(tmp$gamma_name))])
      print(table(tmp$sum_fn_type))# / table(bests$sum_fn_type)[names(table(tmp$sum_fn_type))])
      print(table(tmp$marker))# / table(bests$marker)[names(table(tmp$marker))])
      print(table(tmp$input_type))# / table(bests$input_type)[names(table(tmp$input_type))])
    }
  }

  else if (algorithm = "deconRNASeq") {
    for (tmp in list(tmp5, tmp6, tmp7, tmp8)) {
      print(table(tmp$use_scale) / table(bests$use_scale)[names(table(tmp$use_scale))])
      print(table(tmp$marker) / table(bests$marker)[names(table(tmp$marker))])
    }
  }
}
