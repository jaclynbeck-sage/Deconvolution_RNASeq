library(dplyr)
library(stringr)

source(file.path("functions", "FileIO_HelperFunctions.R"))
source(file.path("functions", "Error_HelperFunctions.R"))

datasets <- c("cain", "lau", "leng", "mathys", "seaRef") #, "seaAD")

granularity <- "broad_class"
bulk_datasets <- c("Mayo", "MSBB", "ROSMAP")
algorithms <- c("DeconRNASeq", "Dtangle", "DWLS", "Music", "HSPE", "Random")

for (bulk_dataset in bulk_datasets) {
  for (algorithm in algorithms) {
    errs_alg <- lapply(datasets, function(dataset) {
      err_files <- Get_ErrorFiles(algorithm, dataset, bulk_dataset, granularity)

      if (length(err_files) == 0) {
        message(str_glue("No data found for {bulk_dataset}/{dataset}/{algorithm}. Skipping..."))
        return(NULL)
      }

      errs_dataset <- lapply(err_files, function(EF) {
        err_list <- readRDS(EF)

        if (length(err_list) == 0) {
          next
        }

        errs_all <- err_list$means$all_tissue
        rownames(err_list$params) <- err_list$means$all_tissue$param_id
        names(err_list$by_sample) <- err_list$means$all_tissue$param_id
        #errs_tissue <- data.frame(err_list$means$by_tissue)

        #pars <- err_list$params %>% select(-reference_data_name, -test_data_name)

        get_best_vals <- function(col_name, errs_df) {
          if (length(grep("cor", col_name)) > 0) {
            top_ind <- which.max(errs_df[,col_name])
          }
          else {
            top_ind <- which.min(errs_df[,col_name])
          }
          return(top_ind)
        }

        cols_test <- setdiff(colnames(errs_all), "param_id")
        best_inds <- sapply(cols_test, get_best_vals, errs_all)
        #best_inds_tissue <- sapply(unique(errs_tissue$tissue_assignments), function(tissue) {
        #  df <- sapply(cols_test, get_best_vals,
        #               subset(errs_tissue, tissue_assignments == tissue))
        #  return(data.frame(t(df), tissue = tissue))
        #})
        #best_inds_tissue <- t(best_inds_tissue)

        errs_sub <- data.frame(param_id = errs_all$param_id[best_inds],
                               best_group = names(best_inds)) %>%
                      group_by(param_id) %>%
                      summarize(metrics = str_c(best_group, collapse = ", "))

        err_list_new <- list(means_all_tissue = subset(err_list$means$all_tissue,
                                                              param_id %in% errs_sub$param_id),
                             means_by_tissue = subset(err_list$means$by_tissue,
                                                      param_id %in% errs_sub$param_id),
                             params = err_list$params[errs_sub$param_id,],
                             by_sample = err_list$by_sample[errs_sub$param_id],
                             metrics = errs_sub)
        # TODO tissue-specific
        return(err_list_new)
      })

      return(Flatten_ErrorList(errs_dataset))
    })

    best_errors <- Flatten_ErrorList(errs_alg)

    if (length(best_errors) == 0) {
      message(str_glue("No data found for {bulk_dataset}/{algorithm}. Skipping..."))
      next
    }

    # Save some information about how many times each parameter set came up as
    # the best, ignoring dataset and normalization/regression
    param_strings <- format(best_errors$params, trim = TRUE) %>%
                        select(-reference_data_name, -test_data_name,
                               -normalization, -regression_method,
                               -total_markers_used) %>%
                        apply(1, paste, collapse = "_")

    param_strings <- data.frame(param_id = names(param_strings),
                                param_string = as.character(param_strings))
    param_strings <- merge(param_strings, best_errors$metrics, by = "param_id")

    param_stats <- param_strings %>% group_by(param_string) %>%
                      summarize(metrics = str_c(unique(metrics), collapse = ", "),
                                param_ids = str_c(unique(param_id), collapse = ", "),
                                count = n())
    best_errors$param_stats <- param_stats
    best_errors$params$algorithm <- algorithm

    saveRDS(best_errors,
            file.path(dir_best_errors,
                      str_glue("best_errors_{bulk_dataset}_{algorithm}_{granularity}.rds")))
  }
}
