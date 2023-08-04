library(dplyr)
library(stringr)

source(file.path("functions", "FileIO_HelperFunctions.R"))

datasets <- c("cain", "lau", "leng", "mathys", "morabito")#, "seaRef") #, "seaAD")

granularity <- "broad"
bulk_datasets <- c("ROSMAP", "Mayo", "MSBB")
algorithms = c("deconRNASeq", "dtangle", "music", "random")

for (dataset in datasets) {
  params <- list()

  for (bulk_dataset in bulk_datasets) {
    for (algorithm in algorithms) {

      err_list <- Load_ErrorList(algorithm, dataset, bulk_dataset, granularity)
      if (length(err_list) == 0) {
        next
      }

      errs <- err_list[["gof_mean"]]

      if (bulk_dataset == "ROSMAP" & "errs_ihc_props" %in% names(err_list)) {
        errs_ihc_props <- lapply(err_list[["errs_ihc_props"]], function(err) {
          err %>% summarize(across(colnames(err), mean))
        })
        errs_ihc_props <- do.call(rbind, errs_ihc_props) %>%
                            dplyr::rename(cor_props = cor,
                                          rMSE_props = rMSE,
                                          mAPE_props = mAPE)

        errs_ihc_pct <- lapply(err_list[["errs_ihc_pct"]], function(err) {
          err %>% summarize(across(colnames(err), mean))
        })
        errs_ihc_pct <- do.call(rbind, errs_ihc_pct) %>%
                            dplyr::rename(cor_pct = cor,
                                          rMSE_pct = rMSE,
                                          mAPE_pct = mAPE)

        errs <- cbind(errs, errs_ihc_props, errs_ihc_pct)
      }

      pars <- do.call(rbind, err_list[["params"]]) %>%
                select(-test_data_name)

      get_best_vals <- function(col_name, errs_df) {
        if (length(grep("cor", col_name)) > 0) {
          top_ind <- which.max(errs_df[,col_name])
        }
        else {
          top_ind <- which.min(errs_df[,col_name])
        }
        return(top_ind)
      }

      bests <- sapply(colnames(errs), get_best_vals, errs)
      bests <- data.frame(name = rownames(errs)[bests],
                          group = names(bests))

      bests$params <- lapply(bests$name, function(X) { as.list(pars[X,]) })

      bests$algorithm <- algorithm
      bests$test_data_name <- bulk_dataset
      params <- append(params, list(bests))
    }
  }

  params <- do.call(rbind, params)

  # For duplicate names, create a new column with a list of each metric
  # and a list of each data type for the duplicates. Also count how many times
  # a particular parameter set was duplicated. Remove the old group,
  # datatype, and rank columns
  params <- params %>% select(-name) %>% group_by(algorithm, params) %>%
              mutate(metrics = str_c(unique(group), collapse = ", "),
                     test_datasets = str_c(unique(test_data_name), collapse = ", "),
                     total = n(), .keep = "unused") %>% distinct()

  saveRDS(params, file = file.path(dir_output,
                                   str_glue("best_params_{dataset}_{granularity}.rds")))
}
