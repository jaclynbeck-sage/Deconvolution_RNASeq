library(dplyr)
library(tidyr)
library(stringr)

source(file.path("functions", "FileIO_HelperFunctions.R"))

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito",
              "seaRef") #, "seaAD")

granularity <- "broad"
datatypes <- c("donors", "training")

for (dataset in datasets) {
  params <- list()

  for (datatype in datatypes) {
    err_list <- Load_ErrorList(dataset, datatype, granularity)

    if (length(err_list) == 0) {
      next
    }

    for (algorithm in names(err_list)) {
      errs <- err_list[[algorithm]][["means"]]

      errs_gof <- err_list[[algorithm]][["gof"]]
      errs_gof <- lapply(names(errs_gof), function(X) {
        errs_gof[[X]]$filter_lvl <- rownames(errs_gof[[X]])
        errs_gof[[X]]$name <- X
        errs_gof[[X]]
      })
      errs_gof <- do.call(rbind, errs_gof)

      pars <- do.call(rbind, err_list[[algorithm]][["params"]]) %>%
                select(-datatype)

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
                          filter_lvl = "0",
                          group = names(bests))

      bests_gof <- sapply(grep("gof", colnames(errs_gof), value = TRUE),
                          get_best_vals, errs_gof)
      bests_gof <- data.frame(name = errs_gof$name[bests_gof],
                              filter_lvl = errs_gof$filter_lvl[bests_gof],
                              group = names(bests_gof))

      bests <- rbind(bests, bests_gof)
      bests$params <- apply(pars[bests$name,], 1, as.list)

      bests$algorithm <- algorithm
      bests$datatype <- datatype
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
                     datatypes = str_c(unique(datatype), collapse = ", "),
                     filter_lvls = str_c(unique(filter_lvl), collapse = ", "),
                     total = n(), .keep = "unused") %>% distinct()

  saveRDS(params, file = file.path(dir_output,
                                   str_glue("best_params_{dataset}_{granularity}.rds")))
}
