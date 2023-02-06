library(dplyr)
library(tidyr)
library(stringr)
source("Filenames.R")

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito",
              "seaRef") #, "seaAD")

datatypes = c("donors", "training")

for (dataset in datasets) {
  params <- list()

  for (datatype in datatypes) {
    err_file <- file.path(dir_output,
                          paste0("errors_", dataset, "_", datatype, "_broad_shortsig.rds"))
    if (!file.exists(err_file)) {
      break
    }

    err_list <- readRDS(err_file)

    for (algorithm in names(err_list)) {
      errs <- err_list[[algorithm]][["means"]]

      bests <- lapply(colnames(errs), FUN = function(X) {
        if (length(grep("cor", X)) > 0) {
          ord <- order(errs[, X], decreasing = TRUE)
        }
        else {
          ord <- order(errs[, X], decreasing = FALSE)
        }

        ord <- ord[1:min(3, length(ord))]
        data.frame(name = rownames(errs[ord,]), group = X, rank = 1:length(ord))
      })

      bests <- do.call(rbind, bests)

      # For now, don't use goodness of fit
      bests <- subset(bests, !grepl("gof", group))

      # Get rank 1 param for each error metric
      top_one <- bests %>% group_by(group) %>% top_n(1, wt = -rank)

      top_one$algorithm <- algorithm
      top_one$datatype <- datatype
      params <- append(params, list(top_one))
    }
  }

  params <- do.call(rbind, params)

  # For duplicate names, create a new column with a list of each metric
  # and a list of each data type for the duplicates. Also count how many times
  # a particular parameter set was duplicated. Remove the old group,
  # datatype, and rank columns
  params <- params %>% group_by(name) %>%
              mutate(metrics = str_c(unique(group), collapse = ", "),
                     datatypes = str_c(unique(datatype), collapse = ", "),
                     total = n(), .keep = "unused") %>%
              select(-rank) %>% distinct()

  saveRDS(params, file = file.path(dir_output, paste0("best_params_", dataset,
                                                      "_broad.rds")))
}
