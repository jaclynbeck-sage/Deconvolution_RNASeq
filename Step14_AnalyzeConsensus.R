library(stringr)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(Hmisc)

source(file.path("functions", "Step14_Analysis_HelperFunctions.R"))

algorithms <- c("CibersortX", "DeconRNASeq", "Dtangle", "DWLS", "HSPE", "Music",
                "Scaden", "Baseline")

granularity <- c("sub_class")

bulk_datasets <- c("Mayo", "MSBB", "ROSMAP")

combined_metadata <- lapply(bulk_datasets, function(B) {
  bulk <- Load_BulkData(B)
  meta <- as.data.frame(colData(bulk)) %>% select(sample, diagnosis, tissue)
})
combined_metadata <- do.call(rbind, combined_metadata)

# Get the full set of all best errors
best_errors <- Get_AllBestErrorsAsDf(bulk_datasets, granularity)

# We are analyzing each tissue separately, so we will use the best parameter
# sets for each tissue instead of the best overall
best_errors <- subset(best_errors, tissue != "All")

best_signatures <- Find_BestSignature(best_errors)

# Only keep the error information for scores from the top signature
best_errors <- merge(best_errors, best_signatures,
                     by = c("tissue", "normalization", "regression_method"))
best_errors <- subset(best_errors, signature == best_signature) %>%
  select(-best_signature)

group_cols <- c("tissue", "reference_data_name", "test_data_name",
                "normalization", "regression_method", "algorithm")
group_cols_toplevel <- setdiff(group_cols, "reference_data_name")

# Best parameters for each reference_data_name for each possible data input
# (bulk dataset, normalization, regression, input type, and algorithm)
best_params <- Find_BestParameters(best_errors, group_cols)

# Best parameters for each possible data input, ignoring reference_data_name
best_params_toplevel <- Find_BestParameters(best_errors, group_cols_toplevel)

# Only keep error information for the best parameters
best_errors <- merge(best_params, best_errors, by = c("tissue", "param_id"),
                     all.y = FALSE)

best_errors_toplevel <- merge(best_params_toplevel, best_errors,
                              by = c("tissue", "param_id"), all.y = FALSE)

saveRDS(list("best_errors_all" = best_errors,
             "best_errors_toplevel" = best_errors_toplevel),
        file.path(dir_analysis, paste0("best_errors_", granularity, ".rds")))

# Get the estimates associated with each parameter ID left in the errors df
best_estimates <- Get_AllBestEstimatesAsDf(bulk_datasets, granularity,
                                           combined_metadata, best_params)

# This is faster than doing a merge with best_params_toplevel w/ all.y = FALSE,
# because best_estimates is so large
best_estimates_toplevel <- lapply(unique(best_params_toplevel$tissue), function(tiss) {
  ests_tmp <- subset(best_estimates, tissue == tiss)
  params_tmp <- subset(best_params_toplevel, tissue == tiss)
  return(subset(ests_tmp, param_id %in% params_tmp$param_id))
})
best_estimates_toplevel <- do.call(rbind, best_estimates_toplevel)

gc()

best_errors$avg_id <- apply(best_errors, 1, function(row) {
  paste(row[group_cols], collapse = "_")
})

best_errors_toplevel$avg_id <- apply(best_errors_toplevel, 1, function(row) {
  paste(row[group_cols_toplevel], collapse = "_")
})

# Average the estimates corresponding to best correlation, best rMSE, and best
# mAPE together for each data input type
avg_list <- Create_AveragesList(best_errors, best_estimates)
avg_list_toplevel <- Create_AveragesList(best_errors_toplevel, best_estimates_toplevel)

saveRDS(list("average_list_all" = avg_list,
             "average_list_toplevel" = avg_list_toplevel),
        file.path(dir_analysis, paste0("averages_lists_", granularity, ".rds")))

# Calculate significance of cell type differences on a tissue-by-tissue basis
mean_props_all <- Get_MeanProps_Significance(avg_list)
mean_props_toplevel <- Get_MeanProps_Significance(avg_list_toplevel)

saveRDS(list("significance_props_all" = mean_props_all,
             "significance_props_toplevel" = mean_props_toplevel),
        file.path(dir_analysis, paste0("significance_lists_", granularity, ".rds")))

# This is really slow so we do it last
best_errors_single <- best_errors

# Mimic the structure of avg_list as returned from Create_AverageList
best_errors_single$avg_id <- paste(best_errors_single$tissue,
                                   best_errors_single$param_id,
                                   sep = "_")
best_estimates_single <- best_estimates %>%
  dplyr::rename(percent_mean = percent) %>%
  mutate(avg_id = paste(tissue, param_id, sep = "_")) %>%
  select(-param_id)

avg_list_single <- lapply(unique(best_errors_single$tissue), function(tiss) {
  return(list("avg_errors" = subset(best_errors_single, tissue == tiss),
              "avg_estimates" = subset(best_estimates_single, tissue == tiss)))
})
names(avg_list_single) <- unique(best_errors_single$tissue)

# Significance of each individual parameter set
mean_props_single <- Get_MeanProps_Significance(avg_list_single)

saveRDS(list("significance_props_all" = mean_props_all,
             "significance_props_toplevel" = mean_props_toplevel,
             "significance_props_single" = mean_props_single),
        file.path(dir_analysis, paste0("significance_lists_", granularity, ".rds")))
