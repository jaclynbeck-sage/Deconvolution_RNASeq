library(stringr)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(Hmisc)

source(file.path("functions", "Step14_Analysis_HelperFunctions.R"))

n_cores <- parallel::detectCores() - 2

granularity <- c("broad_class")

bulk_datasets <- c("Mayo", "MSBB", "ROSMAP")

combined_metadata <- lapply(bulk_datasets, function(B) {
  bulk <- Load_BulkData(B)
  meta <- as.data.frame(colData(bulk)) %>% select(sample, diagnosis, tissue)
})
combined_metadata <- do.call(rbind, combined_metadata)

# Get the full set of all best errors and quality stats
best_errors <- Get_AllBestErrorsAsDf(bulk_datasets, granularity, n_cores)
qstats <- Get_AllQualityStatsAsDf(bulk_datasets, granularity, n_cores)

# We are analyzing each tissue separately, so we will use the best parameter
# sets for each tissue instead of the best overall
best_errors <- subset(best_errors, tissue != "All") %>%
  dplyr::mutate(data_transform = paste(normalization, "+", regression_method),
                data_transform = str_replace(data_transform, "counts", "cpm"),
                data_transform = str_replace(data_transform, "cpm", "counts/cpm"),
                data_transform = str_replace(data_transform, "log_", ""))

# Top-level view of the top 3 results for each error metric, regardless of
# normalization, regression, algorithm, etc. We do NOT include Baseline
# estimates here because guessing all 0's nearly always out-performs everything
# else for rMSE.
err_ranks <- best_errors %>%
  subset(algorithm != "Baseline") %>%
  Rank_Errors(group_cols = "tissue") %>%
  Get_TopRanked(n_top = 3)

# Select a single signature per tissue/data transform against which errors were
# calculated, so all errors for that tissue/transform are directly comparable
#best_signatures <- Find_BestSignature(best_errors)

# Only keep the error information for scores from the top signature
#best_errors <- merge(best_errors, best_signatures,
#                     by = c("tissue", "data_transform"))
#best_errors <- subset(best_errors, signature == best_signature) %>%
#  select(-best_signature)
best_baselines <- subset(best_errors, algorithm == "Baseline" &
                           reference_data_name != "zeros") %>%
  Rank_Errors(group_cols = c("tissue", "data_transform")) %>%
  Get_TopRanked(n_top = 1) %>%
  select(-cor_rank, -rMSE_rank, -mAPE_rank, -mean_rank, -type) %>%
  distinct()

# The same param_id can appear in a grouping with multiple signatures, so this
# selects a single item with the best mean rank per grouping (tissue + data transform)
#best_baselines <- best_baselines %>%
#  select(-cor_rank, -rMSE_rank, -mAPE_rank, -type) %>%
#  slice_min(order_by = mean_rank, n = 1) %>%
#  select(-mean_rank) %>%
#  distinct()

# The signature doesn't matter when all percentages are zeros, so we just subset
# to have one unique param_id per tissue/data transform
best_zeros <- subset(best_errors, reference_data_name == "zeros") %>%
  mutate(signature = NA) %>%
  distinct()

best_errors <- subset(best_errors, algorithm != "Baseline" &
                        signature == reference_data_name)

best_errors <- do.call(rbind, list(best_errors, best_baselines, best_zeros))


# Top 3 results again, limited to only errors from the best signatures
err_ranks_best <- best_errors %>%
  subset(algorithm != "Baseline") %>%
  Rank_Errors(group_cols = "tissue") %>%
  Get_TopRanked(n_top = 3)

saveRDS(list("ranked_errors_all" = err_ranks,
             "ranked_errors_best_signatures" = err_ranks_best),
        file.path(dir_analysis, paste0("ranked_errors_", granularity, ".rds")))

group_cols <- c("tissue", "reference_data_name", "test_data_name",
                "normalization", "regression_method", "algorithm")
group_cols_toplevel <- setdiff(group_cols, "reference_data_name")

# Best parameters for each reference_data_name for each possible data input
# (bulk dataset, normalization, regression, input type, and algorithm)
best_params <- Find_BestParameters(best_errors, group_cols)

# Best parameters for each possible data input, ignoring reference_data_name
best_params_toplevel <- Find_BestParameters(best_errors, group_cols_toplevel)

# Only keep error information for the best parameters
best_errors_toplevel <- merge(best_params_toplevel, best_errors,
                              by = c("tissue", "param_id"), all.y = FALSE)

best_errors <- merge(best_params, best_errors, by = c("tissue", "param_id"),
                     all.y = FALSE)

saveRDS(list("best_errors_all" = best_errors,
             "best_errors_toplevel" = best_errors_toplevel,
             "quality_stats" = qstats),
        file.path(dir_analysis, paste0("best_errors_", granularity, ".rds")))

# Get the estimates associated with each parameter ID left in the errors df.
# best_params and best_params_toplevel may not have 100% overlap so we call
# the get estimates function separately for each
best_estimates <- Get_AllBestEstimatesAsDf(bulk_datasets, granularity,
                                           combined_metadata, best_params,
                                           round(n_cores/2))

best_estimates_toplevel <- Get_AllBestEstimatesAsDf(bulk_datasets, granularity,
                                                    combined_metadata,
                                                    best_params_toplevel,
                                                    round(n_cores/2))

gc()

best_errors <- Rank_Errors(best_errors, group_cols)
best_errors_toplevel <- Rank_Errors(best_errors_toplevel, group_cols_toplevel)

best_errors$avg_id <- unlist(apply(best_errors, 1, function(row) {
  paste(row[group_cols], collapse = "_")
}))

best_errors_toplevel$avg_id <- unlist(apply(best_errors_toplevel, 1, function(row) {
  paste(row[group_cols_toplevel], collapse = "_")
}))

# Average the estimates corresponding to best correlation, best rMSE, and best
# mAPE together for each data input type
avg_list <- Create_AveragesList(best_errors, best_estimates, 2) #round(n_cores/3))
avg_list_toplevel <- Create_AveragesList(best_errors_toplevel,
                                         best_estimates_toplevel,
                                         2)

saveRDS(list("average_list_all" = avg_list,
             "average_list_toplevel" = avg_list_toplevel),
        file.path(dir_analysis, paste0("averages_lists_", granularity, ".rds")))

# Calculate significance of cell type differences on a tissue-by-tissue basis
mean_props_all <- Get_MeanProps_Significance(avg_list, 2)
mean_props_toplevel <- Get_MeanProps_Significance(avg_list_toplevel, 2)

saveRDS(list("significance_props_all" = mean_props_all,
             "significance_props_toplevel" = mean_props_toplevel),
        file.path(dir_analysis, paste0("significance_lists_", granularity, ".rds")))

gc()

# This is really slow so we do it last
best_errors_single <- best_errors %>% ungroup()

# Mimic the structure of avg_list as returned from Create_AverageList
best_errors_single$avg_id <- paste(best_errors_single$tissue,
                                   best_errors_single$param_id,
                                   sep = "_")
best_estimates_single <- best_estimates %>%
  dplyr::rename(percent_mean = percent) %>%
  mutate(avg_id = paste(tissue, param_id, sep = "_")) %>%
  select(-param_id)

rm(best_errors, best_errors_toplevel, best_estimates, best_estimates_toplevel)
gc()

avg_list_single <- lapply(unique(best_errors_single$tissue), function(tiss) {
  return(list("avg_errors" = subset(best_errors_single, tissue == tiss),
              "avg_estimates" = subset(best_estimates_single, tissue == tiss)))
})
names(avg_list_single) <- unique(best_errors_single$tissue)

rm(best_estimates_single)
gc()

# Significance of each individual parameter set
mean_props_single <- Get_MeanProps_Significance(avg_list_single, 2)

saveRDS(list("significance_props_all" = mean_props_all,
             "significance_props_toplevel" = mean_props_toplevel,
             "significance_props_single" = mean_props_single),
        file.path(dir_analysis, paste0("significance_lists_", granularity, ".rds")))
