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

# Get the full set of all best errors
best_errors <- Get_AllBestErrorsAsDf(bulk_datasets, granularity, n_cores)

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
  Get_TopRanked(group_cols = "tissue", n_top = 3, with_mean_rank = TRUE)

best_baselines <- subset(best_errors, algorithm == "Baseline" &
                           reference_data_name != "zeros") %>%
  Get_TopRanked(group_cols = c("tissue", "data_transform"), n_top = 1,
                with_mean_rank = TRUE) %>%
  select(-cor_rank, -rMSE_rank, -mAPE_rank, -mean_rank, -type) %>%
  distinct()

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
  Get_TopRanked(group_cols = "tissue", n_top = 3, with_mean_rank = TRUE)

saveRDS(list("ranked_errors_all" = err_ranks,
             "ranked_errors_best_signatures" = err_ranks_best),
        file.path(dir_analysis, str_glue("ranked_errors_{granularity}.rds")))

group_cols <- c("tissue", Get_ParameterColumnNames())  %>%
  setdiff("reference_input_type") # Ignore input type
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
             "best_errors_toplevel" = best_errors_toplevel),
        file.path(dir_analysis, str_glue("best_errors_{granularity}.rds")))

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
avg_list <- Create_AveragesList(best_errors,
                                best_estimates,
                                group_cols,
                                n_cores)
avg_list_toplevel <- Create_AveragesList(best_errors_toplevel,
                                         best_estimates_toplevel,
                                         group_cols_toplevel,
                                         n_cores)

saveRDS(list("average_list_all" = avg_list,
             "average_list_toplevel" = avg_list_toplevel),
        file.path(dir_analysis,
                  str_glue("averages_lists_{granularity}.rds")))

# Calculate significance of cell type differences on a tissue-by-tissue basis
mean_props_all <- mclapply(avg_list,
                           Get_MeanProps_Significance,
                           group_cols = group_cols,
                           mc.cores = n_cores)
mean_props_toplevel <- mclapply(avg_list_toplevel,
                                Get_MeanProps_Significance,
                                group_cols = group_cols_toplevel,
                                mc.cores = n_cores)

saveRDS(list("significance_props_all" = mean_props_all,
             "significance_props_toplevel" = mean_props_toplevel),
        file.path(dir_analysis,
                  str_glue("significance_lists_{granularity}.rds")))
