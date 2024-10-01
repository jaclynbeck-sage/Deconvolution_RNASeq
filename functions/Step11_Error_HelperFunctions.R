library(dplyr)
library(Matrix)
library(Metrics)

source(file.path("functions", "FileIO_HelperFunctions.R"))

# Goodness-of-fit error calculations -------------------------------------------

# Calculate goodness of fit (GOF) for each sample across 3 metrics (correlation,
# root mean squared error, and mean absolute percent error).
#
# Arguments:
#   meas_expr_cpm - the observed bulk data matrix, on the linear scale (CPM-,
#                   TMM-, or TPM-normalized values depending on the parameters)
#   est_expr - the estimated expression matrix calculated from multiplying a
#              signature matrix by the estimated cell type percentages, on the
#              linear scale
#   param_id - a string uniquely identifying this set of parameters
#
# Returns:
#   a data frame where rows are samples and columns are the error metrics plus
#   some metadata
CalcGOF_BySample <- function(meas_expr_cpm, est_expr, param_id) {
  meas_expr_cpm <- as.matrix(meas_expr_cpm)

  # Using log2 transformation for correlation so distributions are more gaussian
  cor_sample <- sapply(colnames(meas_expr_cpm), function(sample) {
    cor(log2(meas_expr_cpm[, sample] + 1), log2(est_expr[, sample] + 1),
        use = "na.or.complete")
  })

  rmse_sample <- sapply(colnames(meas_expr_cpm), function(sample) {
    rmse(meas_expr_cpm[, sample], est_expr[, sample])
  })

  mape_sample <- sapply(colnames(meas_expr_cpm), function(sample) {
    ok <- meas_expr_cpm[, sample] != 0 | est_expr[, sample] != 0
    smape(meas_expr_cpm[ok, sample], est_expr[ok, sample]) / 2
  })

  gof_by_sample <- data.frame(cor = cor_sample,
                              rMSE = rmse_sample,
                              mAPE = mape_sample,
                              param_id = param_id,
                              sample = names(cor_sample),
                              row.names = names(cor_sample))

  return(gof_by_sample)
}


# Takes the mean goodness-of-fit across all samples, for each error metric
#
# Arguments:
#   gof_by_sample - a data frame where rows are samples and columns are the
#                   error metrics plus some metadata, as returned by CalcGOF_BySample
#   bulk_metadata - a data frame of metadata where rows are samples, and one
#                   column is "tissue"
#   param_id - a string uniquely identifying this set of parameters
#
# Returns:
#   a data frame containing one row for the means across all samples, plus one
#   row for each tissue containing the means across all samples in that tissue,
#   for each signature used in the error calculation
CalcGOF_Means <- function(gof_by_sample, bulk_metadata, param_id) {
  # take the mean of each numerical column across all samples
  gof_means_all <- gof_by_sample %>%
    summarise(across(where(is.numeric), mean),
              param_id = unique(param_id),
              .groups = "drop")

  gof_means_all$tissue <- "All"

  tissue <- bulk_metadata[rownames(gof_by_sample), "tissue"]

  gof_by_sample <- cbind(gof_by_sample, tissue)

  # mean of numerical columns for samples within a given tissue
  gof_means_tissue <- gof_by_sample %>%
    group_by(tissue) %>%
    summarise(across(where(is.numeric), mean),
              param_id = unique(param_id),
              .groups = "drop")

  gof_means <- rbind(gof_means_all, gof_means_tissue[, colnames(gof_means_all)])
  gof_means$signature <- unique(gof_by_sample$signature)

  return(gof_means)
}


# Calculates statistics about goodness-of-fit data:
#   - the mean and SD of estimated percentages for each cell type for a given
#     sample across all estimates in the file
#   - the same mean and SD except only across the top 10-scoring estimates in
#     the file
#   - the mean of these means across all samples
#
# Arguments:
#   all_ests - a melted data frame with columns "param_id", "sample", "celltype",
#              and "percent", representing the estimates for each cell type for
#              each sample
#   bulk_metadata - a data frame of metadata, where rows are bulk samples and
#                   one column is "tissue"
#   gof_means_all - a data frame with N rows per parameter set (where N = [1 +
#                   the number of tissues] * [the number of signatures]), and
#                   columns are the mean error metrics as returned by CalcGOF_Means
#
# Returns:
#   a list of statistics, containing items "by_sample", "all_tissue",
#   "by_tissue" (all data frames), and "top_10_stats" (a list, one entry for
#   each error metric)
CalcEstimateStats <- function(all_ests, bulk_metadata, gof_means_all) {
  # Mean of each cell type percentage over all estimates in the file, per
  # celltype per sample
  estimate_stats_sample <- all_ests %>%
    group_by(celltype, sample) %>%
    summarize(mean_pct = mean(percent),
              sd_pct = sd(percent),
              rel_sd_pct = sd_pct / mean_pct,
              .groups = "drop_last")

  estimate_stats_sample$tissue <- bulk_metadata[estimate_stats_sample$sample, "tissue"]

  # Top 10 scoring parameter sets per tissue, per signature used to
  # calculate the error. Each error metric may have different top 10 lists. If
  # there are less than 10 parameter sets, this just selects all of them.
  top_cor <- gof_means_all %>% group_by(tissue, signature) %>% top_n(10, wt = cor)
  top_rMSE <- gof_means_all %>% group_by(tissue, signature) %>% top_n(10, wt = -rMSE)
  top_mAPE <- gof_means_all %>% group_by(tissue, signature) %>% top_n(10, wt = -mAPE)

  top_10 <- list(cor = top_cor,
                 rMSE = top_rMSE,
                 mAPE = top_mAPE)

  ests_tmp <- all_ests
  ests_tmp$tissue <- bulk_metadata[ests_tmp$sample, "tissue"]

  # Stats for each of the error metrics on the top 10 scoring parameter sets
  # for each tissue
  top_10_stats <- lapply(top_10, function(top_errs) {
    # Calculate stats for each tissue (which may be a specific tissue or "All")
    ests_top_10 <- lapply(unique(top_errs$tissue), function(tiss) {
      top_tmp <- subset(top_errs, tissue == tiss)

      top_stats <- lapply(unique(top_tmp$signature), function(sig) {
        top_sub <- subset(top_tmp, signature == sig)
        ests_filt <- subset(ests_tmp, param_id %in% top_sub$param_id)

        if (tiss != "All") {
          ests_filt <- subset(ests_filt, tissue == tiss)
        }

        # Mean cell type percent for each cell type for each sample, across the
        # top 10 scoring parameter sets
        ests_stats <- ests_filt %>%
          group_by(celltype, sample) %>%
          summarize(mean_pct = mean(percent),
                    sd_pct = sd(percent),
                    rel_sd_pct = sd_pct / mean_pct,
                    .groups = "drop_last")

        ests_stats$tissue <- tiss
        ests_stats$signature <- sig
        return(ests_stats)
      })

      top_stats <- do.call(rbind, top_stats)
    })

    ests_top_10 <- do.call(rbind, ests_top_10)
    return(ests_top_10)
  })

  # average and sd of the *mean* values for each sample, per cell type
  estimate_stats_all <- estimate_stats_sample %>%
    group_by(celltype) %>%
    summarize(mean_pct = mean(mean_pct),
              mean_sd_pct = mean(sd_pct),
              mean_rel_sd_pct = mean(rel_sd_pct),
              .groups = "drop_last")

  # Same but broken down by tissue
  estimate_stats_tissue <- estimate_stats_sample %>%
    group_by(tissue, celltype) %>%
    summarize(mean_pct = mean(mean_pct),
              mean_sd_pct = mean(sd_pct),
              mean_rel_sd_pct = mean(rel_sd_pct),
              .groups = "drop_last")

  return(list("by_sample" = estimate_stats_sample,
              "all_tissue" = estimate_stats_all,
              "by_tissue" = estimate_stats_tissue,
              "top_10_stats" = top_10_stats))
}


# IHC error functions ----------------------------------------------------------

# NOTE: This function has not been used or tested in a long time so it probably
# doesn't work.
CalcError_MeanIHC_AllTissues <- function(err_list) {
  # IHC props and percent errors are per individual, get the mean for each
  # error type across all individuals
  errs_ihc_props <- lapply(err_list[["errs_ihc_props"]], function(err) {
    err %>% summarize(across(colnames(err), ~ mean(.x, na.rm = TRUE)))
  })

  errs_ihc_props <- do.call(rbind, errs_ihc_props) %>%
    dplyr::rename(cor_props = cor,
                  rMSE_props = rMSE,
                  mAPE_props = mAPE)

  errs_ihc_pct <- lapply(err_list[["errs_ihc_pct"]], function(err) {
    err %>% summarize(across(colnames(err), ~ mean(.x, na.rm = TRUE)))
  })

  errs_ihc_pct <- do.call(rbind, errs_ihc_pct) %>%
    dplyr::rename(cor_pct = cor,
                  rMSE_pct = rMSE,
                  mAPE_pct = mAPE)

  return(cbind(errs_ihc_props, errs_ihc_pct))
}


# NOTE: This function has not been used or tested in a long time so it probably
# doesn't work.
CalcError_MeanIHC_ByTissue <- function(err_list, params, bulk_meta) {
  errs_ihc_props <- err_list[["errs_ihc_props"]][names(params)]
  errs_ihc_props <- lapply(names(errs_ihc_props), function(N) {
    err <- errs_ihc_props[[N]]
    err <- cbind(err, tissue = bulk_meta[rownames(err), "tissue"])

    df <- err %>%
      group_by(tissue) %>%
      summarize(across(!contains("tissue"), ~ mean(.x, na.rm = TRUE)))
    df$name <- N

    return(df)
  })

  errs_ihc_props <- do.call(rbind, errs_ihc_props) %>%
    dplyr::rename(cor_props = cor,
                  rMSE_props = rMSE,
                  mAPE_props = mAPE)

  errs_ihc_pct <- err_list[["errs_ihc_pct"]][names(params)]

  errs_ihc_pct <- lapply(names(errs_ihc_pct), function(N) {
    err <- errs_ihc_pct[[N]]
    err <- cbind(err, tissue = bulk_meta[rownames(err), "tissue"])

    df <- err %>%
      group_by(tissue) %>%
      summarize(across(!contains("tissue"), ~ mean(.x, na.rm = TRUE)))
    df$name <- N

    return(df)
  })

  errs_ihc_pct <- do.call(rbind, errs_ihc_pct) %>%
    dplyr::rename(cor_pct = cor,
                  rMSE_pct = rMSE,
                  mAPE_pct = mAPE)

  errs_ihc <- merge(errs_ihc_props, errs_ihc_pct, by = c("name", "tissue"))
  return(errs_ihc)
}
