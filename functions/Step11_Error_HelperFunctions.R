library(dplyr)
library(Matrix)
library(Metrics)

source(file.path("functions", "FileIO_HelperFunctions.R"))

# Loading / data preparation functions -----------------------------------------

Get_Signatures <- function(singlecell_datasets, granularity) {
  all_signatures_cpm <- sapply(singlecell_datasets, function(X) {
    Load_SignatureMatrix(X, granularity, "cpm")
  })
  all_signatures_tmm <- sapply(singlecell_datasets, function(X) {
    Load_SignatureMatrix(X, granularity, "tmm")
  })

  return(list(all_signatures_cpm = all_signatures_cpm,
              all_signatures_tmm = all_signatures_tmm))
}


Get_CommonGenes <- function(bulk_datasets, signatures) {
  common_genes <- lapply(bulk_datasets, function(X) {
    bulk_se <- Load_BulkData(X, output_type = "log_cpm", regression_method = "none")
    rownames(bulk_se)
  })
  common_genes <- append(common_genes, lapply(signatures$all_signatures_cpm, rownames))
  common_genes <- purrr::reduce(common_genes, intersect)

  return(common_genes)
}


Get_FilteredSignatures <- function(signatures, common_genes) {
  filtered_sigs <- lapply(signatures, function(sig_objs) {
    lapply(sig_objs, function(X) {
      X[common_genes, ]
    })
  })

  # Make a list that can be indexed by normalization. TPM normalization uses the
  # CPM signature matrix since single cell data is not normalized to TPM.
  filtered_sigs_final <- list("cpm" = filtered_sigs$all_signatures_cpm,
                              "tmm" = filtered_sigs$all_signatures_tmm,
                              "tpm" = filtered_sigs$all_signatures_cpm)

  return(filtered_sigs_final)
}


# Quality control --------------------------------------------------------------

Too_Many_Zeros <- function(est_pct, granularity, zero_thresh = 0.25) {
  # Major cell types are >10% of the population in the Cain dataset
  if (granularity == "broad_class") {
    major_celltypes <- c("Astrocyte", "Excitatory", "Inhibitory", "Oligodendrocyte")

  } else {
    # sub_class uses cell types >5% of the population, plus the most abundant
    # excitatory and inhibitory subclasses
    major_celltypes <- c("Astrocyte", "L2/3/6 IT", "L4/5 IT", "Pvalb", "Oligodendrocyte")
  }

  zeros <- colSums(est_pct == 0)[major_celltypes]
  zero_thresh <- zero_thresh * nrow(est_pct) # 25% can be zeros

  # Baseline never counts as having "too many zeros" because of having
  # deliberately added zeros in some of the estimates.
  if (algorithm != "Baseline" && any(zeros > zero_thresh)) {
    return(TRUE)
  }

  return(FALSE)
}


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
  # Add the tissue variable for easier downstream analysis
  gof_by_sample <- merge(gof_by_sample,
                         select(bulk_metadata, sample, tissue),
                         by = "sample") %>%
    mutate(tissue = as.character(tissue))

  gof_means <- gof_by_sample %>%
    group_by(param_id, tissue, signature) %>%
    summarize(across(where(is.numeric), mean),
              .groups = "drop") %>%
    as.data.frame()

  return(gof_means)
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
