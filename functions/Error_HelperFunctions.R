library(dplyr)
library(Matrix)
library(Metrics)
library(lme4)

source(file.path("functions", "FileIO_HelperFunctions.R"))

CalcGOF_BySample <- function(meas_expr_cpm, est_expr, param_id) {
  meas_expr_cpm <- as.matrix(meas_expr_cpm)

  # Using log2 transformation for correlation so distributions are more gaussian
  cor_sample <- diag(cor(log2(meas_expr_cpm+1), log2(est_expr+1), use = "na.or.complete"))

  rmse_sample <- sapply(colnames(meas_expr_cpm), function(sample) {
    rmse(meas_expr_cpm[,sample], est_expr[,sample])
  })
  mape_sample <- sapply(colnames(meas_expr_cpm), function(sample) {
    ok <- meas_expr_cpm[,sample] != 0 | est_expr[,sample] != 0
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

CalcGOF_BySample_LM <- function(bulk_dataset_name, covariates, meas_expr_log,
                                est_pct, param_id) {
  covariates_ct <- cbind(covariates[colnames(meas_expr_log),],
                         est_pct[colnames(meas_expr_log),])

  formulas <- Load_ModelFormulas(bulk_dataset_name)

  # Replace all biological variables with the cell type proportions instead
  form <- str_replace(formulas$formula_mixed, "~ diagnosis \\+ tissue",
                      paste(colnames(est_pct), collapse = " + "))
  form <- str_replace(form, "sex \\+ ", "")

  # ROSMAP benefits from a linear mixed effect model
  # For speed we're using a fixed effect model for ROSMAP, but saving this code
  # just in case.
  if (bulk_dataset_name == "UNUSED") {#"ROSMAP") {
    form <- paste("expr ~ 0 +", form)

    lm_est <- sapply(rownames(meas_expr_log), function(gene) {
      expr <- meas_expr_log[gene,]

      res <- lmer(as.formula(form), data = covariates_ct)
      return(fitted(res))
    })
  }
  # Mayo and MSBB can use fixed effects only -- using batch as a random effect
  # causes singular fit for a lot of genes, so we don't do that
  else {
    form <- str_replace(form, " \\+ \\(.*", "") # remove random effect term
    form <- paste("t(meas_expr_log) ~ 0 +", form)
    lm_est <- lm(as.formula(form), data = covariates_ct)
    lm_est <- fitted(lm_est)
  }

  # Reverse the log2(cpm+1) transform
  lm_est <- t(2^lm_est-1)

  return(CalcGOF_BySample(2^meas_expr_log-1, lm_est, param_id))
}


CalcGOF_Means <- function(gof_by_sample, bulk_metadata, param_id) {
  numeric_cols <- setdiff(colnames(gof_by_sample), c("param_id", "sample"))
  gof_means_all <- gof_by_sample %>%
                      summarise(across(all_of(numeric_cols), mean),
                                param_id = unique(param_id))

  tissue_assignments <- bulk_metadata[rownames(gof_by_sample), "tissue"]

  gof_by_sample <- cbind(gof_by_sample, tissue_assignments)
  gof_means_tissue <- gof_by_sample %>% group_by(tissue_assignments) %>%
                          summarise(across(all_of(numeric_cols), mean),
                                    param_id = unique(param_id))

  return(list("all_tissue" = gof_means_all,
              "by_tissue" = gof_means_tissue))
}


CalcEstimateStats <- function(all_ests, bulk_metadata) {
  estimate_stats_sample <- all_ests %>% group_by(celltype, sample) %>%
                              summarize(mean_pct = mean(percent),
                                        sd_pct = sd(percent),
                                        rel_sd_pct = sd_pct / mean_pct)

  # average and sd of the *mean* values for each sample
  estimate_stats_all <- estimate_stats_sample %>% group_by(celltype) %>%
                          summarize(mean_pct = mean(mean_pct),
                                    mean_sd_pct = mean(sd_pct),
                                    mean_rel_sd_pct = mean(rel_sd_pct))

  estimate_stats_sample$tissue <- bulk_metadata[estimate_stats_sample$sample, "tissue"]

  estimate_stats_tissue <- estimate_stats_sample %>% group_by(tissue, celltype) %>%
                              summarize(mean_pct = mean(mean_pct),
                                        mean_sd_pct = mean(sd_pct),
                                        mean_rel_sd_pct = mean(rel_sd_pct))

  return(list("by_sample" = estimate_stats_sample,
              "all_tissue" = estimate_stats_all,
              "by_tissue" = estimate_stats_tissue))
}


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


CalcError_MeanIHC_ByTissue <- function(err_list, params, bulk_meta) {
  errs_ihc_props <- err_list[["errs_ihc_props"]][names(params)]
  errs_ihc_props <- lapply(names(errs_ihc_props), function(N) {
    err <- errs_ihc_props[[N]]
    err <- cbind(err, tissue = bulk_meta[rownames(err),"tissue"])
    df <- err %>% group_by(tissue) %>%
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
    err <- cbind(err, tissue = bulk_meta[rownames(err),"tissue"])
    df <- err %>% group_by(tissue) %>%
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


CalcError_MeanByTissue <- function(bulk_dataset, err_list, params, bulk_meta) {
  # TODO add metadata to error df in calculate errors function
  errs_by_sample <- err_list[["gof_sample"]][names(params)]
  errs_by_sample <- lapply(errs_by_sample, function(df) {
    cbind(df, tissue = bulk_meta[rownames(df),"tissue"])
  })

  errs_by_tissue <- lapply(names(errs_by_sample), function(E) {
    df <- errs_by_sample[[E]] %>% group_by(tissue) %>%
      summarize(across(c("cor", "rMSE", "mAPE"), ~ mean(.x, na.rm = TRUE)))
    df$name <- E
    return(df)
  })
  errs_by_tissue <- do.call(rbind, errs_by_tissue)

  if (bulk_dataset == "ROSMAP" & "errs_ihc_props" %in% names(err_list)) {
    errs_ihc <- CalcError_MeanIHC_ByTissue(err_list, params, bulk_meta)
    errs_by_tissue <- merge(errs_by_tissue, errs_ihc, by = c("name", "tissue"))
  }

  return(errs_by_tissue)
}


Filter_Params <- function(err_list, best_params) {
  params <- lapply(err_list[["params"]], function(X) {
    as.list(X %>% select(-test_data_name))
  })

  params_use <- params %in% best_params

  return(params[params_use])
}


Get_AllBestParamsAsDf <- function(reference_datasets, granularity) {
  best_params_all <- lapply(reference_datasets, function(ref_dataset) {
    data <- readRDS(file.path(dir_output, str_glue("best_params_{ref_dataset}_{granularity}.rds")))
    data$reference_data_name <- ref_dataset
    return(data)
  })
  best_params_all <- do.call(rbind, best_params_all)
  return(best_params_all)
}


Get_AllErrorsAsDf <- function(bulk_datasets, reference_datasets, algorithms, granularity, best_params) {
  # Read in all errors for each dataset & data type,
  # put in one big dataframe
  errs_all <- lapply(bulk_datasets, function(bulk_dataset) {
    bulk_se <- Load_BulkData(bulk_dataset)
    bulk_meta <- colData(bulk_se)

    errs_r <- lapply(reference_datasets, function(ref_dataset) {
      errs_a <- lapply(algorithms, function(alg) {
        err_list <- Load_ErrorList(alg, ref_dataset, bulk_dataset, granularity)

        if (length(err_list) == 0) {
          return(NULL)
        }

        params <- Filter_Params(err_list, best_params)
        param_ids <- sapply(params, paste, collapse = " ")

        errs <- err_list[["gof_mean"]]

        if (bulk_dataset == "ROSMAP" & "errs_ihc_props" %in% names(err_list)) {
          errs_ihc <- CalcError_MeanIHC_AllTissues(err_list)
          errs <- cbind(errs, errs_ihc)
        }

        errs <- errs[names(params),]
        errs$tissue <- "All"
        errs$name <- rownames(errs)

        errs_by_tissue <- CalcError_MeanByTissue(bulk_dataset, err_list, params, bulk_meta)

        errs <- rbind(errs, errs_by_tissue)

        if (!("cor_props" %in% colnames(errs))) {
          errs$cor_props <- NA
          errs$rMSE_props <- NA
          errs$mAPE_props <- NA
          errs$cor_pct <- NA
          errs$rMSE_pct <- NA
          errs$mAPE_pct <- NA
        }

        errs$method <- alg
        errs$reference_data_name <- ref_dataset
        errs$test_data_name <- bulk_dataset
        errs$param_id <- param_ids[errs$name]

        return(errs)
      })

      if (all(sapply(errs_a, is.null)) | length(errs_a) == 0) {
        return(NULL)
      }

      return(do.call(rbind, errs_a))
    })

    errs_r <- do.call(rbind, errs_r)
    return(errs_r)
  })

  errs_all <- do.call(rbind, errs_all)

  return(errs_all)
}


Get_AllEstimatesAsDf <- function(ref_dataset, bulk_dataset, algorithms, granularity, params_keep) {
  ests_alg <- list()

  for (algorithm in algorithms) {
    if (algorithm == "random") {
      next
    }
    est_field <- est_fields[[algorithm]]

    ests <- Load_AlgorithmOutputList(algorithm, ref_dataset, bulk_dataset, granularity)
    ests <- ests[names(ests) %in% params_keep]

    if (length(ests) == 0) {
      next
    }

    params <- lapply(ests, function(X) {
      as.list(X[["params"]] %>% select(-test_data_name))
    })

    param_ids <- sapply(params, paste, collapse = " ")

    ests_melt <- lapply(ests, "[[", est_field)

    #if (bulk_dataset == "ROSMAP") {
    #  ests_melt <- lapply(ests_melt, ConvertToROSMAPCelltypes, remove_unused = FALSE)
    #}

    ests_melt <- lapply(names(ests_melt), FUN = function(X) {
      ests_melt[[X]] <- as.data.frame(ests_melt[[X]])
      ests_melt[[X]]$name <- X
      ests_melt[[X]]$sample <- rownames(ests_melt[[X]])
      ests_melt[[X]]$param_id <- param_ids[X]
      return(ests_melt[[X]])
    })
    ests_melt <- do.call(rbind, ests_melt)
    ests_melt <- melt(ests_melt) %>% dplyr::rename(celltype = variable,
                                                   pct_est = value)

    #ests_melt$name <- str_replace(ests_melt$name, "_ROSMAP", "")
    ests_melt$algorithm <- algorithm

    ests_alg[[algorithm]] <- ests_melt
  }

  ests_alg <- do.call(rbind, ests_alg)

  return(ests_alg)
}
