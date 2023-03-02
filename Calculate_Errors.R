library(Matrix)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(Metrics)
library(stringr)
library(dplyr)
library(sparseMatrixStats)

source(file.path("functions", "General_HelperFunctions.R"))
source(file.path("functions", "FileIO_HelperFunctions.R"))

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito",
              "seaRef") #, "seaAD")

granularity = "broad"
datatype = "donors"

est_fields = list("dtangle" = "estimates",
                  "music" = "Est.pctRNA.weighted",
                  #"music2" = "Est.pctRNA",
                  "hspe" = "estimates",
                  "deconRNASeq" = "out.all")

algorithms <- names(est_fields)

err_list <- list()

# Goodness of fit
gof <- function(meas_expr_cpm, est_pct, sig_matrix_cpm) {
  est_expr <- t(est_pct %*% t(sig_matrix_cpm))

  stats <- list()

  meas_expr_cpm <- as.matrix(meas_expr_cpm)
  ok <- meas_expr_cpm != 0 | est_expr != 0

  stats$gof_cor_subject <- mean(diag(cor(meas_expr_cpm, est_expr, use = "na.or.complete")))
  stats$gof_rMSE <- rmse(meas_expr_cpm, est_expr)
  stats$gof_mAPE <- smape(meas_expr_cpm[ok], est_expr[ok])

  return(as.data.frame(stats))
}


for (dataset in datasets) {
  pseudobulk <- Load_Pseudobulk(dataset, datatype, granularity, output_type = "cpm")

  pctRNA <- as.matrix(metadata(pseudobulk)[["pctRNA"]])
  pseudobulk_cpm <- assay(pseudobulk, "counts")

  sig_matrix <- Load_SignatureMatrix(dataset, granularity)

  #### For a single donor (donor 6):
  #sig.matrix.single <- sapply(levels(metadata$broadcelltype), FUN = function(ct) {
  #  cells <- subset(metadata, broadcelltype == ct & donor == colnames(pseudobulk)[6])
  #  if (length(cells) == 0) {
  #    return(rep(0, nrow(pseudobulk)))
  #  }
  #  rowMeans(counts(sce)[,rownames(cells)])
  #})

  #sig.matrix.single <- sweep(sig.matrix.single, 2, colSums(sig.matrix.single), "/") * 1e6

  #t(pct[donor,] %*% t(sig.matrix.single)) will give gene counts == pseudobulk.cpm[,donor]

  ####
  for (algorithm in algorithms) {
    print(str_glue("Calculating errors for {dataset}: {algorithm}"))
    est_field <- est_fields[[algorithm]]

    deconv_list <- Load_AlgorithmOutputList(algorithm, dataset, datatype, granularity)

    # If the file didn't exist, we haven't run the algorithm for this set of
    # params. Skip it.
    if (length(deconv_list) == 0) {
      next
    }

    errs <- list()
    errs_by_celltype <- list()
    errs_by_subject <- list()
    errs_gof <- list()
    params <- list()

    for (ii in names(deconv_list)) {
      tmp <- deconv_list[[ii]][[est_field]]
      tmp <- tmp[rownames(pctRNA), colnames(pctRNA)]

      if (any(is.na(tmp))) {
        print(paste("Param set", algorithm, ii, "has NA values. Skipping..."))
        next
      }

      run_by_celltype <- function(err_fun) {
        sapply(colnames(pctRNA), FUN = function(X) {
          ok <- pctRNA[,X] != 0 | tmp[,X] != 0
          err_fun(pctRNA[ok,X], tmp[ok,X])
        })
      }
      run_by_subject <- function(err_fun) {
        sapply(rownames(pctRNA), FUN = function(X) {
          ok <- pctRNA[X,] != 0 | tmp[X,] != 0
          err_fun(pctRNA[X,ok], tmp[X,ok])
        })
      }

      errs_by_celltype[[ii]] <- data.frame("cor" = diag(cor(tmp, pctRNA, use = "na.or.complete")),
                                           "rMSE" = run_by_celltype(rmse),
                                           "mAPE" = run_by_celltype(smape))


      errs_by_subject[[ii]] <- data.frame("cor" = diag(cor(t(tmp), t(pctRNA), use = "na.or.complete")),
                                          "rMSE" = run_by_subject(rmse),
                                          "mAPE" = run_by_subject(smape))

      # mAPE can't be calculated if both prediction and truth are 0
      ok <- pctRNA != 0 | tmp != 0

      errs[[ii]] <- data.frame("cor_celltype" = mean(diag(cor(tmp, pctRNA, use = "na.or.complete"))),
                               "cor_subject" = mean(diag(cor(t(tmp), t(pctRNA), use = "na.or.complete"))),
                               "rMSE" = rmse(pctRNA, tmp),
                               "mAPE" = smape(pctRNA[ok], tmp[ok]))

      # Goodness of fit
      gof_list <- list()
      for (filter_lvl in 0:3) {
        if (filter_lvl < 3) {
          sig_filt <- FilterSignature(sig_matrix, filter_lvl)
          pseudobulk_filt <- pseudobulk_cpm[rownames(sig_filt), ]
          gof_list[[paste(filter_lvl, -1)]] <- gof(pseudobulk_filt, tmp, sig_filt)
        }
        else {
          for (filt_percent in c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.75, 1.0,
                                 10, 50, 100, 200, 500)) {
            sig_filt <- FilterSignature(sig_matrix, filter_lvl, dataset, "broad", filt_percent)
            pseudobulk_filt <- pseudobulk_cpm[rownames(sig_filt), ]
            gof_list[[paste(filter_lvl, filt_percent)]] <- gof(pseudobulk_filt, tmp, sig_filt)
          }
        }
      }
      errs_gof[[ii]] <- do.call(rbind, gof_list)

      params[[ii]] <- deconv_list[[ii]]$params
      print(str_glue("\t{ii}"), paste(deconv_list[[ii]]$params, collapse = " "))
    }

    err_list[[algorithm]] <- list("means" = do.call(rbind, errs),
                                  "by_celltype" = errs_by_celltype,
                                  "by_subject" = errs_by_subject,
                                  "gof" = errs_gof,
                                  "params" = params)
    print("Done")
  }

  saveRDS(err_list, file = file.path(dir_errors,
                                     str_glue("errors_{dataset}_{datatype}_{granularity}.rds")))
}
