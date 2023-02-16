library(Matrix)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(Metrics)
library(preprocessCore)
library(stringr)
library(dplyr)
library(sparseMatrixStats)
library(scuttle)

source("Filenames.R")
source(file.path("functions", "General_HelperFunctions.R"))

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito",
              "seaRef") #, "seaAD")

datatype = "donors"

est_fields = list("dtangle" = "estimates",
                  "music_wt" = "Est.prop.weighted",
                  "music_nnls" = "Est.prop.allgene",
                  #"music2" = "Est.prop",
                  "hspe" = "estimates",
                  "deconRNASeq" = "Est.prop")

convert <- list("dtangle" = FALSE,
                "music_wt" = TRUE,
                "music_nnls" = TRUE,
                "music2" = TRUE,
                "hspe" = FALSE,
                "deconRNASeq" = FALSE)

algorithms <- names(est_fields)

err_list <- list()

# Goodness of fit
gof <- function(meas.expr.cpm, est.pct, sig.matrix.cpm) {
  est.expr <- t(est.pct %*% t(sig.matrix.cpm))

  stats <- list()

  meas.expr.cpm <- as.matrix(meas.expr.cpm)
  ok <- meas.expr.cpm != 0 | est.expr != 0

  #meas.expr.log2 <- log2(meas.expr.cpm+1)
  #est.expr.log2 <- log2(est.expr + 1)

  stats$gof.cor_subject <- mean(diag(cor(meas.expr.cpm, est.expr)))
  stats$gof.rMSE <- rmse(meas.expr.cpm, est.expr)
  #stats$gof.mAE <- mae(meas.expr.cpm, est.expr)
  stats$gof.mAPE <- smape(meas.expr.cpm[ok], est.expr[ok])

  return(as.data.frame(stats))
}


for (dataset in datasets) {
  se <- readRDS(file.path(dir_pseudobulk, paste0("pseudobulk_", dataset, "_",
                                                 datatype, "_broadcelltypes.rds")))

  pseudobulk <- assays(se)[["counts"]]
  pseudobulk.cpm <- calculateCPM(pseudobulk)
  pctRNA <- as.matrix(metadata(se)[["pctRNA"]])

  A <- readRDS(file.path(dir_input, str_glue("{dataset}_A_matrix.rds")))
  A <- A[["A_broad"]]

  sig.matrix <- readRDS(file.path(dir_input, str_glue("{dataset}_signature.rds")))
  sig.matrix <- sig.matrix[["sig_broad"]]

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

  # TODO figure out filter level for goodness of fit
  #sig.matrix <- FilterSignature(sig.matrix, filter_level = 2)

  #pseudobulk.filt <- pseudobulk.cpm[rownames(sig.matrix), ]

  # Unclear if this is necessary -- probably no?
  #qn <- cbind(sig.matrix, pseudobulk.filt)
  #qn <- normalize.quantiles(as.matrix(qn), copy = FALSE)

  #sig.matrix <- as.matrix(qn[,1:ncol(sig.matrix)])
  #pseudobulk.filt <- as.matrix(qn[,-c(1:ncol(sig.matrix))])

  for (algorithm in algorithms) {
    print(str_glue("Calculating errors for {dataset}: {algorithm}"))
    est.field <- est_fields[[algorithm]]

    alg_name <- algorithm
    if (algorithm == "music_wt" | algorithm == "music_nnls") {
      alg_name <- "music"
    }

    params_file <- file.path(dir_params_lists,
                             paste0(alg_name, "_list_",  dataset, "_",
                                    datatype, "_broad.rds"))

    if (!file.exists(params_file)) {
      next
    }

    deconv_list <- readRDS(params_file)

    # TODO temporary for dtangle & hspe
    if (algorithm == "dtangle" | algorithm == "hspe") {
      params_file2 <- file.path(dir_params_lists,
                                str_glue("{alg_name}_list_{dataset}_{datatype}_broad_input_pseudobulk.rds"))
      if (file.exists(params_file2)) {
        deconv_list2 <- readRDS(params_file2)
        deconv_list <- append(deconv_list, deconv_list2)
        rm(deconv_list2)
      }
    }
    errs <- list()
    errs_by_celltype <- list()
    errs_by_subject <- list()
    errs_gof <- list()

    ###step 1: find correlations between predictions and pseudobulk proportions
    for (ii in names(deconv_list)) {
      if (!(any(is.na(deconv_list[[ii]][[est.field]])))) {  ###exclude all predictions that return any NA values on the pseudobulk
        tmp <- deconv_list[[ii]][[est.field]]
        tmp <- tmp[rownames(pctRNA), colnames(pctRNA)]

        # MuSiC outputs percent of cells instead of percent RNA. Convert to
        # percent RNA using the estimated A matrix
        if (convert[[algorithm]] == TRUE) {
          tmp <- sweep(tmp, 2, A, "*")
          tmp <- sweep(tmp, 1, rowSums(tmp), "/")
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

        errs_by_celltype[[ii]] <- data.frame("cor" = diag(cor(tmp, pctRNA)),
                                             "rMSE" = run_by_celltype(rmse),
                                             #"mAE" = run_by_celltype(mae),
                                             "mAPE" = run_by_celltype(smape))


        errs_by_subject[[ii]] <- data.frame("cor" = diag(cor(t(tmp), t(pctRNA))),
                                            "rMSE" = run_by_subject(rmse),
                                            #"mAE" = run_by_subject(mae),
                                            "mAPE" = run_by_subject(smape))

        # mAPE can't be calculated if both prediction and truth are 0
        ok <- pctRNA != 0 | tmp != 0

        errs[[ii]] <- data.frame("cor_celltype" = mean(diag(cor(tmp, pctRNA))),
                                 "cor_subject" = mean(diag(cor(t(tmp), t(pctRNA)))),
                                 #"cor_overall" = cor(as.vector(tmp), as.vector(pctRNA)),
                                 "rMSE" = rmse(pctRNA, tmp),
                                 #"mAE" = mae(pctRNA, tmp),
                                 "mAPE" = smape(pctRNA[ok], tmp[ok]))

        # Goodness of fit
        gof_list <- list()
        for (filter_lvl in 0:3) {
          if (filter_lvl < 3) {
            sig.filt <- FilterSignature(sig.matrix, filter_lvl)
            pseudobulk.filt <- pseudobulk.cpm[rownames(sig.filt), ]
            gof_list[[paste(filter_lvl, 1)]] <- gof(pseudobulk.filt, tmp, sig.filt)
          }
          else {
            for (filt_percent in c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.75, 1.0)) {
              sig.filt <- FilterSignature(sig.matrix, filter_lvl, dataset, "broad", filt_percent)
              pseudobulk.filt <- pseudobulk.cpm[rownames(sig.filt), ]
              gof_list[[paste(filter_lvl, filt_percent)]] <- gof(pseudobulk.filt, tmp, sig.filt)
            }
          }
        }
        errs_gof[[ii]] <- do.call(rbind, gof_list)
      }
      print(str_glue("\t{ii}"))
    }

    err_list[[algorithm]] <- list("means" = do.call(rbind, errs),
                                  "by_celltype" = errs_by_celltype,
                                  "by_subject" = errs_by_subject,
                                  "gof" = errs_gof)
    print("Done")
  }

  saveRDS(err_list, file = file.path(dir_output,
                                     paste0("errors_", dataset, "_",
                                            datatype, "_broad.rds")))
}
