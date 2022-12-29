library(Matrix)
library(SummarizedExperiment)
library(Metrics)
library(preprocessCore)

source("Filenames.R")

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito") #,
#"seaRef") #, "seaAD")
datatype = "training"

est_fields = list("dtangle" = "estimates",
                  "music" = "Est.prop.weighted", # TODO analyze Est.prop.allgene too
                  #"music2" = "Est.prop",
                  #"hspe" = "estimates",
                  "deconRNASeq" = "Est.prop")

algorithms <- names(est_fields)

err_list <- list()

# Goodness of fit
gof <- function(meas.expr.cpm, est.prop, sig.matrix.cpm) {
  est.expr <- t(est.prop %*% t(sig.matrix))

  stats <- list()

  stats$gof.cor <- mean(diag(cor(log2(meas.expr.cpm + 1), log2(est.expr + 1))))
  stats$gof.rMSE <- rmse(meas.expr.cpm, est.expr)
  stats$gof.mAE <- mae(meas.expr.cpm, est.expr)
  stats$gof.mAPE <- smape(meas.expr.cpm, est.expr)

  return(as.data.frame(stats))
}


for (dataset in datasets) {
  se <- readRDS(file.path(dir_pseudobulk, paste0("pseudobulk_", dataset, "_",
                                                 datatype, "_broadcelltypes.rds")))

  pseudobulk <- assays(se)[["counts"]]
  pseudobulk.cpm <- scuttle::calculateCPM(pseudobulk)
  propval <- as.matrix(colData(se))

  sce <- readRDS(file.path(dir_input, paste0(dataset, "_sce.rds")))
  sc.cpm <- scuttle::calculateCPM(counts(sce))
  metadata <- colData(sce)

  # Mean expression of all genes for each cell type -- in cpm
  sig.list <- lapply(levels(metadata$broadcelltype), FUN = function(ct) {
    cells <- subset(metadata, broadcelltype == ct)
    rowMeans(sc.cpm[,rownames(cells)])
  })

  names(sig.list) <- levels(metadata$broadcelltype)
  sig.matrix <- do.call(cbind, sig.list)

  # Filter for genes where at least one cell type expresses at > 1 cpm
  ok <- which(rowSums(sig.matrix >= 1) > 0)
  sig.matrix <- sig.matrix[ok, ]

  pseudobulk.filt <- pseudobulk.cpm[rownames(sig.matrix), ]

  # Unclear if this is necessary
  qn <- cbind(sig.matrix, pseudobulk.filt)
  qn <- normalize.quantiles(as.matrix(qn), copy = FALSE)

  sig.matrix <- as.matrix(qn[,1:ncol(sig.matrix)])
  pseudobulk.filt <- as.matrix(qn[,-c(1:ncol(sig.matrix))])

  for (algorithm in algorithms) {
    est.field <- est_fields[[algorithm]]
    deconv_list <- readRDS(file.path(dir_params_lists,
                                     paste0(algorithm, "_list_",  dataset, "_",
                                            datatype, "_broad.rds")))
    errs <- list()
    errs_by_celltype <- list()
    errs_by_subject <- list()

    ###step 1: find correlations between predictions and pseudobulk proportions
    for (ii in names(deconv_list)) {
      if (!(any(is.na(deconv_list[[ii]][[est.field]])))) {  ###exclude all predictions that return any NA values on the pseudobulk
        tmp <- deconv_list[[ii]][[est.field]]
        tmp <- tmp[rownames(propval), colnames(propval)]

        # mAPE can't be calculated if both prediction and truth are 0
        ok <- propval != 0 | tmp != 0

        run_by_celltype <- function(err_fun) {
          sapply(colnames(propval), FUN = function(X) {
            err_fun(propval[,X], tmp[,X])
          })
        }
        run_by_subject <- function(err_fun) {
          sapply(rownames(propval), FUN = function(X) {
            err_fun(propval[X,], tmp[X,])
          })
        }

        errs_by_celltype[[ii]] <- data.frame("cor" = diag(cor(tmp,propval)),
                                             "rMSE" = run_by_celltype(rmse),
                                             "mAE" = run_by_celltype(mae),
                                             "mAPE" = run_by_celltype(smape))


        errs_by_subject[[ii]] <- data.frame("cor" = diag(cor(t(tmp), t(propval))),
                                            "rMSE" = run_by_subject(rmse),
                                            "mAE" = run_by_subject(mae),
                                            "mAPE" = run_by_subject(smape))

        errs[[ii]] <- data.frame("cor_celltype" = mean(diag(cor(tmp,propval))),
                                 "cor_subject" = mean(diag(cor(t(tmp), t(propval)))),
                                 "cor_overall" = cor(as.vector(tmp), as.vector(propval)),
                                 "rMSE" = rmse(propval, tmp),
                                 "mAE" = mae(propval, tmp),
                                 "mAPE" = smape(propval[ok], tmp[ok]))

        # Ignore goodness of fit for testing
        #errs[[ii]] <- cbind(errs[[ii]], gof(pseudobulk.filt, tmp, sig.matrix))
      }
    }

    err_list[[algorithm]] <- list("means" = do.call(rbind, errs),
                                  "by_celltype" = errs_by_celltype,
                                  "by_subject" = errs_by_subject)
    print(c(dataset, algorithm))
  }

  saveRDS(err_list, file = file.path(dir_output,
                                     paste0("errors_", dataset, "_",
                                            datatype, "_broad.rds")))
}
