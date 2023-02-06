library(Matrix)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(Metrics)
library(preprocessCore)
library(stringr)
library(dplyr)
library(sparseMatrixStats)

source("Filenames.R")
source(file.path("functions", "General_HelperFunctions.R"))

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito",
              "seaRef") #, "seaAD")

datatype = "training"

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
  est.expr <- t(est.pct %*% t(sig.matrix))

  stats <- list()

  meas.expr.cpm <- as.matrix(meas.expr.cpm)
  ok <- meas.expr.cpm != 0 | est.expr != 0

  stats$gof.cor <- mean(diag(cor(log2(meas.expr.cpm + 1), log2(est.expr + 1))))
  stats$gof.rMSE <- rmse(meas.expr.cpm, est.expr)
  stats$gof.mAE <- mae(meas.expr.cpm, est.expr)
  stats$gof.mAPE <- smape(meas.expr.cpm[ok], est.expr[ok])

  return(as.data.frame(stats))
}


for (dataset in datasets) {
  se <- readRDS(file.path(dir_pseudobulk, paste0("pseudobulk_", dataset, "_",
                                                 datatype, "_broadcelltypes.rds")))

  pseudobulk <- assays(se)[["counts"]]
  pseudobulk.cpm <- scuttle::calculateCPM(pseudobulk)
  pctRNA <- as.matrix(metadata(se)[["pctRNA"]])

  sce <- readRDS(file.path(dir_input, paste0(dataset, "_sce.rds")))

  sc.cpm <- scuttle::calculateCPM(counts(sce))
  metadata <- colData(sce)

  A <- CalculateA(sce, metadata$donor, metadata$broadcelltype)

  # Mean expression of all genes for each cell type -- in cpm
  # TODO should this be mean cpms? or mean counts then converted to cpms?
  sig.matrix <- CalculateSignature(sc.cpm, metadata$donor, metadata$broadcelltype)

  # standard deviation
  #sd.list <- lapply(levels(metadata$broadcelltype), FUN = function(ct) {
  #  cells <- subset(metadata, broadcelltype == ct)
  #  rowSds(sc.cpm[,rownames(cells)], useNames = TRUE)
  #})

  #names(sd.list) <- levels(metadata$broadcelltype)
  #sd.matrix <- do.call(cbind, sd.list)

  #### For a single donor (donor 6):
  #sig.matrix.single <- sapply(levels(metadata$broadcelltype), FUN = function(ct) {
  #  cells <- subset(metadata, broadcelltype == ct & donor == colnames(pseudobulk)[6])
  #  if (length(cells) == 0) {
  #    return(rep(0, nrow(pseudobulk)))
  #  }
  #  rowMeans(sc.cpm[,rownames(cells)])
  #})

  #names(sig.list.single) <- levels(metadata$broadcelltype)
  #sig.matrix.single <- do.call(cbind, sig.list.single)

  #sig.matrix.single <- sweep(sig.matrix.single, 2, colSums(sig.matrix.single), "/") * 1e6

  #t(pct[donor,] %*% t(sig.matrix.single)) will give gene counts == pseudobulk.cpm[,donor]

  ####

  # TODO is this the best way to pick genes? We obviously can't use this when
  # running on real bulk data...
  #res <- t(pctRNA %*% t(sig.matrix))
  #errs <- sapply(rownames(res), function(gene) {
  #  ok <- pseudobulk.cpm[gene,] != 0 | res[gene,] != 0
  #  smape(pseudobulk.cpm[gene, ok], res[gene, ok])
  #})

  #good.genes <- names(errs)[which(errs <= 0.2)] # Genes with less than 10% error on avg (smape ~= 2*error)
  #sig.matrix <- sig.matrix[good.genes,]
  # Filter for genes where at least one cell type expresses at > 1 cpm
  ok <- which(rowSums(sig.matrix >= 1) > 0)
  sig.matrix <- sig.matrix[ok, ]

  pseudobulk.filt <- pseudobulk.cpm[rownames(sig.matrix), ]

  # Unclear if this is necessary -- probably yes?
  #qn <- cbind(sig.matrix, pseudobulk.filt)
  #qn <- normalize.quantiles(as.matrix(qn), copy = FALSE)

  #sig.matrix <- as.matrix(qn[,1:ncol(sig.matrix)])
  #pseudobulk.filt <- as.matrix(qn[,-c(1:ncol(sig.matrix))])

  for (algorithm in algorithms) {
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
    errs <- list()
    errs_by_celltype <- list()
    errs_by_subject <- list()

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

        # mAPE can't be calculated if both prediction and truth are 0
        ok <- pctRNA != 0 | tmp != 0

        run_by_celltype <- function(err_fun) {
          sapply(colnames(pctRNA), FUN = function(X) {
            err_fun(pctRNA[,X], tmp[,X])
          })
        }
        run_by_subject <- function(err_fun) {
          sapply(rownames(pctRNA), FUN = function(X) {
            err_fun(pctRNA[X,], tmp[X,])
          })
        }

        errs_by_celltype[[ii]] <- data.frame("cor" = diag(cor(tmp, pctRNA)),
                                             "rMSE" = run_by_celltype(rmse),
                                             "mAE" = run_by_celltype(mae),
                                             "mAPE" = run_by_celltype(smape))


        errs_by_subject[[ii]] <- data.frame("cor" = diag(cor(t(tmp), t(pctRNA))),
                                            "rMSE" = run_by_subject(rmse),
                                            "mAE" = run_by_subject(mae),
                                            "mAPE" = run_by_subject(smape))

        errs[[ii]] <- data.frame("cor_celltype" = mean(diag(cor(tmp, pctRNA))),
                                 "cor_subject" = mean(diag(cor(t(tmp), t(pctRNA)))),
                                 "cor_overall" = cor(as.vector(tmp), as.vector(pctRNA)),
                                 "rMSE" = rmse(pctRNA, tmp),
                                 "mAE" = mae(pctRNA, tmp),
                                 "mAPE" = smape(pctRNA[ok], tmp[ok]))

        # Goodness of fit
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
