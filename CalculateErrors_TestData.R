library(Matrix)
library(SingleCellExperiment)
library(Metrics)
library(dplyr)

source("Filenames.R")
source(file.path("functions", "General_HelperFunctions.R"))

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito",
              "seaRef") #, "seaAD")

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
  sce <- readRDS(file.path(dir_input, paste(dataset, "sce.rds", sep = "_")))
  metadata <- colData(sce)

  # Convert gene names to Ensembl IDs
  genes <- rowData(sce)
  rownames(sce) <- genes[rownames(sce), "Ensembl.ID"]

  sce_cpm <- scuttle::calculateCPM(counts(sce))

  # Specific to ROSMAP for now
  bulk <- read.table(file.path(dir_input, "ROSMAP_DLPFC_Counts.tsv"), sep = "\t",
                     header = TRUE)
  rownames(bulk) <- bulk$ensembl_gene_id
  bulk <- bulk[,-1]

  keepgene <- intersect(rownames(sce),rownames(bulk))

  bulk_cpm <- scuttle::calculateCPM(bulk)

  bulk_cpm <- bulk_cpm[keepgene,]
  sce_cpm <- sce_cpm[keepgene,]

  sig.matrix <- CalculateSignature(sce_cpm, metadata$donor, metadata$broadcelltype)
  ok <- which(rowSums(sig.matrix >= 1) > 0)
  sig.matrix <- sig.matrix[ok, ]

  bulk.filt <- bulk_cpm[rownames(sig.matrix),]

  # TODO test whether this helps with GOF
  #qn <- cbind(sig.matrix, bulk.filt)
  #qn <- normalize.quantiles(as.matrix(qn), copy = FALSE)

  #sig.matrix <- as.matrix(qn[,1:ncol(sig.matrix)])
  #bulk.filt <- as.matrix(qn[,-c(1:ncol(sig.matrix))])

  for (algorithm in algorithms) {
    est.field <- est_fields[[algorithm]]

    alg_name <- algorithm
    if (algorithm == "music_wt" | algorithm == "music_nnls") {
      alg_name <- "music"
    }

    params_file <- file.path(dir_rosmap,
                             paste0(alg_name, "_list_",  dataset,
                                    "_broad_ROSMAP.rds"))

    if (!file.exists(params_file)) {
      next
    }

    deconv_list <- readRDS(params_file)
    errs <- list()

    ###step 1: find correlations between predictions and pseudobulk proportions
    for (ii in names(deconv_list)) {
      if (!(any(is.na(deconv_list[[ii]][[est.field]])))) {  ###exclude all predictions that return any NA values on the pseudobulk
        tmp <- deconv_list[[ii]][[est.field]]
        tmp <- tmp[colnames(bulk.filt), colnames(sig.matrix)]

        # MuSiC outputs percent of cells instead of percent RNA. Convert to
        # percent RNA using the estimated A matrix
        if (convert[[algorithm]] == TRUE) {
          A <- readRDS(file.path(dir_input, paste0(dataset, "_A_matrix.rds")))
          tmp <- ConvertPropCellsToPctRNA(tmp, A[["A_broad"]])
        }

        errs[[ii]] <- gof(bulk.filt, tmp, sig.matrix)
      }
    }

    err_list[[algorithm]] <- do.call(rbind, errs)
    print(c(dataset, algorithm))
  }

  saveRDS(err_list, file = file.path(dir_output,
                                     paste0("errors_", dataset,
                                            "_broad_ROSMAP.rds")))
}
