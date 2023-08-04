library(Matrix)
library(SummarizedExperiment)
library(Metrics)
library(dplyr)
library(scuttle)

source(file.path("functions", "General_HelperFunctions.R"))
source(file.path("functions", "FileIO_HelperFunctions.R"))

datasets <- c("cain", "lau", "leng", "mathys", "morabito")#, #"seaRef") #, "seaAD")
granularity <- "broad"
bulk_dataset <- "ROSMAP"

random_types <- c("uniform", "mean_singlecell")

err_list <- list()

# Goodness of fit
gof <- function(meas_expr_cpm, est_pct, sig_matrix_cpm) {
  est_expr <- t(est_pct %*% t(sig_matrix_cpm))

  gof_means <- list()

  meas_expr_cpm <- as.matrix(meas_expr_cpm)

  # TODO this should probably be log2(cpm) for correlation
  cor_subject <- diag(cor(meas_expr_cpm, est_expr, use = "na.or.complete"))
  rmse_subject <- sapply(colnames(meas_expr_cpm), function(sample) {
    rmse(meas_expr_cpm[,sample], est_expr[,sample])
  })
  mape_subject <- sapply(colnames(meas_expr_cpm), function(sample) {
    ok <- meas_expr_cpm[,sample] != 0 | est_expr[,sample] != 0
    smape(meas_expr_cpm[ok, sample], est_expr[ok, sample]) / 2
  })

  gof_by_subject <- data.frame(cor = cor_subject,
                               rMSE = rmse_subject,
                               mAPE = mape_subject,
                               row.names = names(cor_subject))

  gof_means <- data.frame(cor = mean(cor_subject),
                          rMSE = mean(rmse_subject),
                          mAPE = mean(mape_subject))

  return(list("means" = gof_means, "by_subject" = gof_by_subject))
}


for (dataset in datasets) {
  sig_matrix <- Load_SignatureMatrix(dataset, granularity)
  #comp_genes <- readRDS(file.path(dir_input, "compartment_specific_genes.rds"))
  #sig_matrix <- sig_matrix[!(rownames(sig_matrix) %in% comp_genes$Symbol),]

  bulk <- Load_BulkData(bulk_dataset, output_type = "cpm")
  bulk_cpm <- assay(bulk, "counts")

  keepgene <- intersect(rownames(sig_matrix), rownames(bulk_cpm))

  # TODO is it correct to re-calculate CPM?
  sig_matrix <- calculateCPM(sig_matrix[keepgene,])
  bulk_cpm <- calculateCPM(bulk_cpm[keepgene,])

  # TODO this needs to be a util function
  if (bulk_dataset == "ROSMAP") {
    ihc_props <- as.matrix(read.csv(file_rosmap_ihc_proportions, row.names = 1))
    #rownames(ihc_props) <- make.names(rownames(ihc_props))
    A <- Load_AvgLibSize(dataset, granularity)
    A2 <- A

    # Some datasets are missing some vascular types
    for (col in c("Endo", "Peri", "VLMC")) {
      if (!(col %in% names(A2))) {
        A2 <- c(A2, 0)
        names(A2)[length(A2)] <- col
      }
    }

    # TODO this isn't quite right and the A matrix needs to be re-processed
    A2["Neuro"] <- mean(A2[c("Exc", "Inh")])
    A2["Oligo"] <- mean(A2[c("Oligo","OPC")])
    A2["Endo"] <- max(A2[c("Endo", "Peri", "VLMC")])

    A2 <- A2[colnames(ihc_props)]
    A2 <- A2 / sum(A2)
    ihc_pct <- ConvertPropCellsToPctRNA(ihc_props, A2)

    specimen2projid <- subset(colData(bulk), donor %in% rownames(ihc_props))
  }

  errs_gof_mean <- list()
  errs_gof_subject <- list()
  errs_ihc_props <- list()
  errs_ihc_pct <- list()
  params <- list()

  sig_filt <- FilterSignature(sig_matrix, 3, dataset, granularity, 1,
                              marker_type = "dtangle",
                              marker_subtype = "diff",
                              marker_input_type = "pseudobulk")

  keepgene <- intersect(rownames(sig_filt), rownames(bulk))
  bulk_filt <- bulk_cpm[keepgene,]
  sig_filt <- sig_filt[keepgene,]

  for (rand_type in random_types) {
    # Randomly generate proportions based on uniform distribution
    if (rand_type == "uniform") {
      est_pct <- matrix(runif(ncol(bulk_filt) * ncol(sig_filt)),
                        nrow = ncol(bulk_filt))
      est_pct <- sweep(est_pct, 1, rowSums(est_pct), "/") # rows sum to 1
    }
    # Use the mean of the single cell dataset (not random)
    else if (rand_type == "mean_singlecell") {
      pb <- Load_Pseudobulk(dataset, "donors", granularity)
      pcts <- colMeans(pb@metadata[["pctRNA"]])
      est_pct <- matrix(rep(pcts/sum(pcts), ncol(bulk_filt)), nrow = ncol(bulk_filt),
                        byrow = TRUE)
    }

    rownames(est_pct) <- colnames(bulk_filt)
    colnames(est_pct) <- colnames(sig_filt)

    gof_list <- gof(bulk_filt, est_pct, sig_filt)

    param_name <- str_glue("{rand_type}_{dataset}_{bulk_dataset}_{granularity}")

    errs_gof_mean[[param_name]] <- gof_list[["means"]]
    errs_gof_subject[[param_name]] <- gof_list[["by_subject"]]

    if (bulk_dataset == "ROSMAP") {
      est_pct_tmp <- est_pct

      # Some datasets are missing one or more of these vascular types
      for (col in c("Endo", "Peri", "VLMC")) {
        if (!(col %in% colnames(est_pct))) {
          est_pct_tmp <- cbind(est_pct_tmp, rep(0, nrow(est_pct_tmp)))
          colnames(est_pct_tmp)[ncol(est_pct_tmp)] <- col
        }
      }

      est_pct2 <- as.data.frame(est_pct_tmp) %>%
        mutate(Neuro = Exc + Inh,
               Oligo = Oligo + OPC,
               Endo = Endo + Peri + VLMC) %>%
        select(colnames(ihc_props))
      est_pct2 <- as.matrix(est_pct2[rownames(specimen2projid),])

      errs_ihc_props[[param_name]] <- data.frame("cor" = sapply(rownames(est_pct2), function(row) {
        projid <- as.character(specimen2projid[row,"donor"])
        cor(est_pct2[row,], ihc_props[projid,],
            use = "na.or.complete")
      }),
      "rMSE" = sapply(rownames(est_pct2), function(row) {
        projid <- as.character(specimen2projid[row,"donor"])
        rmse(est_pct2[row,], ihc_props[projid,])
      }),
      "mAPE" = sapply(rownames(est_pct2), function(row) {
        projid <- as.character(specimen2projid[row,"donor"])
        ok <- est_pct2[row,] != 0 & ihc_props[projid,] != 0
        smape(ihc_props[projid,ok], est_pct2[row,ok]) / 2
      }))

      errs_ihc_pct[[param_name]] <- data.frame("cor" = sapply(rownames(est_pct2), function(row) {
        projid <- as.character(specimen2projid[row,"donor"])
        cor(est_pct2[row,], ihc_pct[projid,],
            use = "na.or.complete")
      }),
      "rMSE" = sapply(rownames(est_pct2), function(row) {
        projid <- as.character(specimen2projid[row,"donor"])
        rmse(est_pct2[row,], ihc_pct[projid,])
      }),
      "mAPE" = sapply(rownames(est_pct2), function(row) {
        projid <- as.character(specimen2projid[row,"donor"])
        ok <- est_pct2[row,] != 0 & ihc_pct[projid,] != 0
        smape(ihc_pct[projid,ok], est_pct2[row,ok]) / 2
      }))
    }

    params[[param_name]] <- data.frame(reference_data_name = dataset,
                                       test_data_name = bulk_dataset,
                                       granularity = granularity,
                                       random_type = rand_type)
  }

  err_list <- list("gof_mean" = do.call(rbind, errs_gof_mean),
                   "gof_subject" = errs_gof_subject,
                   "errs_ihc_props" = errs_ihc_props,
                   "errs_ihc_pct" = errs_ihc_pct,
                   "params" = params)

  Save_ErrorList(err_list, "random", dataset, bulk_dataset, granularity)
  print("Done")
}

