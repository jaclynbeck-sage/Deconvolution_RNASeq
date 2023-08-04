library(Matrix)
library(SingleCellExperiment)
library(Metrics)
library(dplyr)
library(scuttle)
library(stringr)

source(file.path("functions", "General_HelperFunctions.R"))
source(file.path("functions", "Error_HelperFunctions.R"))

datasets <- c("cain", "lau", "leng", "mathys", "morabito", "seaRef") #, "seaAD")
datasets <- c("mathys")
granularity <- "broad"
bulk_dataset <- "Mayo"

est_fields = list("dtangle" = "estimates",
                  "music" = "Est.pctRNA.weighted",
                  "hspe" = "estimates",
                  "deconRNASeq" = "out.all")

algorithms <- names(est_fields)
algorithms <- c("deconRNASeq")

err_list <- list()

for (dataset in datasets) {
  #sig_matrix <- Load_SignatureMatrix(dataset, granularity)
  #comp_genes <- readRDS(file.path(dir_input, "compartment_specific_genes.rds"))
  #sig_matrix <- sig_matrix[!(rownames(sig_matrix) %in% comp_genes$Symbol),]

  data <- Load_AlgorithmInputData(dataset, bulk_dataset, granularity,
                                  reference_input_type = "signature",
                                  output_type = "cpm") # TODO normalization

  signature <- data$reference
  bulk_cpm <- assay(data$test, "counts")

  var_genes <- apply(bulk_cpm, 1, var)
  highly_variable <- names(sort(var_genes, decreasing = TRUE))[1:5000]

  # TODO is it correct to re-calculate CPM?
  #sig_matrix <- calculateCPM(sig_matrix[keepgene,])
  #bulk_cpm <- calculateCPM(bulk_cpm[keepgene,])

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

  for (algorithm in algorithms) {
    print(str_glue("Calculating errors for {dataset} / {bulk_dataset}: {algorithm}"))
    est_field <- est_fields[[algorithm]]

    deconv_list <- Load_AlgorithmOutputList(algorithm, dataset, bulk_dataset,
                                            granularity, "cpm") # TODO normalization

    # If the file didn't exist, we haven't run the algorithm for this set of
    # params. Skip it.
    if (length(deconv_list) == 0) {
      next
    }

    params <- do.call(rbind, lapply(deconv_list, "[[", "params"))
    params <- params[,c("marker_type", "marker_subtype", "marker_input_type")] %>%
                  distinct()
    params <- subset(params, marker_type != "None")

    common_markers <- c()
    for (row in 1:nrow(params)) {
      mkrs <- Load_Markers(dataset, granularity,
                           marker_type = params$marker_type[row],
                           marker_subtype = params$marker_subtype[row],
                           input_type = params$marker_input_type[row])
      common_markers <- c(common_markers, unlist(mkrs))
    }

    common_markers <- table(common_markers)
    common_markers <- names(common_markers)[common_markers >= nrow(params)/2]
    common_markers <- intersect(common_markers, rownames(signature))

    lm_list <- list()

    # Calculate error for each parameter set's results
    for (param_id in names(deconv_list)) {
      est_pct <- deconv_list[[param_id]][[est_field]]
      est_pct <- est_pct[colnames(bulk_cpm), colnames(signature)]

      if (any(is.na(est_pct))) {
        message(str_glue("Param set {param_id} has NA values. Skipping..."))
        next
      }

      zeros <- colSums(est_pct) == 0
      if (granularity == "broad" & any(zeros)) {
        cts <- paste(colnames(est_pct)[which(zeros)], collapse = ", ")
        msg <- str_glue(paste0("Param set '{param_id}' has all 0 estimates for",
                               " cell type(s) [{cts}]. Skipping..."))
        message(msg)
        next
      }

      #mod <- lm(t(bulk_cpm) ~ 0 + ., data = data.frame(est_pct))
      #fitted <- t(mod$fitted.values)

      #cor_subject <- diag(cor(bulk_cpm[common_markers,], fitted, use = "na.or.complete"))
      #rmse_subject <- sapply(colnames(bulk_cpm[common_markers,]), function(sample) {
      #  rmse(bulk_cpm[common_markers,sample], fitted[,sample])
      #})
      #mape_subject <- sapply(colnames(bulk_cpm[common_markers,]), function(sample) {
      #  ok <- bulk_cpm[common_markers,sample] != 0 | fitted[,sample] != 0
      #  smape(bulk_cpm[common_markers,sample][ok], fitted[ok, sample]) / 2
      #})

      #nn <- lapply(rownames(bulk_cpm), function(M) {
      #  nnls(est_pct, bulk_cpm[M,])
      #})

      #sig_nn <- t(sapply(nn, "[[", "x"))
      #colnames(sig_nn) <- colnames(signature)
      #rownames(sig_nn) <- rownames(bulk_cpm)

      #fitted <- t(do.call(cbind, lapply(nn, "[[", "fitted")))
      #rownames(fitted) <- rownames(bulk_cpm)

      #cor_subject <- diag(cor(bulk_cpm, fitted, use = "na.or.complete"))
      #rmse_subject <- sapply(colnames(bulk_cpm), function(sample) {
      #  rmse(bulk_cpm[,sample], fitted[,sample])
      #})
      #mape_subject <- sapply(colnames(bulk_cpm), function(sample) {
      #  ok <- bulk_cpm[,sample] != 0 | fitted[,sample] != 0
      #  smape(bulk_cpm[ok, sample], fitted[ok, sample]) / 2
      #})

      #gof_by_subject <- data.frame(cor = cor_subject,
      #                             rMSE = rmse_subject,
      #                             mAPE = mape_subject,
      #                             param_id = param_id,
      #                             subject = names(cor_subject),
      #                             row.names = names(cor_subject))
      #gof_means <- CalcGOF_Means(gof_by_subject, colData(data$test), param_id)

      #lm_list[[param_id]] <- list("signature_estimate" = sig_nn,
      #                            "gof_by_subject" = gof_by_subject,
      #                            "gof_means" = gof_means)

      # For param sets that test on the whole or partially-filtered signature,
      # pick a reasonable subset of genes instead of calculating error on all of
      # them, since most genes don't provide useful information.
      #markers <- common_markers #deconv_list[[param_id]]$markers
      params <- deconv_list[[param_id]]$params
      #if (params$marker_type == "None") {
      #  markers <- common_markers
      #}
      #else {
      #  markers <- Load_Markers(dataset, granularity,
      #                          marker_type = params$marker_type,
      #                          marker_subtype = params$marker_subtype,
      #                          input_type = params$marker_input_type)
      #  markers <- unlist(markers)
      #  markers <- markers[markers %in% rownames(signature)]
      #}
      markers <- rownames(signature)
      #markers <- highly_variable
      #if (length(markers) > 5000) {
      #  sig_filt <- FilterSignature(signature, 3, dataset, granularity, 1,
      #                              marker_type = "dtangle",
      #                              marker_subtype = "diff",
      #                              marker_input_type = "pseudobulk")
      #  markers <- rownames(sig_filt)
      #}
      sig_filt <- signature[markers,]
      bulk_filt <- bulk_cpm[markers,]

      gof_by_subject <- CalcGOF_BySubject(bulk_filt, est_pct, sig_filt, param_id)
      gof_means <- CalcGOF_Means(gof_by_subject, colData(data$test), param_id)

      deconv_list[[param_id]]$gof_by_subject <- gof_by_subject
      deconv_list[[param_id]]$gof_means_all <- gof_means$all_tissue
      deconv_list[[param_id]]$gof_means_by_tissue <- gof_means$by_tissue

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

        errs_ihc_props[[param_id]] <- data.frame("cor" = sapply(rownames(est_pct2), function(row) {
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

        errs_ihc_pct[[param_id]] <- data.frame("cor" = sapply(rownames(est_pct2), function(row) {
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

      #params[[param_id]] <- deconv_list[[param_id]]$params
      print(str_glue("\t{param_id}"), paste(deconv_list[[param_id]]$params, collapse = " "))
    }

    gof_means_all <- lapply(deconv_list, "[[", "gof_means_all")
    gof_means_tissue <- lapply(deconv_list, "[[", "gof_means_by_tissue")
    params <- lapply(deconv_list, "[[", "params")

    err_list <- list("means" = list("all_tissue" = do.call(rbind, gof_means_all),
                                    "by_tissue" = do.call(rbind, gof_means_tissue)),
                     "params" = do.call(rbind, params),
                     "by_subject" = lapply(deconv_list, "[[", "gof_by_subject"))

    if (bulk_dataset == "ROSMAP") {
      ihc_means_all <- lapply(deconv_list, "[[", "ihc_means_all")
      ihc_means_tissue <- lapply(deconv_list, "[[", "ihc_means_by_tissue")
      #err_list$ihc = list("means") TODO
    }

    Save_ErrorList(err_list, algorithm, dataset, bulk_dataset, granularity, "cpm") # TODO normalization
    print("Done")
  }
}

