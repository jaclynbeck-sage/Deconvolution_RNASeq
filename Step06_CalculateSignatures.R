# This script calculates cell-type gene signatures, using the broad and sub
# class cell type assignments output from mapping in Step 04. It also calculates
# the average cell size (average counts per cell per cell type, normalized) for
# use with error comparison.
library(SingleCellExperiment)
library(omnideconv)

source(file.path("functions", "General_HelperFunctions.R"))

datasets <- all_singlecell_datasets()

run_cibersort_signatures <- TRUE
run_cibersort_batch_correct <- TRUE

for (dataset in datasets) {
  # Calculate the "A" matrix that is needed to convert propCells to pctRNA
  A_broad <- CalculateA(dataset, "broad_class")
  A_sub <- CalculateA(dataset, "sub_class")

  saveRDS(list("A_broad_class" = A_broad, "A_sub_class" = A_sub),
          file = file.path(dir_signatures, str_glue("{dataset}_A_matrix.rds")))

  # Calculate a signature for each cell type. This matrix includes all genes in
  # the data set and isn't filtered at this point.
  signatures <- lapply(c("cpm", "tmm"), function(normalization) {
    pb <- Load_PseudobulkPureSamples(dataset, "broad_class", normalization)
    sig_broad <- scuttle::summarizeAssayByGroup(pb, ids = pb$celltype, statistics = "mean")

    pb <- Load_PseudobulkPureSamples(dataset, "sub_class", normalization)
    sig_sub <- scuttle::summarizeAssayByGroup(pb, ids = pb$celltype, statistics = "mean")

    return(list("broad_class" = assay(sig_broad, "mean"),
                "sub_class" = assay(sig_sub, "mean")))
  })

  names(signatures) <- c("cpm", "tmm")

  # Optionally skip calling CibersortX to have it calculate signature matrices
  if (!run_cibersort_signatures) {
    saveRDS(signatures,
            file = file.path(dir_signatures, str_glue("{dataset}_signature.rds")))
    next
  }

  # Calculate the CibersortX signature for each cell type
  # Very memory and disk space intensive
  cx_signatures <- lapply(c("broad_class", "sub_class"), function(granularity) {
    # omnideconv will do this automatically but doesn't clear out the memory
    # after casting to a dense matrix and doesn't expose their write function,
    # so this is a re-implementation.
    sce <- Load_SingleCell(dataset, granularity, "counts")

    sig_dir <- file.path(dir_cibersort, str_glue("{dataset}_{granularity}"))
    dir.create(sig_dir, showWarnings = FALSE, mode = "777")

    f_name <- Save_SingleCellToCibersort(sce, sig_dir)
    celltype_anno <- Cells_To_Cibersort(sce$celltype)
    rm(sce)
    gc()

    # Assumes that environment variables CIBERSORT_EMAIL and CIBERSORT_TOKEN
    # exist. Works on Linux/Mac, not tested on Windows.
    set_cibersortx_credentials(email = Sys.getenv("CIBERSORT_EMAIL"),
                               token = Sys.getenv("CIBERSORT_TOKEN"))

    sig <- build_model_cibersortx(f_name,
                                  cell_type_annotations = celltype_anno,
                                  container = "docker",
                                  input_dir = sig_dir,
                                  output_dir = sig_dir,
                                  verbose = TRUE)

    # Fix cell type names and gene names, which were modified to work with
    # CibersortX. The signature comes out with genes sorted alphabetically so we
    # maintain that ordering.
    sig <- Cibersort_Celltypes_To_Default(sig, colnames(signatures$cpm[[granularity]]))
    sig <- Cibersort_Genes_To_Default(sig, sort(rownames(signatures$cpm[[granularity]])))

    # Cleanup finished docker container
    Cleanup_Cibersort_Docker()

    return(sig)
  })

  names(cx_signatures) <- c("broad_class", "sub_class")
  signatures[["cibersortx"]] <- cx_signatures

  saveRDS(signatures,
          file = file.path(dir_signatures, str_glue("{dataset}_signature.rds")))
}

if (run_cibersort_batch_correct) {
  # Calculate batch-corrected signatures for CibersortX -- *extremely* time- and
  # memory-intensive

  # Only use more than one core if there is a lot of memory available. Assumes
  # 8 GB RAM per CPU
  #n_cores <- min(parallel::detectCores() / 8, 1) # TODO multi-thread

  params_singlecell <- expand.grid(
    reference_data_name = all_singlecell_datasets(),
    granularity = c("broad_class", "sub_class"),
    stringsAsFactors = FALSE
  )

  params_bulk <- expand.grid(
    test_data_name = setdiff(all_bulk_datasets(), "Mayo_CBE"), # Skip Mayo CBE
    normalization = c("cpm", "tpm"),
    regression_method = c("none", "edger", "lme", "combat"),
    stringsAsFactors = FALSE
  )

  set_cibersortx_credentials(email = Sys.getenv("CIBERSORT_EMAIL"),
                             token = Sys.getenv("CIBERSORT_TOKEN"))

  for (S in 1:nrow(params_singlecell)) {
    dataset <- params_singlecell$reference_data_name[S]
    granularity <- params_singlecell$granularity[S]

    # doesn't matter if it's CPM or counts so we just use CPM
    sce <- Load_SingleCell(dataset, granularity, "cpm")

    sig_dir <- file.path(dir_cibersort, str_glue("{dataset}_{granularity}_batch"))
    dir.create(sig_dir, showWarnings = FALSE, mode = "777")

    f_name <- Save_SingleCellToCibersort(sce, sig_dir)

    signature <- Load_SignatureMatrix(dataset, granularity, "cpm")

    # Save original row and column names for replacement afterward
    orig_celltypes <- colnames(signature)
    orig_genes <- rownames(signature)

    signature <- Dimnames_To_Cibersort(signature)

    for (B in 1:nrow(params_bulk)) {
      params <- cbind(params_singlecell[S,], params_bulk[B, ])
      file_id <- paste(params, collapse="_")

      out_dir <- file.path(sig_dir, file_id)
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      new_f_name <- file.path(out_dir, basename(f_name))
      file.copy(f_name, new_f_name)

      bulk <- Load_BulkData(params$test_data_name,
                            normalization = params$normalization,
                            regression_method = params$regression_method)
      bulk <- as.matrix(assay(bulk, "counts"))
      bulk <- Dimnames_To_Cibersort(bulk)

      out <- omnideconv::deconvolute_cibersortx(
        bulk, as.matrix(signature),
        single_cell_object = basename(new_f_name), # filename or NULL
        cell_type_annotations = colnames(signature), # dummy variable, not used but has to have a non-null value
        rmbatch_S_mode = TRUE,
        verbose = TRUE,
        container = "docker",
        input_dir = out_dir,
        output_dir = out_dir,
        qn = FALSE,
        absolute = FALSE,
        label = file_id
      )

      sig_file <- list.files(out_dir,
                             pattern = str_glue("{file_id}.*sigmatrix_Adjusted.txt"),
                             full.names = TRUE)

      if (length(sig_file) != 0) {
        new_filename <- file.path(dir_cibersort_corrected_signatures,
                                  str_glue("CIBERSORTx_{file_id}_sigmatrix_Adjusted.txt"))

        file.copy(sig_file[1], new_filename)

        # Cleanup finished docker container
        Cleanup_Cibersort_Docker()
      } else {
        message(str_glue("No signature matrix found for {file_id} in {out_dir}"))
      }
    }
  }
}
