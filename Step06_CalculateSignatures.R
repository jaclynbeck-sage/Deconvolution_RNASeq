# This script calculates cell-type gene signatures, using the broad and sub
# class cell type assignments output from mapping in Step 04. It also calculates
# the average cell size (average counts per cell per cell type, normalized) for
# use with error comparison.
library(SingleCellExperiment)
library(omnideconv)
library(parallel)

source(file.path("functions", "General_HelperFunctions.R"))
source(file.path("functions", "CibersortX_HelperFunctions.R"))

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

    # Column and rownames of sce get changed to Cibersort-compatible values in
    # the save function
    f_name <- Save_SingleCellToCibersort(sce, sig_dir)
    celltype_anno <- Cells_To_Cibersort(sce$celltype) # Make sure this field is compatible
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

    # CibersortX doesn't subtract the "1" pseudocount after projecting back to
    # linear space, so we do that here
    if (min(sig) == 1) {
      sig <- sig - 1
    }

    # Cleanup finished docker container
    Cleanup_Cibersort_Docker()

    return(sig)
  })

  names(cx_signatures) <- c("broad_class", "sub_class")
  signatures[["cibersortx"]] <- cx_signatures

  saveRDS(signatures,
          file = file.path(dir_signatures, str_glue("{dataset}_signature.rds")))
}


# Calculate batch-corrected signatures for CibersortX --------------------------

# This is *extremely* time- and memory-intensive. Takes over a day to run with 8
# cores and 64 GB of RAM.
#
# CibersortX treats the data like a dense matrix, which takes a lot of time and
# memory to load and process. To avoid this, we create a smaller data set by
# combining multiple cells of the same type into "metacells" to get about 10K
# cells total, maintaining the overall proportion of different cell types. The
# result has >0.99 correlation with the results from running on the full
# single-cell matrix, which is good given that there is a random element to
# Cibersort's algorithm and the outputs are never exactly the same.
if (run_cibersort_batch_correct) {
  cores <- round(parallel::detectCores() / 4) # Assumes 8 GB RAM per core
  cluster_type <- "FORK"
  cluster_outfile <- "cibersort_batch_correct.txt"

  singlecell_datasets <- all_singlecell_datasets()

  params_bulk <- expand.grid(
    test_data_name = setdiff(all_bulk_datasets(), "Mayo_CBE"), # Skip Mayo CBE
    granularity = c("broad_class", "sub_class"),
    normalization = c("cpm", "tpm"),
    regression_method = c("none", "edger", "lme", "combat"),
    stringsAsFactors = FALSE
  )

  for (dataset in singlecell_datasets) {
    sig_dir <- file.path(dir_cibersort, str_glue("{dataset}_batch"))
    dir.create(sig_dir, showWarnings = FALSE, mode = "777")

    # See if broad and sub class files are already there from a previous run. If
    # they are, we don't need to re-aggregate the data, as it will be identical.
    metacell_files <- list(broad_class = c(), sub_class = c())
    metacell_files$broad_class <- list.files(sig_dir,
                                             pattern = "broad_class_sample",
                                             full.names = TRUE)
    metacell_files$sub_class <- list.files(sig_dir,
                                           pattern = "sub_class_sample",
                                           full.names = TRUE)

    # Otherwise, aggregate cells into metacells based on sub_class
    if (!all(lengths(metacell_files) == 1)) {
      set.seed(sageRNAUtils::string_to_seed(str_glue("{dataset}_batch")))

      # CibersortX converts everything to CPM and adds re-sampled cells together
      # after CPM conversion, so we will use CPM too
      sce <- Load_SingleCell(dataset, "sub_class", "cpm")

      # How many of each cell type we need to end up with to get ~10K cells
      fracs <- table(sce$celltype)
      fracs <- round(fracs / sum(fracs) * 10000)

      # How many cells to combine per metacell
      batch_size <- ceiling(ncol(sce) / sum(fracs))
      print(paste(dataset, "batch size:", batch_size))

      agg <- lapply(names(fracs), function(ct) {
        print(paste0(ct, ": ", fracs[ct], " metacells"))
        sce_sub <- sce[, sce$celltype == ct]

        batches <- rep(1:fracs[ct], times = batch_size)
        batches <- batches[1:ncol(sce_sub)]

        # Shuffle the cells before grouping them just in case there's some pattern
        # in the original order
        shuffled <- sample(colnames(sce_sub), ncol(sce_sub), replace = FALSE)
        colData(sce_sub)[shuffled, "batch"] <- batches

        metacells <- scuttle::aggregateAcrossCells(sce_sub, ids = sce_sub$batch,
                                                   statistics = "sum")
        colnames(metacells) <- as.character(metacells$celltype)
        counts(metacells) <- scuttle::calculateCPM(metacells) # Converts to matrix
        return(metacells)
      })

      agg <- do.call(cbind, agg) # Combine everything into a single SCE object

      # This file will be for sub_class resolution
      metacell_f_name <- Save_SingleCellToCibersort(agg, sig_dir)
      sub_filename <- file.path(sig_dir, paste0("sub_class_", basename(metacell_f_name)))
      file.rename(metacell_f_name, sub_filename)

      metacell_files$sub_class <- sub_filename

      # Change the column names to broad class cell types
      agg$celltype <- agg$broad_class
      colnames(agg) <- agg$broad_class

      metacell_f_name <- Save_SingleCellToCibersort(agg, sig_dir)
      broad_filename <- file.path(sig_dir, paste0("broad_class_", basename(metacell_f_name)))
      file.rename(metacell_f_name, broad_filename)

      metacell_files$broad_class <- broad_filename

      rm(sce, agg)
      gc()
    }

    # Loop over each bulk data set, using the same set of 10K metacells as input
    # for each one
    cl <- makeCluster(cores, type = cluster_type, outfile = cluster_outfile)

    parLapply(cl, 1:nrow(params_bulk), function(B) {
      params <- cbind("reference_data_name" = dataset, params_bulk[B, ])

      file_id <- paste(params, collapse="_")
      print(file_id)

      # See if the signature file already exists
      adj_sig_file <- list.files(dir_cibersort_corrected_signatures, pattern = file_id)

      # File already exists, no need to re-run
      if (length(adj_sig_file) == 1) {
        message(str_glue("Found signature file for {file_id}. Skipping..."))
        return(NULL)
      }

      out_dir <- file.path(sig_dir, file_id)
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      metacell_f_name <- metacell_files[[params$granularity]]

      # omnideconv expects the file to be named this way
      metacell_file_copy <- file.path(out_dir, "sample_file_for_cibersort.tsv")
      file.copy(metacell_f_name, metacell_file_copy)

      # Read bulk data
      bulk <- Load_BulkData(params$test_data_name,
                            normalization = params$normalization,
                            regression_method = params$regression_method)
      bulk <- as.matrix(assay(bulk, "counts"))
      bulk <- Dimnames_To_Cibersort(bulk)

      # Read the default signature for this granularity
      signature <- Load_SignatureMatrix(dataset, params$granularity, "cpm")
      signature <- Dimnames_To_Cibersort(signature, TRUE)

      # We need to run the CibersortX docker container but need to stop waiting
      # for it to finish after about 10 minutes. We do this because it takes
      # between 10-30 minutes to generate the signature matrix file but can take
      # several hours afterward to deconvolute the bulk data, and we don't need
      # the deconvoluted data and don't need to wait for that long. However,
      # omnideconv does not have a setting for calling the docker command with a
      # timeout and R.utils::withTimeout() can't interrupt system() commands, so
      # I have replicated some of the omnideconv::deconvolute_cibersort()
      # function here in order to call the docker command with a 10 minute
      # timeout. Then, this script polls the output folder once per minute until
      # it finds a completed signature matrix file.

      sig_file <- file.path(out_dir, "signature_matrix.txt")
      readr::write_tsv(data.frame(NAME = rownames(signature), signature),
                       sig_file)

      bulk_file <- file.path(out_dir, "mixture_file_for_cibersort.txt")
      readr::write_tsv(data.frame(Gene = rownames(bulk), bulk),
                       bulk_file)

      cmd <- Cibersort_Batch_Correct_Command(file_id, out_dir, bulk_file,
                                             sig_file, metacell_file_copy)

      # Start the docker container, wait for 10 minutes, and return to R. The
      # container will keep running but this script can continue polling below.
      code <- system(cmd, timeout = 600)
      print(code)

      # Poll for the signature matrix file and wait until it exists. Loop until
      # we find the file, or until we've been checking for over 2 hours (120
      # times x 1 minute each poll)
      file_found <- FALSE
      n_iter <- 0

      while (!file_found & n_iter < 120) {
        adj_sig_file <- list.files(out_dir,
                                   pattern = str_glue("{file_id}.*sigmatrix_Adjusted.txt"),
                                   full.names = TRUE)
        if (length(adj_sig_file) != 0) {
          file_found <- TRUE
        }

        n_iter <- n_iter + 1
        print(n_iter)

        # Poll once a minute. Even if we found the file, we still wait another
        # minute just to make sure the whole file is written out to disk
        Sys.sleep(60)
      }

      if (length(adj_sig_file) != 0) {
        new_filename <- file.path(dir_cibersort_corrected_signatures,
                                  str_glue("CIBERSORTx_{file_id}_sigmatrix_Adjusted.txt"))

        # We are keeping the adjusted signature exactly as output by CibersortX,
        # without fixing the column or row names. Otherwise we'd just have to
        # convert them again when we run the full algorithm.
        file.copy(adj_sig_file[1], new_filename)

      } else {
        message(str_glue("No signature matrix found for {file_id} in {out_dir}"))
      }

      # Remove the copy of the single cell reference file and the signature/bulk
      # files
      file.remove(metacell_file_copy)
      file.remove(sig_file)
      file.remove(bulk_file)

      # The docker container is likely still running on its own but we try
      # cleaning up anyway. Since there's no way to tell which docker
      # container was spawned by this process, we stop and remove any
      # container that has been running for over 2 hours, which is the limit
      # we've set above for polling for the signature file. If it's been
      # shorter than that, some future iteration of this loop will catch the
      # container once it's been long enough.
      Cleanup_Cibersort_Docker(timeout = 2, timeout_units = "hours")
    })

    stopCluster(cl)
  }
}
