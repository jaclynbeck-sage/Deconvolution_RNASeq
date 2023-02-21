library(dtangleSparse) # dtangle with my mods to make it sparse matrix-friendly
library(Matrix)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(stringr)
library(dplyr)

source(file.path("functions",  "FileIO_HelperFunctions.R"))

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito")#,
              #"seaRef") #, "seaAD")

datatypes <- c("donors", "training")

params <- expand.grid(dataset = datasets,
                      granularity = c("broad"), #, "fine"),
                      input_type = c("singlecell", "pseudobulk"),
                      stringsAsFactors = FALSE) %>% arrange(dataset)

for (P in 1:nrow(params)) {
  dataset <- params$dataset[P]
  granularity <- params$granularity[P]
  input_type <- params$input_type[P]

  if (input_type == "singlecell") {
    input_obj <- Load_SingleCell(dataset, granularity, output_type = "logcpm")
  }
  else { # Input is pseudobulk pure samples
    input_obj <- Load_PseudobulkPureSamples(dataset, granularity,
                                            output_type = "logcpm")
  }

  metadata <- colData(input_obj)
  input_mat <- assay(input_obj, "counts")

  # Clear up as much memory as possible
  rm(input_obj)
  gc()

  celltypes <- levels(metadata$celltype)
  pure_samples <- lapply(celltypes, function(ct) {
    which(metadata$celltype == ct)
  })
  names(pure_samples) <- celltypes

  # TEMPORARY: dtangle code below will not work with a DelayedArray.
  # The seaRef dataset will fit in memory all at once, so this converts it
  # to a sparse matrix. The seaAD data set will NOT fit so this won't work on it.
  if (is(input_mat, "DelayedArray")) {
    input_mat <- as(input_mat, "CsparseMatrix")
  }

  # To save on memory, we load the input once and run on both donor and training
  # data, and on all permutations of marker parameters, rather than splitting
  # this into many separate threads that each have to load the input data.
  for (datatype in datatypes) {
    pseudobulk <- Load_Pseudobulk(dataset, datatype, granularity, "logcpm")
    bulk_mat <- assays(pseudobulk)[["counts"]]

    # Clear up as much memory as possible
    rm(pseudobulk)
    gc()

    # These SHOULD have the same rownames, but just in case.
    keepgene <- intersect(rownames(input_mat), rownames(bulk_mat))

    marker_methods <- c("ratio", "diff")
    if (input_type == "pseudobulk") {
      marker_methods <- c("ratio", "diff", "p.value", "regression")
    }

    params_run <- expand.grid(marker_method = marker_methods,
                              gamma_name = c("auto", 0.5),
                              sum_fn_type = c("mean", "median"),
                              n_markers = c(0.01, 0.02, 0.05, 0.1, 0.2,
                                            0.5, 0.75, 1.0),
                              stringsAsFactors = FALSE)

    # Early testing showed we don't need to test higher percentages for
    # method = "ratio", since this method returns the full set of genes.
    params_run <- subset(params_run, marker_method != "ratio" | n_markers <= 0.2)

    # Pre-combine matrices so this isn't repeatedly done on every dtangle call.
    # Input data must be first so indices in pure_samples are correct.
    Y <- t(cbind(input_mat[keepgene,],
                 bulk_mat[keepgene,]))
    gc()

    dtangle_list <- list()

    for (R in 1:nrow(params_run)) {
      marker_method <- params_run$marker_method[R]
      gamma_name <- params_run$gamma_name[R]
      sum_fn_type <- params_run$sum_fn_type[R]
      n_markers <- params_run$n_markers[R]

      gamma <- NULL
      if (gamma_name != "auto") {
        gamma <- as.numeric(gamma_name)
      }

      sum_fn <- mean
      if (sum_fn_type == "median") {
        sum_fn <- median
      }

      markers <- Load_DtangleMarkers(dataset, granularity, input_type,
                                     marker_method)
      n_markers_name <- n_markers

      # dtangle doesn't interpret "1" as 100%, so we need to input a list of
      # the length of each marker set instead
      if (n_markers == 1) {
        n_markers_name <- "all"
        n_markers <- lengths(markers$L)
      }

      name <- str_glue(paste0("{dataset}_{granularity}_{datatype}_",
                              "{input_type}_{R}"))

      result <- dtangle(Y = Y,
                        pure_samples = pure_samples,
                        data_type = "rna-seq",
                        gamma = gamma, # If gamma is not NULL, it will override data_type argument
                        n_markers = n_markers,
                        markers = markers,
                        summary_fn = sum_fn)

      # Only keep results for pseudobulk samples
      result$estimates <- result$estimates[colnames(bulk_mat), ]

      # Add the params we used to generate this run
      result$params <- cbind(params[P,], "datatype" = datatype, params_run[R,])

      dtangle_list[[name]] <- result

      gc()
      print(paste(result$params, collapse = "  "))
    } # End params_run loop

    # Next iteration will start with new data, remove the old data
    rm(Y)
    gc()

    # Save the completed list
    print("Saving final list...")
    Save_DtangleOutputList(dtangle_list, dataset, datatype, granularity, input_type)
  } # end datatypes loop
}
