library(dtangleSparse) # dtangle with my mods to make it sparse matrix-friendly
library(Matrix)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(scuttle)
library(stringr)
library(tidyr)

source("Filenames.R")

marker_file_format <- paste0("dtangle_markers_{sndata}_{cellclasstype}",
                             "_normalization_{normtype}_method_{marker_meth}",
                             "_gamma_{gamma_name}.rds")

list_name_format <- paste0("{sndata}_{cellclasstype}_method_{marker_meth}",
                           "_gamma_{gamma_name}_summaryfn_{sum_fn_type}",
                           "_nmarkers_{n_markers_name}_normalization_{normtype}",
                           "_ROSMAP")

cellclasstype <- "broad" ###either "fine" or "broad"

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito",
              "seaRef") #, "seaAD")

###load bulk and snRNA-seq data###
for (sndata in datasets) {
  sce <- readRDS(file.path(dir_input, paste(sndata, "sce.rds", sep = "_")))
  meta <- colData(sce)

  # Convert gene names to Ensembl IDs
  genes <- rowData(sce)
  rownames(sce) <- genes[rownames(sce), "Ensembl.ID"]

  sce_cpm <- calculateCPM(counts(sce))

  sc_mats <- list()

  sc_mats[["logcpm"]] <- sce_cpm

  if (is(sce_cpm, "DelayedArray") | is(sce_cpm, "matrix")) {
    sc_mats[["logcpm"]] <- log2(sce_cpm + 1)
  }
  else { # Sparse matrix
    sc_mats[["logcpm"]]@x <- log2(sce_cpm@x + 1)
  }

  if (is(sc_mats[["logcpm"]], "DelayedArray")) {
    sc_mats[["logcpm"]] <- as(sc_mats[["logcpm"]], "dgCMatrix")
  }

  if (cellclasstype == "broad") {
    broadtypes <- levels(meta$broadcelltype)
    pure_samples <- list()
    for (ii in broadtypes) {
      pure_samples[[ii]] <- which(meta$broadcelltype == ii)
    }
  }

  if (cellclasstype == "fine") {
    finetypes <- levels(meta$subcluster)
    pure_samples <- list()
    for (ii in finetypes) {
      pure_samples[[ii]] <- which(meta$subcluster == ii)
    }
  }

  # Specific to ROSMAP for now
  bulk <- read.table(file_rosmap, header = TRUE, row.names = 1)
  bulk_cpm <- calculateCPM(bulk)

  bulk_mats <- list()
  bulk_mats[["logcpm"]] <- log2(bulk_cpm + 1)

  keepgene <- intersect(rownames(sce), rownames(bulk))

  # Clear up as much memory as possible
  rm(meta, sce, sce_cpm, bulk, bulk_cpm)
  gc()

  best_params <- readRDS(file.path(dir_output, paste0("best_params_", sndata,
                                                      "_", cellclasstype, ".rds")))
  best_params <- subset(best_params, algorithm == "dtangle")
  params <- best_params %>%
    extract(name, c("marker_meth", "gamma_name", "sum_fn_type", "n_markers", "normtype"),
            "method_([a-z]+)_gamma_([[:alnum:]]+)_summaryfn_([a-z]+)_nmarkers_([0-9\\.|all]+)_normalization_([[:alnum:]]+)")

  dtangle_list <- list()

  for (R in 1:nrow(params)) {
    marker_meth <- params$marker_meth[R]
    gamma_name <- params$gamma_name[R]
    sum_fn_type <- params$sum_fn_type[R]
    n_markers <- as.numeric(params$n_markers[R])
    normtype <- params$normtype[R]

    n_markers_name <- n_markers

    if (gamma_name == "auto") {
      gamma <- NULL
    }
    else {
      gamma <- as.numeric(gamma_name)
    }

    sum_fn <- mean
    if (sum_fn_type == "median") {
      sum_fn <- median
    }

    markers <- readRDS(file.path(dir_markers, str_glue(marker_file_format)))
    markers <- markers$L

    # Convert marker gene names to Ensembl IDs.
    # Remove markers that aren't in both data sets. The original marker list
    # for each cell is a named vector of indices corresponding to rownames,
    # where the names of the vector are the genes. Since we've removed genes
    # from the matrices, the indicies are no longer valid so these get removed
    # and replaced with a vector of gene names, which works just as well.
    for (cell in names(markers)) {
      old <- names(markers[[cell]])
      new <- genes[old, "Ensembl.ID"]
      markers[[cell]] <- intersect(keepgene, new)
    }

    Y <- t(cbind(sc_mats[[normtype]][keepgene,],
                 bulk_mats[[normtype]][keepgene,]))
    gc()

    name <- str_glue(list_name_format)

    result <- dtangle(Y = Y,
                      pure_samples = pure_samples,
                      data_type = "rna-seq",
                      gamma = gamma, # If gamma is not NULL, it will override data_type argument
                      n_markers = n_markers,
                      markers = markers,
                      summary_fn = sum_fn)

    # Only keep results for pseudobulk samples
    result$estimates <- result$estimates[colnames(bulk_mats[[normtype]]), ]

    dtangle_list[[name]] <- result

    gc()
    print(name)

    # Next iteration will start with new data, remove the old data
    rm(Y)
    gc()
  } # End params loop

  print("Saving final list...")
  saveRDS(dtangle_list, file = file.path(dir_output,
                                         paste0("dtangle_list_", sndata,
                                                "_", cellclasstype, "_ROSMAP.rds")))
}
