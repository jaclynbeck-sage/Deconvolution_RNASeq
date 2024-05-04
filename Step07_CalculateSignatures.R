# This script calculates cell-type gene signatures, using the broad and sub
# class cell type assignments output from mapping in Step 04. It also calculates
# the average cell size (average counts per cell per cell type, normalized) for
# use with error comparison.
library(SingleCellExperiment)
library(omnideconv)

source(file.path("functions", "FileIO_HelperFunctions.R"))
source(file.path("functions", "General_HelperFunctions.R"))

datasets <- c("cain", "lau", "leng", "mathys", "seaRef")

for (dataset in datasets) {
  # Calculate the "A" matrix that is needed to convert propCells to pctRNA
  A_broad <- CalculateA(dataset, "broad_class")
  A_sub <- CalculateA(dataset, "sub_class")

  saveRDS(list("A_broad_class" = A_broad, "A_sub_class" = A_sub),
          file = file.path(dir_signatures, str_glue("{dataset}_A_matrix.rds")))

  # Calculate a signature for each cell type. This matrix includes all genes in
  # the data set and isn't filtered at this point.
  signatures <- lapply(c("cpm", "tmm"), function(output_type) {
    sig_broad <- CalculateSignature(dataset, "broad_class", output_type, geom_mean = TRUE)
    sig_sub <- CalculateSignature(dataset, "sub_class", output_type, geom_mean = TRUE)
    return(list("broad_class" = sig_broad, "sub_class" = sig_sub))
  })

  names(signatures) <- c("cpm", "tmm")

  # Calculate the CibersortX signature for each cell type
  # Very memory and disk space intensive
  cx_signatures <- lapply(c("broad_class", "sub_class"), function(granularity) {
    # omnideconv will do this automatically but doesn't clear out the memory
    # after casting to a dense matrix and doesn't expose their write function,
    # so this is a re-implementation.
    sce <- Load_SingleCell(dataset, granularity, "counts")
    f_name <- Save_SingleCellToCibersort(sce, dataset, granularity)
    celltypes <- sce$celltype
    rm(sce)
    gc()

    # Assumes that environment variables CIBERSORT_EMAIL and CIBERSORT_TOKEN
    # exist. Works on Linux/Mac, not tested on Windows.
    set_cibersortx_credentials(email = Sys.getenv("CIBERSORT_EMAIL"),
                               token = Sys.getenv("CIBERSORT_TOKEN"))

    sig <- build_model_cibersortx(f_name,
                                  cell_type_annotations = as.character(celltypes),
                                  container = "docker",
                                  input_dir = dir_cibersort,
                                  output_dir = dir_cibersort,
                                  verbose = TRUE)

    # Cell types are out of order and some genes have had "-" replaced by ".",
    # this undoes that
    sig <- sig[, levels(celltypes)]
    rownames(sig) <- str_replace(rownames(sig), "\\.", "-")

    # Cleanup finished docker container
    system("docker rm $(docker ps -a -q --filter ancestor=cibersortx/fractions)")
    gc()

    return(sig)
  })

  names(cx_signatures) <- c("broad_class", "sub_class")
  signatures[["cibersortx"]] <- cx_signatures

  saveRDS(signatures,
          file = file.path(dir_signatures, str_glue("{dataset}_signature.rds")))
}
