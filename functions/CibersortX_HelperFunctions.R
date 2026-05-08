# CibersortX save function -----------------------------------------------------

# Save_SingleCellToCibersort: CibersortX requires passing files to its docker
# container, so the single cell data needs to be written as a dense matrix
# in a specific format. The file is saved to the directory that is shared
# with the CibersortX docker container. omnideconv will do this automatically
# if you pass it a matrix instead of a filename, but it doesn't clear out the
# memory after casting to a dense matrix and doesn't expose its write function,
# so this is a re-implementation.
#
# Arguments:
#   sce - a SingleCellExperiment object, which must have a column for "celltype"
#   output_dir - the path to the directory where the data should be stored
#
# Returns:
#   f_name - the filename of the newly-created file
Save_SingleCellToCibersort <- function(sce, output_dir) {
  # omnideconv requires that the file be named this way
  f_name <- file.path(output_dir, "sample_file_for_cibersort.tsv")

  sce <- Dimnames_To_Cibersort(sce)

  # Header (column names) is "GeneSymbol" followed by the cell type labels for
  # each cell
  cts <- as.data.frame(t(c("GeneSymbol", as.character(sce$celltype))))
  vroom_write(cts, file = f_name, delim = "\t", col_names = FALSE)

  # Write in blocks of 500 genes to avoid loading the whole matrix into memory
  # at once, and to speed this process up considerably.
  block_size <- 500
  n_blocks <- ceiling(nrow(sce) / block_size)

  for (block in 1:n_blocks) {
    start <- (block-1) * block_size + 1
    end <- min(block * block_size, nrow(sce))

    df <- as.data.frame(as.matrix(counts(sce[start:end,])))
    df <- cbind(rownames(df), df)
    vroom_write(df, file = f_name, delim = "\t", col_names = FALSE,
                append = TRUE,
                num_threads = parallel::detectCores()-1)
  }

  rm(df)
  gc()

  return(f_name)
}


# Cibersort re-formatting functions --------------------------------------------

# Replace spaces with "_" and change any other special characters to "."
Genes_To_Cibersort <- function(genes) {
  str_replace_all(as.character(genes), " ENSG", "_ENSG") |> make.names()
}

# Replace "-", " ", and "." with "_", and change any other special characters to
# ".". Optionally, also lower-case the names so sorting is consistent between
# R and C.
Cells_To_Cibersort <- function(celltypes, lowercase = FALSE) {
  res <- str_replace_all(as.character(celltypes), "[-/ \\.]", "_") |> make.names()
  if (lowercase) {
    res <- str_to_lower(res)
  }
  return(res)
}

# Shortcut function to reformat both the gene names and the cell type names for
# CibersortX. Additionally, if a SingleCellExperiment is passed in, this also
# adjusts the "celltype" metadata field. The column names of "obj" don't need to
# be celltypes, they can also be sample names. The reformatting is the same for
# both.
Dimnames_To_Cibersort <- function(obj, col_lower_case = FALSE) {
  # Remove special characters from colnames, remove the space from gene names
  colnames(obj) <- Cells_To_Cibersort(colnames(obj), lowercase = col_lower_case)
  rownames(obj) <- Genes_To_Cibersort(rownames(obj))

  if (is(obj, "SingleCellExperiment")) {
    # Remove special characters from cell type names, and convert them to all
    # lowercase to avoid string sorting issues between R and C.
    obj$celltype <- Cells_To_Cibersort(obj$celltype, lowercase = TRUE)
  }

  return(obj)
}

# Convert CibersortX-formatted cell type names back to their original values.
# Assumes the column names of "obj" are cell types
Cibersort_Celltypes_To_Default <- function(obj, correct_celltype_names) {
  # Cell types come out of CibersortX out of order and lower-case. We put them
  # back in the right order and replace with the correct names.
  cx_names <- Cells_To_Cibersort(correct_celltype_names, lowercase = TRUE)

  stopifnot(all(colnames(obj) %in% cx_names))

  obj <- obj[, cx_names]
  colnames(obj) <- correct_celltype_names
  return(obj)
}

# Sample names get formatted with Cells_To_Cibersort for CibersortX, as they are
# originally the column names of the input bulk data. We replace the sample
# names in the output estimate matrix with the correct names. This function is
# nearly identical to Celltypes_To_Default, however sample names are row names
# in the estimates matrix instead of column names.
Cibersort_Samples_To_Default <- function(obj, correct_sample_names) {
  cx_names <- Cells_To_Cibersort(correct_sample_names, lowercase = FALSE)

  stopifnot(all(rownames(obj) %in% cx_names))

  obj <- obj[cx_names, ]
  rownames(obj) <- correct_sample_names
  return(obj)
}

# Convert CibersortX-formatted gene names back to their original values.
# Assumes the row names of "obj" are genes.
Cibersort_Genes_To_Default <- function(obj, correct_gene_names) {
  cx_genes <- Genes_To_Cibersort(correct_gene_names)

  # There may be fewer than the full number of genes so we have to keep track of
  # the real names of the genes that are kept
  names(cx_genes) <- correct_gene_names
  cx_genes <- cx_genes[cx_genes %in% rownames(obj)]

  stopifnot(length(cx_genes) == nrow(obj))

  obj <- obj[cx_genes, ]
  rownames(obj) <- names(cx_genes)
  return(obj)
}


# Docker container control -----------------------------------------------------

# Find all exited docker containers of type "cibersortx/fractions" and remove
# them with "docker rm". Optionally, this function will also stop any running
# containers that have been up for longer than <timeout>, which then get removed
# afterward. If a timeout is set, "timeout_units" should be "minutes", "hours",
# or "days".
Cleanup_Cibersort_Docker <- function(timeout = NULL, timeout_units = NULL) {
  # If a timeout is set, stop any docker containers that have been running for
  # longer than the timeout
  if (!is.null(timeout) && !is.null(timeout_units)) {
    find_cmd <- paste0(
      "docker ps -a --filter ancestor=cibersortx/fractions --format ",
      "'{{.ID}} {{.Names}} {{.RunningFor}}' | ",
      "awk '/", timeout_units, "/ && $3 >= ", timeout, " { print $1 }'"
    )
    system(paste0("docker stop $(", find_cmd, ")"), ignore.stderr = TRUE)
  }

  # Remove any exited docker containers
  find_cmd <- paste0(
    "docker ps -a -q --filter ancestor=cibersortx/fractions ",
    "--filter status=exited"
  )
  system(paste0("docker rm $(", find_cmd, ")"), ignore.stderr = TRUE)
  gc()
}

Cibersort_Batch_Correct_Command <- function(file_id, out_dir, bulk_filename,
                                            sig_filename, ref_filename) {
  str_glue(
    "docker run -v {out_dir}:/src/data:z -v ",
    "{out_dir}:/src/outdir:z cibersortx/fractions ",
    "--single_cell TRUE --verbose FALSE ",
    "--username {Sys.getenv('CIBERSORT_EMAIL')} ",
    "--token {Sys.getenv('CIBERSORT_TOKEN')} ",
    "--mixture {basename(bulk_filename)} ",
    "--sigmatrix {basename(sig_filename)} ",
    "--perm 0 --label {file_id} --rmbatchBmode FALSE --rmbatchSmode TRUE ",
    "--sourceGEPs {basename(sig_filename)} --QN FALSE --absolute FALSE ",
    "--abs_method sig.score --refsample {basename(ref_filename)}")
}
