# This function creates one pseudobulk sample per cell type per donor in the
# single cell data set, and creates the necessary metadata for the proportion of
# cells and percent RNA in each sample for both broad and fine cell types. The
# pseudobulk sets are then written to files as SummarizedExperiment objects.
#
# This script sacrifices readability for speed, as matrix multiplication is
# >10x faster than using rowSums() at this scale. Below is how things work:
#
# 1. We create the model matrix "y" where rows are cells and columns are
#    cell type + donor combinations.
#    Entries in the matrix are 1 if that cell belongs to that donor/celltype,
#    0 otherwise:
#            Ast_Donor1  Ast_Donor2  Mic_Donor1 ...
#     Cell1  1           0           0
#     Cell2  1           0           0
#     Cell3  0           1           0
#     Cell4  0           0           1
#
# 2. Multiplying this y matrix with the counts matrix adds the counts of all
#    cells for each cell type/donor combination to create the pseudobulk sample:
#        Cell1 Cell2 Cell3 Cell4            Ast_Donor1  Ast_Donor2  Mic_Donor1             Ast_Donor1  Ast_Donor2  Mic_Donor1
# Gene1  0     1     1     0         Cell1  1           0           0               Gene1  1           1           0
# Gene2  1     10    5     10     x  Cell2  1           0           0           =   Gene2  11          5           10
# Gene3  10    5     0     1         Cell3  0           1           0               Gene3  15          0           1
# Gene4  2     0     2     5         Cell4  0           0           1               Gene4  2           2           5
#
# Arguments:
#   singlecell_counts - a gene x cell matrix of counts
#   metadata - a cell x feature dataframe describing the cells. Must contain
#              columns "broadcelltype" and "subcluster", corresponding to
#              the broad and fine cell type assignments, respectively, for each
#              cell
#   dataset - the name of the dataset
#   dir_pseudobulk - the directory to write the pseudobulk files to
#
# Returns: nothing

CreatePseudobulk_PureSamplesByDonor <- function(singlecell_counts, metadata, dataset, dir_pseudobulk) {

  for (cell_class in c("broad", "fine")) {
    if (cell_class == "broad") {
      metadata$celltype <- metadata$broadcelltype
    }
    else { # "fine"
      metadata$celltype <- metadata$subcluster
    }

    metadata$celltypedonor <- paste(metadata$celltype, metadata$donor, sep = "_")
    metadata$celltypedonor <- factor(metadata$celltypedonor)

    y <- model.matrix(~0 + celltypedonor, data = metadata)
    colnames(y) <- str_replace(colnames(y), "celltypedonor", "puresample_")

    counts <- singlecell_counts %*% y
    counts <- as(counts, "CsparseMatrix") # For cases where counts is a DelayedMatrix

    propCells <- table(metadata$celltypedonor, metadata$celltype) # TODO check this works
    rownames(propCells) <- paste0("puresample_", rownames(propCells))
    propCells <- sweep(propCells, 1, rowSums(propCells), "/")

    # Since these are pure samples, pctRNA and propCells = 1 where the cell type
    # matches the pure sample.
    pctRNA <- propCells

    pseudobulk <- SummarizedExperiment(assays = SimpleList(counts = counts),
                                       metadata = list("propCells" = propCells,
                                                       "pctRNA" = pctRNA))

    file_name <- file.path(dir_pseudobulk,
                           str_glue("pseudobulk_{dataset}_puresamplesbydonor_{cell_class}celltypes.rds"))
    saveRDS(pseudobulk, file = file_name)
  }
}
