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
#
# Returns: nothing

source(file.path("functions", "FileIO_HelperFunctions.R"))

CreatePseudobulk_PureSamples <- function(singlecell_counts, metadata, dataset) {

  for (granularity in c("broad", "fine")) {
    if (granularity == "broad") {
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
    counts <- as(counts, "matrix")

    pb_meta <- as.data.frame(metadata) %>%
                  select(donor, diagnosis, celltypedonor) %>% distinct() %>%
                  mutate(donor = paste0("puresample_", celltypedonor),
                         celltype = str_replace(celltypedonor, "_.*", "")) %>%
                  select(-celltypedonor)
    pb_meta$celltype <- factor(pb_meta$celltype)
    pb_meta$diagnosis <- factor(pb_meta$diagnosis)
    rownames(pb_meta) <- pb_meta$donor
    pb_meta <- pb_meta[colnames(counts),]
    pb_meta$tmm_factors <- calcNormFactors(counts, method = "TMMwsp")

    propCells <- table(metadata$celltypedonor, metadata$celltype)
    rownames(propCells) <- paste0("puresample_", rownames(propCells))
    propCells <- sweep(propCells, 1, rowSums(propCells), "/")

    # Since these are pure samples, pctRNA and propCells = 1 where the cell type
    # matches the pure sample.
    pctRNA <- propCells

    pseudobulk <- SummarizedExperiment(assays = SimpleList(counts = counts),
                                       colData = pb_meta,
                                       metadata = list("propCells" = propCells,
                                                       "pctRNA" = pctRNA))

    Save_PseudobulkPureSamples(pseudobulk, dataset, granularity)
  }
}
