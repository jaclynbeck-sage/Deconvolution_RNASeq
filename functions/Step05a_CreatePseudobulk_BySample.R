# This function creates one pseudobulk sample per sample in the single cell
# data set, and creates the necessary metadata for the proportion of cells and
# percent RNA in each sample for both broad and fine cell types. The pseudobulk
# sets are then written to files as SummarizedExperiment objects.
#
# This script sacrifices readability for speed, as matrix multiplication is
# >10x faster than using rowSums() at this scale. Below is how things work:
#
# 1. We create the model matrix "y" where rows are cells and columns are samples.
#    Entries in the matrix are 1 if that cell belongs to that sample, 0 otherwise:
#            Sample1 Sample2
#     Cell1  1       0
#     Cell2  1       0
#     Cell3  0       1
#
# 2. Multiplying this y matrix with the counts matrix adds the counts of all
#    cells in a sample to create the pseudobulk sample:
#        Cell1 Cell2 Cell3            Sample1 Sample2           Sample1 Sample2
# Gene1  0     1     1         Cell1  1       0          Gene1  1       1
# Gene2  1     10    5     x   Cell2  1       0      =   Gene2  11      5
# Gene3  10    5     0         Cell3  0       1          Gene3  15      0
# Gene4  2     0     2                                   Gene4  2       2
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

source("Filenames.R")
source(file.path("functions", "FileIO_HelperFunctions.R"))
source(file.path("functions", "General_HelperFunctions.R"))

CreatePseudobulk_BySample <- function(singlecell_counts, metadata, dataset) {

  y <- model.matrix(~0 + sample, data = metadata)
  counts <- singlecell_counts %*% y

  # colnames end up as "sample<#>" because of model.matrix. Remove the "sample".
  colnames(counts) <- str_replace(colnames(counts), "sample", "")
  counts <- as(counts, "matrix")

  pb_meta <- as.data.frame(metadata) %>% select(sample, diagnosis) %>% distinct()
  pb_meta$diagnosis <- factor(pb_meta$diagnosis)
  rownames(pb_meta) <- pb_meta$sample
  pb_meta <- pb_meta[colnames(counts),]
  pb_meta$tmm_factors <- calcNormFactors(counts, method = "TMMwsp")

  propCells_broad <- table(metadata$sample, metadata$broadcelltype)
  propCells_broad <- sweep(propCells_broad, 1, rowSums(propCells_broad), "/")

  pctRNA_broad <- CalculatePercentRNA(singlecell_counts, metadata$sample,
                                      metadata$broadcelltype)

  pseudobulk <- SummarizedExperiment(assays = SimpleList(counts = counts),
                                     colData = pb_meta,
                                     metadata = list("propCells" = propCells_broad,
                                                     "pctRNA" = pctRNA_broad))
  Save_Pseudobulk(pseudobulk, dataset, "sc_samples", "broad_class")

  # The counts for the fine cell types pseudobulk set are the same, only the
  # metadata changes. But we create it as a separate file to make looping
  # easier further down the pipeline.
  propCells_fine <- table(metadata$sample, metadata$subcluster)
  propCells_fine <- sweep(propCells_fine, 1, rowSums(propCells_fine), "/")

  pctRNA_fine <- CalculatePercentRNA(singlecell_counts, metadata$sample,
                                     metadata$subcluster)

  pseudobulk_fine <- SummarizedExperiment(assays = SimpleList(counts = counts),
                                          colData = pb_meta,
                                          metadata = list("propCells" = propCells_fine,
                                                          "pctRNA" = pctRNA_fine))
  Save_Pseudobulk(pseudobulk_fine, dataset, "sc_samples", "sub_class")
}
