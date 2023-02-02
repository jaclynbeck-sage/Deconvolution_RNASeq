CreatePseudobulk_ByDonor <- function(singlecell_counts, metadata, dataset, dir_pseudobulk) {
  # This is SIGNIFICANTLY faster than calling rowSums on individual donor sets

  y <- model.matrix(~0 + donor, data = metadata)
  counts <- singlecell_counts %*% y
  counts <- as(counts, "CsparseMatrix") # For cases where counts is a DelayedMatrix

  # colnames end up as "donor<#>" because of model.matrix. Remove the "donor".
  colnames(counts) <- str_replace(colnames(counts), "donor", "")

  propCells <- table(metadata$donor, metadata$broadcelltype)
  propCells <- propCells / rowSums(propCells)

  pctRNA <- CalculatePercentRNA(singlecell_counts, metadata$donor, metadata$broadcelltype)

  pseudobulk <- SummarizedExperiment(assays = SimpleList(counts = counts),
                                     metadata = list("propCells" = propCells,
                                                     "pctRNA" = pctRNA))

  saveRDS(pseudobulk, file = file.path(dir_pseudobulk,
                                       paste0("pseudobulk_", dataset,
                                              "_donors_broadcelltypes.rds")))

  # TODO Is there a smooth way to put both broad and fine cell types in one object
  # without requiring a lot of reshaping/processing when it's read in from a file?
  propCells_fine <- table(metadata$donor, metadata$subcluster)
  propCells_fine <- propCells_fine / rowSums(propCells_fine)

  pctRNA_fine <- CalculatePercentRNA(singlecell_counts, metadata$donor, metadata$subcluster)

  pseudobulk_fine <- SummarizedExperiment(assays = SimpleList(counts = counts),
                                          metadata = list("propCells" = propCells_fine,
                                                          "pctRNA" = pctRNA_fine))

  saveRDS(pseudobulk_fine, file = file.path(dir_pseudobulk,
                                            paste0("pseudobulk_", dataset,
                                                   "_donors_finecelltypes.rds")))

}
