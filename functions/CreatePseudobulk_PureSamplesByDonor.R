CreatePseudobulk_PureSamplesByDonor <- function(singlecell_counts, metadata, dataset, dir_pseudobulk) {
  # This is SIGNIFICANTLY faster than calling rowSums on individual donor sets

  metadata$celltypedonor <- factor(paste(metadata$broadcelltype, metadata$donor,
                                         sep = "_"))

  y <- model.matrix(~0 + celltypedonor, data = metadata)
  colnames(y) <- str_replace(colnames(y), "celltypedonor", "puresample_")

  counts <- singlecell_counts %*% y
  counts <- as(counts, "CsparseMatrix") # For cases where counts is a DelayedMatrix

  propCells <- table(metadata$celltypedonor, metadata$broadcelltype)
  rownames(propCells) <- paste0("puresample_", rownames(propCells))
  propCells <- propCells / rowSums(propCells)

  # Since these are pure samples, pctRNA and propCells = 1 where the cell type
  # matches the pure sample.
  pctRNA <- propCells

  pseudobulk <- SummarizedExperiment(assays = SimpleList(counts = counts),
                                     metadata = list("propCells" = propCells,
                                                     "pctRNA" = pctRNA))

  saveRDS(pseudobulk, file = file.path(dir_pseudobulk,
                                       paste0("pseudobulk_", dataset,
                                              "_puresamplesbydonor_broadcelltypes.rds")))

  # Fine cell types -- this may not work well for analysis
  metadata$celltypedonor <- factor(paste(metadata$subcluster, metadata$donor, sep = "_"))

  y <- model.matrix(~0 + celltypedonor, data = metadata)
  colnames(y) <- str_replace(colnames(y), "celltypedonor", "puresample_")

  counts <- singlecell_counts %*% y
  counts <- as(counts, "CsparseMatrix") # For cases where counts is a DelayedMatrix

  propCells <- table(metadata$celltypedonor, metadata$subcluster)
  rownames(propCells) <- paste0("puresample_", rownames(propCells))
  propCells <- propCells / rowSums(propCells)

  pctRNA <- propCells

  pseudobulk <- SummarizedExperiment(assays = SimpleList(counts = counts),
                                     metadata = list("propCells" = propCells,
                                                     "pctRNA" = pctRNA))

  saveRDS(pseudobulk, file = file.path(dir_pseudobulk,
                                       paste0("pseudobulk_", dataset,
                                              "_puresamplesbydonor_finecelltypes.rds")))
}
