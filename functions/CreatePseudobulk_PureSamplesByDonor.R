CreatePseudobulk_PureSamplesByDonor <- function(singlecell_counts, metadata, dataset, dir_pseudobulk) {
  # This is SIGNIFICANTLY faster than calling rowSums on individual donor sets

  metadata$celltypedonor <- factor(paste(metadata$broadcelltype, metadata$donor, sep = "_"))

  y = model.matrix(~0 + celltypedonor, data = metadata)
  colnames(y) = str_replace(colnames(y), "celltypedonor", "puresample_")

  counts = singlecell_counts %*% y

  props <- table(metadata$celltypedonor, metadata$broadcelltype)
  rownames(props) <- paste0("puresample_", rownames(props))
  # necessary to get the correct shape for Summarized Experiment to convert to DataFrame
  props <- as.data.frame.matrix(props / rowSums(props))

  pseudobulk <- SummarizedExperiment(assays = SimpleList(counts = counts),
                                     colData = props)

  saveRDS(pseudobulk, file = file.path(dir_pseudobulk,
                                       paste0("pseudobulk_", dataset, "_puresamplesbydonor_broadcelltypes.rds")))

  # Fine cell types -- this may not work well for analysis
  metadata$celltypedonor <- factor(paste(metadata$subcluster, metadata$donor, sep = "_"))

  y = model.matrix(~0 + celltypedonor, data = metadata)
  colnames(y) = str_replace(colnames(y), "celltypedonor", "puresample_")

  counts = singlecell_counts %*% y

  props <- table(metadata$celltypedonor, metadata$broadcelltype)
  rownames(props) <- paste0("puresample_", rownames(props))
  # necessary to get the correct shape for Summarized Experiment to convert to DataFrame
  props <- as.data.frame.matrix(props / rowSums(props))

  pseudobulk <- SummarizedExperiment(assays = SimpleList(counts = counts),
                                     colData = props)

  saveRDS(pseudobulk, file = file.path(dir_pseudobulk,
                                       paste0("pseudobulk_", dataset, "_puresamplesbydonor_finecelltypes.rds")))
}
