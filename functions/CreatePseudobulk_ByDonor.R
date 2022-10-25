CreatePseudobulk_ByDonor <- function(singlecell_counts, metadata, dataset, dir_pseudobulk) {
  # This is SIGNIFICANTLY faster than calling rowSums on individual donor sets

  y = model.matrix(~0 + donor, data = metadata)
  counts = singlecell_counts %*% y

  props <- table(metadata$donor, metadata$broadcelltype)
  rownames(props) <- paste0("donor", rownames(props))
  # necessary to get the correct shape for Summarized Experiment to convert to DataFrame
  props <- as.data.frame.matrix(props / rowSums(props))

  pseudobulk <- SummarizedExperiment(assays = SimpleList(counts = counts),
                                     colData = props)

  saveRDS(pseudobulk, file = file.path(dir_pseudobulk,
                                       paste0("pseudobulk_", dataset, "_bydonor_broadcelltypes.rds")))

  # TODO Is there a smooth way to put both broad and fine cell types in one object
  # without requiring a lot of reshaping/processing when it's read in from a file?
  props_fine <- table(metadata$donor, metadata$subcluster)
  rownames(props_fine) <- paste0("donor", rownames(props_fine))
  props_fine <- as.data.frame.matrix(props_fine / rowSums(props_fine))

  pseudobulk_fine <- SummarizedExperiment(assays = SimpleList(counts = counts),
                                          colData = props_fine)

  saveRDS(pseudobulk_fine, file = file.path(dir_pseudobulk,
                                            paste0("pseudobulk_", dataset, "_bydonor_finecelltypes.rds")))

}
