CreatePseudobulk_ByDonor <- function(singlecell_counts, metadata, dataset, dir_pseudobulk) {
  # This is SIGNIFICANTLY faster than calling rowSums on individual donor sets

  y <- model.matrix(~0 + donor, data = metadata)
  counts <- singlecell_counts %*% y
  counts <- as(counts, "dgCMatrix") # For cases where counts is a DelayedMatrix

  # colnames end up as "donordonor<#>" because of model.matrix. Remove the extra "donor".
  colnames(counts) <- str_replace(colnames(counts), "donor", "")

  props <- table(metadata$donor, metadata$broadcelltype)

  # necessary to get the correct shape for Summarized Experiment to convert to DataFrame
  props <- as.data.frame.matrix(props / rowSums(props))

  pseudobulk <- SummarizedExperiment(assays = SimpleList(counts = counts),
                                     colData = props)

  saveRDS(pseudobulk, file = file.path(dir_pseudobulk,
                                       paste0("pseudobulk_", dataset,
                                              "_donors_broadcelltypes.rds")))

  # TODO Is there a smooth way to put both broad and fine cell types in one object
  # without requiring a lot of reshaping/processing when it's read in from a file?
  props_fine <- table(metadata$donor, metadata$subcluster)
  props_fine <- as.data.frame.matrix(props_fine / rowSums(props_fine))

  pseudobulk_fine <- SummarizedExperiment(assays = SimpleList(counts = counts),
                                          colData = props_fine)

  saveRDS(pseudobulk_fine, file = file.path(dir_pseudobulk,
                                            paste0("pseudobulk_", dataset,
                                                   "_donors_finecelltypes.rds")))

}
