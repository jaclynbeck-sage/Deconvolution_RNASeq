# This script maps single cell datasets onto a reference dataset (Cain 2020), in
# order to get consistent broad and fine cell-type assignments across all the
# data sets. We use the broad_class and sub_class definitions from the Cain
# data set.
# Note: This script requires > 64 GB of RAM for some mappings.
library(Seurat)
library(SingleCellExperiment)
library(stringr)
library(dplyr)

source(file.path("functions", "FileIO_HelperFunctions.R"))
source(file.path("functions", "Step04_Mapping_HelperFunctions.R"))

query_datasets <- c("lau", "leng", "mathys", "seaRef")
reference_dataset <- "cain"

ref_seurat <- Load_MapReference(reference_dataset, "broad_class")
ref_seurat_subclass <- Load_MapReference(reference_dataset, "sub_class")

# Fix labels in the Cain pre-processed data ------------------------------------

# Endothelial cells, pericytes, and VLMCs were grouped under "Vascular" in
# the integrated reference. This just transfers that change onto the original
# pre-processed data and saves the object in the post-processed folder.

sce <- Load_PreprocessedData(reference_dataset, remove_excluded = TRUE)

# The cells should already be in the same order, but just in case
cell_order <- match(sce$cell_id, ref_seurat$cell_id)
sce$broad_class <- ref_seurat$broad_class[cell_order]
sce$broad_class <- factor(sce$broad_class)

subclass_meta <- lapply(ref_seurat_subclass, function(S) {
  return(select(S@meta.data, cell_id, broad_class, sub_class))
})
subclass_meta <- do.call(rbind, subclass_meta)

cell_order <- match(sce$cell_id, subclass_meta$cell_id)
sce$sub_class <- as.character(subclass_meta$sub_class[cell_order])
sce$sub_class <- factor(sce$sub_class)

Save_SingleCell(reference_dataset, sce)

rm(sce)
gc()


# Map the other data sets to the reference -------------------------------------

for (query_dataset in query_datasets) {
  message(str_glue("Mapping dataset {query_dataset}..."))

  query_sce <- Load_PreprocessedData(query_dataset, remove_excluded = TRUE)

  # Set up one seurat object per sample
  query_seurat <- CreateSeuratObject(counts(query_sce),
                                     meta.data = data.frame(colData(query_sce)))
  query_seurat <- SplitObject(query_seurat, split.by = "sample")

  rm(query_sce)
  gc()


  ## Map each individual sample to reference -----------------------------------

  query_seurat <- lapply(query_seurat, function(query) {
    samp_name <- unique(query$sample)
    print(str_glue("{query_dataset} sample {samp_name}..."))

    Map_Cells(query, ref_seurat,
              dims_use = 1:30,
              map_cols = list(broad_class = "broad_class"),
              recompute_residuals = TRUE)
  })

  query_mapped <- Merge_Queries(query_seurat)

  # Save mapping info for later inspection. We don't need to save the entire
  # Seurat object, just the parts that were mapped.
  to_save <- list(metadata = query_mapped@meta.data,
                  reductions = query_mapped@reductions)

  saveRDS(to_save,
          file.path(dir_tmp, str_glue("mapped_{query_dataset}_broad_class.rds")))

  rm(query_seurat, query_mapped)
  gc()


  ## Map sub classes -----------------------------------------------------------

  # Uses all samples together due to low cell counts per sample, rather than
  # splitting by sample

  query_sce <- Load_PreprocessedData(query_dataset, remove_excluded = TRUE)

  cell_order <- match(query_sce$cell_id, to_save$metadata$cell_id)
  metadata <- to_save$metadata[cell_order, ]

  query_seurat <- CreateSeuratObject(counts(query_sce),
                                     meta.data = metadata)
  rm(query_sce)
  gc()

  query_seurat <- SplitObject(query_seurat, split.by = "predicted.broad_class")

  query_seurat <- lapply(names(query_seurat), function(X) {
    Map_Cells(query_seurat[[X]], ref_seurat_subclass[[X]],
              dims_use = 1:30,
              map_cols = list(sub_class = "sub_class"),
              recompute_residuals = TRUE)
  })

  query_mapped <- Merge_Queries(query_seurat)

  to_save <- list(metadata = query_mapped@meta.data,
                  reductions = query_mapped@reductions)
  saveRDS(to_save,
          file.path(dir_tmp, str_glue("mapped_{query_dataset}_sub_class.rds")))

  rm(query_seurat, query_mapped)
  gc()


  ## Save the final result -----------------------------------------------------

  query_sce <- Load_PreprocessedData(query_dataset, remove_excluded = TRUE)

  cell_order <- match(query_sce$cell_id, to_save$metadata$cell_id)
  metadata <- to_save$metadata[cell_order, ]

  query_sce$original_broad_class <- query_sce$broad_class
  query_sce$original_sub_class <- query_sce$sub_class

  query_sce$broad_class <- factor(metadata$predicted.broad_class)
  query_sce$sub_class <- factor(metadata$predicted.sub_class)

  Save_SingleCell(query_dataset, query_sce)

  rm(query_sce)
  gc()
}

# Notes on initial sub_class mapping using all subclasses present in cain:
# lau has all subtypes
# lau has low (<20): Ast.5**, End.3**, End.5, End.6*

# leng has all subtypes
# leng has low (<20): Ast.5, End.3**, End.6

# mathys has no: End.3, End.6
# mathys has low (<20): Ast.4, End.2*, End.4, End.5**, Mic.5

# seaRef has no: Ast.4, Ast.5, End.2, End.3, End.6, Mic.3, Mic.5
# seaRef has low (<20): Ast.2**, End.4*, End.5**, Mic.4, Mic.5**, Oligo.4**

# Low or missing in all data sets: End.3, End.6
# Low or missing in 3 of 4 datasets: Ast.5, End.5
# Low or missing in 2 of 4 datasets: Ast.4, End.2, End.4, Mic.5
# Low or missing in seaRef only: Mic.4, Oligo.4

# After these results, I decided that the endothelial subclusters should
# be combined since most endothelial clusters were missing from at least two
# datasets.

# Mapping after combination of endothelial clusters:

# lau has all subtypes
# lau has low (<20): Ast.5**

# leng has all subtypes
# leng has low (<20): Ast.5

# mathys has all subtypes
# mathys has low (<20): Ast.4, Mic.5*

# seaRef is missing: Ast.4, Ast.5, Mic.3, Mic.5
# seaRef has low (<20): Ast.2**, Mic.4*, Oligo.4**

# These results led me to remove glial subtypes as well, resulting in only
# Excitatory, Inhibitory, and Vascular cells having subtypes.
