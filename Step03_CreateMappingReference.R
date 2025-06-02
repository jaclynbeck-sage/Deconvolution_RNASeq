# This script creates references for cell type mapping from the Cain 2020
# single cell dataset. This dataset has a good number of cells in each broad
# class, including microglia and vascular cells, and has well-annotated
# subclasses.
#
# We create one reference for broad cell type mapping and another reference
# for sub class mapping. For broad class mapping, the dataset is normalized
# using SCTransform() on the whole dataset at once.
#
# For sub class mapping, a reference is made for each broad class by itself
# so that the variable features pick up markers for rarer cell types within
# the broad class, rather than picking up differences between broad classes.
# For each broad class, the cells belonging to that broad class are SCTransformed
# as a group.
# Examination of UMAPs of each SCTransformed broad class shows fairly good
# grouping by subclass, so that should be sufficient for mapping.
#
# Note: At this time only Excitatory, Inhibitory, and Vascular cell types are
# broken out into subclasses, as the glial cells (Astrocytes, Microglia, OPCs,
# Oligodendrocytes) do not have well-defined subclasses.
#
# Note: This requires > 64 GB of RAM to do the SCTransform.

library(Seurat)
library(SingleCellExperiment)
library(Matrix)
library(stringr)
library(dplyr)

source(file.path("functions", "FileIO_HelperFunctions.R"))

reference_dataset <- "cain"

sce <- Load_PreprocessedData(reference_dataset, remove_excluded = TRUE)

# Create a "Vascular" broad class for endo, peri, and VLMCs --------------------

metadata <- colData(sce)
metadata$broad_class <- as.character(metadata$broad_class)
metadata$sub_class <- as.character(metadata$sub_class)

# Re-group endothelial cells, pericytes, and VLMCs under "Vascular" at the
# broad_class level. Remove endothelial subtypes.
endos <- metadata$broad_class == "Endothelial"
metadata$sub_class[endos] <- "Endothelial"

vascular <- metadata$broad_class %in% c("Endothelial", "Pericyte", "VLMC")
metadata$broad_class[vascular] <- "Vascular"

# Remove glial subclasses due to poorly-defined borders between subclasses and
# low or missing subclasses in mapped datasets
for (ct in c("Astrocyte", "Microglia", "Oligodendrocyte")) {
  metadata$sub_class[metadata$broad_class == ct] <- ct
}

metadata$broad_class <- factor(metadata$broad_class)
metadata$sub_class <- factor(metadata$sub_class)

colData(sce) <- metadata


# Broad cell types -------------------------------------------------------------

## SCTransform Seurat object ---------------------------------------------------

seurat <- CreateSeuratObject(counts(sce),
                             meta.data = data.frame(colData(sce)))

n_cells <- round(ncol(seurat) * 0.1)
seurat <- SCTransform(seurat,
                      ncells = n_cells,
                      method = "glmGamPoi",
                      do.correct.umi = FALSE,
                      conserve.memory = TRUE)

seurat <- RunPCA(seurat) %>%
  RunUMAP(dims = 1:30, return.model = TRUE)

metadata_fixed <- seurat@meta.data

seurat <- DietSeurat(seurat, layers = c("data", "scale.data"),
                     assays = "SCT",
                     dimreducs = c("pca", "umap"))

Save_MapReference(reference_dataset, seurat, "broad_class")


# Fine cell types --------------------------------------------------------------

sce <- Load_PreprocessedData(reference_dataset, remove_excluded = TRUE)

seurat <- CreateSeuratObject(counts(sce),
                             meta.data = metadata_fixed)
seurat <- SplitObject(seurat, split.by = "broad_class")

# Number of PCs to use was determined by visual inspection of ElbowPlot
pcs <- list("Astrocyte" = 30,
            "Excitatory" = 30,
            "Inhibitory" = 30,
            "Microglia" = 15,
            "Oligodendrocyte" = 20,
            "OPC" = 15,
            "Vascular" = 20)

seurat <- lapply(seurat, function(S) {
  bc <- unique(as.character(S$broad_class))

  S <- SCTransform(S, method = "glmGamPoi", do.correct.umi = FALSE, conserve.memory = TRUE) %>%
    RunPCA() %>%
    RunUMAP(dims = 1:pcs[[bc]], return.model = TRUE)

  S <- DietSeurat(S, counts = FALSE, data = TRUE,
                  scale.data = TRUE, assays = "SCT",
                  dimreducs = c("pca", "umap"))

  return(S)
})

# This is a list of Seurat objects instead of a single object
Save_MapReference(reference_dataset, seurat, "sub_class")

rm(seurat, sce)
gc()
