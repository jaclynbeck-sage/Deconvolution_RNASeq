# This script creates references for cell type mapping from the Cain 2020
# single cell dataset. This dataset has a good number of cells in each broad
# class, including microglia and vascular cells, and has well-annotated
# subclasses.
#
# We create one reference for broad cell type mapping and another reference
# for sub class mapping. For broad class mapping, each of the 24 samples is
# SCTransformed individually and then integrated together into one data set.
# This dataset is too large to reasonably find anchors and integrate across
# all combinations of samples. Instead, we use 8 samples as a 'reference',
# which are first integrated together, and then all remaining samples are
# integrated into that reference.
#
# For sub class mapping, a reference is made for each broad class by itself
# so that the variable features pick up markers for rarer cell types within
# the broad class, rather than picking up differences between broad classes.
# Each broad class (except for excitatory neurons) has, on average, < 1500 cells
# per sample, so we do not integrate these for sub class mapping. Integrating
# them causes massive overcorrection in the non-neuronal cell types, and the
# subclasses are no longer distinct. Instead, broad classes are SCTransformed
# only. Examination of UMAPs of each SCTransformed broad class shows fairly good
# grouping by subclass, so that should be sufficient for mapping.
#
# Note: Oligodendrocytes are mostly clustered by subclass but the subclasses do
# not have clear boundaries between them, as noted in the paper, so it is
# unclear whether these subclasses will provide good mappings.
#
# Note: This requires nearly 128 GB of RAM to integrate the dataset.

library(Seurat)
library(SingleCellExperiment)
library(Matrix)
library(stringr)
library(dplyr)

source(file.path("functions", "FileIO_HelperFunctions.R"))

reference_dataset <- "cain"
do_correction <- FALSE

sce <- Load_PreprocessedData(reference_dataset, remove_excluded = TRUE)

##### Create a "Vascular" broad class for endo, peri, and VLMCs #####

metadata <- colData(sce)
metadata$broad_class <- as.character(metadata$broad_class)
metadata$sub_class <- as.character(metadata$sub_class)

# Re-group endothelial cells, pericytes, and VLMCs under "Vascular" at the
# broad_class level
endos <- metadata$broad_class == "Endothelial"
metadata$sub_class[endos] <- "Endothelial"

vascular <- metadata$broad_class %in% c("Endothelial", "Pericyte", "VLMC")
metadata$broad_class[vascular] <- "Vascular"

metadata$broad_class <- factor(metadata$broad_class)
metadata$sub_class <- factor(metadata$sub_class)

colData(sce) <- metadata


########## Broad cell types ##########

##### SCTransform Seurat object #####

seurat <- CreateSeuratObject(counts(sce),
                             meta.data = data.frame(colData(sce)))

n_cells <- round(ncol(seurat) * 0.1)
seurat <- SCTransform(seurat, ncells = n_cells, method = "glmGamPoi",
                      do.correct.umi = FALSE)

seurat <- RunPCA(seurat) %>%
            RunUMAP(dims = 1:30, return.model = TRUE) %>%
            FindNeighbors(dims = 1:30) %>%
            FindClusters(resolution = 0.1)


##### Correct cell type labels #####

if (do_correction) {
  # Assign cluster labels based on the majority cell type in each cluster
  confusion <- table(seurat$seurat_clusters, seurat$broad_class)
  clust_labels <- sapply(rownames(confusion), function(R) {
    colnames(confusion)[which.max(confusion[R,])]
  })

  seurat$new_labels <- clust_labels[seurat$seurat_clusters]

  mismatches <- sum(seurat$new_labels != seurat$broad_class)
  mm_pct <- round(mismatches / ncol(seurat) * 100, digits = 2)
  print(paste0(str_glue("{mismatches} cells ({mm_pct}%) are mis-labeled at the ",
                        "broad class level and will be corrected in the map ",
                        "reference.")))

  # Update cell type assignments to match the clustering
  seurat$original_labels <- seurat$broad_class
  seurat$broad_class <- factor(seurat$new_labels)
}

metadata_fixed <- seurat@meta.data

seurat <- DietSeurat(seurat, counts = FALSE, data = TRUE,
                     scale.data = TRUE, assays = "SCT",
                     dimreducs = c("pca", "umap"))

Save_MapReference(reference_dataset, seurat, "broad_class")


########## Fine cell types ##########

sce <- Load_PreprocessedData(reference_dataset, remove_excluded = TRUE)

seurat <- CreateSeuratObject(counts(sce),
                             meta.data = metadata_fixed)
seurat <- SplitObject(seurat, split.by = "broad_class")

pcs <- list("Astrocyte" = 30,
            "Excitatory" = 30,
            "Inhibitory" = 30,
            "Microglia" = 15,
            "Oligodendrocyte" = 20,
            "OPC" = 15,
            "Vascular" = 20)

res <- list("Astrocyte" = 0.6,
            "Excitatory" = 0.2,
            "Inhibitory" = 0.1,
            "Microglia" = 0.4,
            "Oligodendrocyte" = 0.6,
            "OPC" = 0.1,
            "Vascular" = 0.1)

seurat <- lapply(seurat, function(S) {
  bc <- unique(as.character(S$broad_class))

  S <- SCTransform(S, method = "glmGamPoi", do.correct.umi = FALSE) %>%
          RunPCA() %>%
          RunUMAP(dims = 1:pcs[[bc]], return.model = TRUE) %>%
          FindNeighbors(dims = 1:pcs[[bc]]) %>%
          FindClusters(resolution = res[[bc]])

  if (do_correction) {
    # Assign cluster labels based on majority cell type in each cluster
    confusion <- table(S$seurat_clusters, as.character(S$sub_class))
    clust_labels <- sapply(rownames(confusion), function(R) {
      colnames(confusion)[which.max(confusion[R,])]
    })

    S$new_labels <- clust_labels[S$seurat_clusters]
    mismatches <- sum(S$new_labels != S$sub_class)
    mm_pct <- round(mismatches / ncol(S) * 100, digits = 2)
    print(paste0(str_glue("{bc}: {mismatches} cells ({mm_pct}%) are mis-labeled ",
                          "at the sub class level and will be corrected in the ",
                          "map reference.")))
    S$original_labels <- S$sub_class
    S$sub_class <- S$new_labels
  }

  S <- DietSeurat(S, counts = FALSE, data = TRUE,
                  scale.data = TRUE, assays = "SCT",
                  dimreducs = c("pca", "umap"))

  return(S)
})

# This is a list of Seurat objects instead of a single object
Save_MapReference(reference_dataset, seurat, "sub_class")

rm(seurat, sce)
gc()
