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

sce <- Load_PreprocessedData(reference_dataset)

##### Create a "Vascular" broad class for endo, peri, and VLMCs #####

metadata <- colData(sce)
metadata$broad_class <- as.character(metadata$broad_class)

# Re-group endothelial cells, pericytes, and VLMCs under "Vascular" at the
# broad_class level
vascular <- metadata$broad_class %in% c("Endothelial", "Pericyte", "VLMC")
metadata$broad_class[vascular] <- "Vascular"

metadata$broad_class <- factor(metadata$broad_class)

colData(sce) <- metadata


########## Broad cell types ##########

##### Choose reference samples for anchors #####

# 1 sample from each diagnosis/sex combination, selecting by highest number of
# cells in each group. This gives 8 reference samples.

covariates <- Load_Covariates(reference_dataset)

meta_covar <- as.data.frame(merge(metadata, covariates, by = "sample"))
cell_counts <- meta_covar %>% group_by(msex, diagnosis, sample) %>%
                summarise(counts = n()) %>%
                summarise(max_sample = sample[which.max(counts)])

ref_samples <- cell_counts$max_sample


##### Prep Seurat object #####

seurat <- CreateSeuratObject(counts(sce),
                             meta.data = data.frame(colData(sce)))
seurat <- SplitObject(seurat, split.by = "sample")

# Make sure the reference samples have the right indices into the Seurat list
ref_samples_index <- match(ref_samples, names(seurat))
print(ref_samples_index)


##### SCTransform each sample's data and find anchors #####

# Use the default 3000 variable features. Corrected UMI are necessary for
# anchor finding so they cannot be removed yet.
seurat <- lapply(seurat, SCTransform, method = "glmGamPoi")

features <- SelectIntegrationFeatures(object.list = seurat,
                                      nfeatures = 3000)
seurat <- PrepSCTIntegration(object.list = seurat,
                                 anchor.features = features)

# Find anchors between all samples and the 8 reference samples
anchors <- FindIntegrationAnchors(object.list = seurat,
                                  anchor.features = features,
                                  normalization.method = "SCT",
                                  dims = 1:30,
                                  reference = ref_samples_index)

# Reduce size in memory by getting rid of RNA assay and unused SCT counts matrix
anchors@object.list <- lapply(anchors@object.list, function(obj) {
  obj <- DietSeurat(obj, counts = FALSE, data = TRUE, scale.data = TRUE,
                    assays = c("SCT"))
  return(obj)
})

# Save in case of crashing
saveRDS(anchors, file.path(dir_tmp, str_glue("{reference_dataset}_anchors.rds")))

rm(seurat, sce, features)
gc()


##### Integrate and get the UMAP #####

ref_seurat <- IntegrateData(anchorset = anchors,
                            normalization.method = "SCT",
                            dims = 1:30)

ref_seurat <- DietSeurat(ref_seurat, counts = FALSE, data = TRUE,
                         scale.data = TRUE, assays = "integrated",
                         dimreducs = c("pca", "umap"))
gc()

ref_seurat <- RunPCA(ref_seurat)

# 25 dims determined by inspection of ElbowPlot. "return.model" is necessary
# for mapping.
ref_seurat <- RunUMAP(ref_seurat, dims = 1:25,
                      return.model = TRUE)

Save_MapReference(reference_dataset, ref_seurat, "broad_class")

rm(ref_seurat, anchors)
gc()


########## Fine cell types ##########

sce <- Load_PreprocessedData(reference_dataset)
colData(sce) <- metadata

seurat <- CreateSeuratObject(counts(sce),
                             meta.data = data.frame(colData(sce)))
seurat <- SplitObject(seurat, split.by = "broad_class")


##### SCTransform and get the UMAP #####

# SCTransform only, no integration for sub classes, to avoid over-correction.
seurat <- lapply(seurat, SCTransform,
                 method = "glmGamPoi",
                 do.correct.umi = FALSE)

seurat <- lapply(seurat, RunPCA)

# Dimensions determined by inspection of ElbowPlot for each broad class
pcs <- list("Astrocyte" = 1:30,
            "Excitatory" = 1:30,
            "Inhibitory" = 1:30,
            "Microglia" = 1:15,
            "Oligodendrocyte" = 1:20,
            "OPC" = 1:15,
            "Vascular" = 1:20)

seurat_names <- names(seurat)
seurat <- lapply(seurat_names, function(N) {
  RunUMAP(seurat[[N]], dims = pcs[[N]], return.model = TRUE)
})

names(seurat) <- seurat_names

# This is a list of Seurat objects instead of a single object
Save_MapReference(reference_dataset, seurat, "sub_class")

rm(seurat, sce)
gc()
