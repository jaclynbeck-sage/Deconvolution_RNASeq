# This requires nearly 128 GB of RAM and sometimes tries to use more
library(Seurat)
library(SingleCellExperiment)
library(Matrix)
library(Azimuth)
library(stringr)

source(file.path("functions", "FileIO_HelperFunctions.R"))

sce <- Load_SingleCell("seaRef", granularity = "broad_class", output_type = "counts")
donors <- levels(sce$donor)

# For memory conservation purposes, use the donor with the most cells as
# the reference for anchor selection / integration as suggested by
# https://satijalab.org/seurat/articles/integration_large_datasets
n_cells <- table(sce$donor)
ref_donor <- which.max(n_cells)

metadata <- as.data.frame(colData(sce))

#### SCTransform each donor's data

for (donor in donors) {
  donor_metadata <- metadata[metadata$donor == donor,]

  seurat <- CreateSeuratObject(counts(sce)[,as.character(donor_metadata$cell_id)],
                               meta.data = donor_metadata,
                               min.cells = 3)

  seurat <- SCTransform(seurat, method="glmGamPoi", variable.features.n = 5000)

  # Save in case of crashing
  saveRDS(seurat, file.path(dir_tmp, str_glue("sct_{donor}.rds")))

  rm(seurat, donor_metadata)
  gc()
}

rm(sce)
gc()

#### Integrate all donors

seurat_list <- lapply(donors, function(donor) {
  readRDS(file.path(dir_tmp, str_glue("sct_{donor}.rds")))
})

features <- SelectIntegrationFeatures(object.list = seurat_list,
                                      nfeatures = 5000)
seurat_list <- PrepSCTIntegration(object.list = seurat_list,
                                  anchor.features = features)

anchors <- FindIntegrationAnchors(
  object.list = seurat_list,
  anchor.features = features,
  normalization.method = "SCT",
  dims = 1:30,
  reference = ref_donor
)

gc()

# Reduce size by getting rid of RNA assay and unused SCT counts matrix
anchors@object.list <- lapply(anchors@object.list, function(obj) {
  obj <- DietSeurat(obj, counts = FALSE, data = TRUE, scale.data = TRUE,
                    assays = c("SCT"))
  return(obj)
})
saveRDS(anchors, file.path(dir_tmp, "integration_anchors.rds"))

rm(seurat_list, features)
gc()

integrated <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT",
  dims = 1:30
)

gc()

integrated <- RunPCA(integrated, verbose = FALSE)
integrated <- RunUMAP(integrated, dims = 1:20)
saveRDS(integrated, file.path(dir_tmp, "seaRef_integrated_5000genes.rds"))

Idents(integrated) <- integrated$sub_class
integrated$broad_class <- factor(integrated$broad_class)
integrated$fine_cluster <- factor(integrated$fine_cluster)
integrated$sub_class <- factor(integrated$sub_class)

# Downsample to 2000 cells per cell type (at the sub_class level)
set.seed(12345)
ref <- subset(integrated, downsample = 2000)
gc()

ref <- RunUMAP(
  object = ref,
  reduction = "pca",
  dims = 1:30,
  return.model = TRUE
)
gc()

ref <- AzimuthReference(
  object = ref,
  refUMAP = "umap",
  refDR = "pca",
  refAssay = "integrated",
  metadata = c("broad_class", "sub_class", "fine_cluster"),
  dims = 1:50,
  k.param = 31,
  reference.version = "1.0.0"
)

# The files have to be named "idx.annoy" and "ref.Rds". The "R" has to be
# capitalized or Azimuth complains
SaveAnnoyIndex(ref[["refdr.annoy.neighbors"]],
               file = file.path(dir_azimuth_reference, "idx.annoy"))
saveRDS(ref, file = file.path(dir_azimuth_reference, "ref.Rds"))
