# Load metadata and counts matrices from Synapse and create a SingleCellExperiment
# object that includes metadata, raw counts, and mappings from
# Ensembl ID -> gene symbol.
library(Matrix)
library(SingleCellExperiment)

source("Filenames.R")
source(file.path("functions", "Preprocess_HelperFunctions.R"))
source(file.path("functions", "General_HelperFunctions.R"))

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito",
              "seaRef")#, "seaAD")

##### Download files #####

for (dataset in datasets) {
  print(str_glue("Creating data set for {dataset}..."))
  files <- DownloadData(dataset)

  ##### Read in metadata file #####

  metadata <- ReadMetadata(dataset, files)
  colnames(metadata) <- c("cellid", "donor", "diagnosis", "broadcelltype", "subcluster")
  rownames(metadata) <- metadata$cellid

  # Final modifications to metadata
  metadata$broadcelltype <- factor(metadata$broadcelltype)
  metadata$subcluster <- factor(metadata$subcluster)
  metadata$donor <- factor(metadata$donor)
  metadata$diagnosis <- factor(metadata$diagnosis)


  ##### Read in matrix of counts #####

  counts <- ReadCounts(dataset, files, metadata)

  # Remove genes that are expressed in less than 3 cells
  ok <- rowSums(counts > 0) >= 3
  counts <- counts[ok,]

  # Make sure metadata is in the same order as counts
  metadata <- metadata[colnames(counts), ]

  # Convert gene symbol to Ensembl ID.
  genes <- GeneSymbolToEnsembl(dataset, files, rownames(counts))
  genes <- genes[rownames(counts),]


  ##### Create SingleCellExperiment object and save #####

  if (!all(metadata$cellid == colnames(counts))) {
    stop("*** Cells in the metadata and counts matrix are in different orders.
         Double-check the data. ***")
  }

  sce <- SingleCellExperiment(assays = list(counts = counts),
                              colData = metadata, rowData = genes)

  # If the counts matrix is a DelayedArray (i.e. from the SEA-AD files), the
  # sce file will contain a pointer to the original data file rather than
  # writing the full data to disk again
  saveRDS(sce, file = file.path(dir_input, str_glue("{dataset}_sce.rds")))

  # Calculate the "A" matrix that is needed to convert propCells to pctRNA
  A_broad <- CalculateA(sce, metadata$donor, metadata$broadcelltype)
  A_fine <- CalculateA(sce, metadata$donor, metadata$subcluster)

  saveRDS(list("A_broad" = A_broad, "A_fine" = A_fine),
          file = file.path(dir_input, str_glue("{dataset}_A_matrix.rds")))

  # Calculate a signature for each cell type. This matrix includes all genes in
  # the data set and isn't filtered at this point.
  sig_broad <- CalculateSignature(sce, metadata$donor, metadata$broadcelltype)
  sig_fine <- CalculateSignature(sce, metadata$donor, metadata$subcluster)

  saveRDS(list("sig_broad" = sig_broad, "sig_fine" = sig_fine),
          file = file.path(dir_input, str_glue("{dataset}_signature.rds")))
  print("Done")
}
