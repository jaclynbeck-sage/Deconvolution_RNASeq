# Load metadata and counts matrices from Synapse and create either a
# SingleCellExperiment (single cell data) or a Summarized experiment (bulk data)
# object that includes metadata, raw counts, and mappings from
# Ensembl ID -> gene symbol.
library(Matrix)
library(SummarizedExperiment)
library(SingleCellExperiment)

source(file.path("functions", "Step01_Preprocess_HelperFunctions.R"))
source(file.path("functions", "General_HelperFunctions.R"))

##### Setup #####

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito",
              "seaRef", "seaAD", # Single cell
              "Mayo", "MSBB", "ROSMAP") # Bulk

is_singlecell <- function(dataset) {
  return(!is_bulk(dataset))
}

is_bulk <- function(dataset) {
  return(dataset %in% c("Mayo", "MSBB", "ROSMAP"))
}

##### Download files and process data #####

for (dataset in datasets) {
  print(str_glue("Creating data set for {dataset}..."))
  files <- DownloadData(dataset)

  ##### Read in metadata file #####

  metadata <- ReadMetadata(dataset, files)

  if (is_singlecell(dataset)) {
    colnames(metadata) <- c("cellid", "donor", "diagnosis", "broadcelltype", "subcluster")
    rownames(metadata) <- metadata$cellid
  }
  else { # bulk
    colnames(metadata) <- c("donor", "diagnosis", "tissue", "sex")
    rownames(metadata) <- metadata$donor
  }

  # Final modifications to metadata
  for (col in colnames(metadata)) {
    metadata[,col] = factor(metadata[,col])
  }

  ##### Read in matrix of counts #####

  counts <- ReadCounts(dataset, files, metadata)

  # Remove genes that are expressed in less than 3 cells (or samples)
  ok <- rowSums(counts > 0) >= 3
  counts <- counts[ok,]

  ##### Convert gene names and adjust for mitochondrial genes #####

  # Bulk only -- Convert bulk data Ensembl IDs to gene symbols. Some symbols are
  # not unique so we use the first entry in the list for duplicate symbols.
  # All single cell data has symbols as row names already.
  if (is_bulk) {
    genes <- EnsemblIdToGeneSymbol(dataset, files, rownames(counts))
    genes <- genes[rownames(counts),]

    dupes <- duplicated(genes$Symbol)
    genes <- genes[!dupes, ]

    counts <- counts[rownames(genes),]
    rownames(counts) <- genes[rownames(counts), "Symbol"]
  }

  # Remove mitochondrial genes, + remove samples with > 50% mito genes by count.
  # (single cell will all have <50% mito genes, but bulk has higher percentages)
  mt_genes <- grepl("^MT-", rownames(counts))
  pct <- colSums(counts[mt_genes,]) / colSums(counts)
  counts <- counts[!mt_genes, pct <= 0.5]

  # Make sure metadata has the same samples and is in the same order as counts
  metadata <- metadata[colnames(counts), ]


  ##### Bulk data -- create SummarizedExperiment and save #####
  if (is_bulk(dataset)) {
    se <- SummarizedExperiment(assays = list(counts = counts),
                               colData = metadata)

    saveRDS(se, file = file.path(dir_input, str_glue("{dataset}_se.rds")))
  }

  ##### Single cell data -- create SingleCellExperiment object and save #####
  else {
    sce <- SingleCellExperiment(assays = list(counts = counts),
                                colData = metadata)

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
  }

  print("Done")
}
