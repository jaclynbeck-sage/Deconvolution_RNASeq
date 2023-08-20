# Load metadata and counts matrices from Synapse and create either a
# SingleCellExperiment (single cell data) or a Summarized experiment (bulk data)
# object that includes metadata, raw counts, and mappings from
# Ensembl ID -> gene symbol. Gene names and cell type labels are also
# standardized across data sets.
library(Matrix)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(edgeR)

source(file.path("functions", "Step02_Preprocess_HelperFunctions.R"))
source(file.path("functions", "General_HelperFunctions.R"))

##### Setup #####

datasets <- c("cain", "lau", "leng", "mathys", #"morabito",
              "seaRef", #"seaAD", # Single cell
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

  metadata_list <- ReadMetadata(dataset, files)
  metadata <- metadata_list$metadata

  # For data privacy reasons, each data set's covariates are saved in a separate
  # file from the main SingleCellExperiment object and are not uploaded anywhere
  covariates <- metadata_list$covariates
  saveRDS(covariates, file.path(dir_covariates,
                                str_glue("{dataset}_covariates.rds")))

  if (is_singlecell(dataset)) {
    colnames(metadata)[1:5] <- c("cell_id", "donor", "diagnosis", "broad_class",
                                   "sub_class")

    # seaRef and seaAD are special cases and have a 'fine_cluster' field
    if (ncol(metadata) > 5) {
      colnames(metadata)[6] <- "fine_cluster"
    }

    rownames(metadata) <- metadata$cell_id
  }
  else { # bulk
    colnames(metadata) <- c("donor", "specimenID", "diagnosis", "tissue", "sex")
    rownames(metadata) <- metadata$specimenID
  }


  ##### Read in matrix of counts #####

  counts <- ReadCounts(dataset, files, metadata)

  # Remove genes that are expressed in less than 3 cells (or samples)
  ok <- rowSums(counts > 0) >= 3
  counts <- counts[ok,]


  ##### Convert gene names #####

  # Bulk only -- Convert bulk data Ensembl IDs to gene symbols.
  if (is_bulk(dataset)) {
    genes <- EnsemblIdToGeneSymbol(rownames(counts))
  }
  # Single cell only -- update potentially outdated gene symbols to the most
  # current version possible, and get the matching Ensembl IDs too.
  else {
    genes <- UpdateGeneSymbols(dataset, rownames(counts))
  }

  # Some symbols are not unique so we use the first entry in the list for
  # duplicate symbols, which is what Seurat does and is probably consistent with
  # most of the single cell data sets.
  dupes <- duplicated(genes$hgnc_symbol)
  genes <- genes[!dupes, ]

  # Puts 'genes' in the same gene order as counts so the counts matrix doesn't
  # get rearranged unnecessarily
  genes <- genes[intersect(rownames(counts), rownames(genes)),]

  # Assign rownames to be the canonical hgnc symbol, applies to both bulk and
  # single cell
  counts <- counts[rownames(genes),]
  rownames(counts) <- genes$hgnc_symbol
  rownames(genes) <- genes$hgnc_symbol

  ##### Adjust for and remove mitochondrial genes #####

  # Remove samples with > 50% mito genes by count. (single cell will mostly have
  # <10% mito genes, but bulk has higher percentages)
  mt_genes <- grepl("^MT-", rownames(counts))
  pct <- colSums(counts[mt_genes,]) / colSums(counts)
  counts <- counts[!mt_genes, pct <= 0.5] # Remove mito genes
  genes <- genes[rownames(counts),]

  ##### Final modifications to metadata #####

  # Make sure metadata has the same samples and is in the same order as counts
  metadata <- metadata[colnames(counts), ]

  if (is_singlecell(dataset)) {
    metadata$broadcelltype <- RemapCelltypeNames(metadata$broadcelltype)
  }

  for (col in colnames(metadata)) {
    metadata[,col] = factor(metadata[,col])
  }

  # TMM normalization factors -- unfortunately will convert to dense matrix
  if (is_singlecell(dataset)) {
    tmm <- calcNormFactors(counts, method = "TMMwsp")
  }
  else {
    tmm <- calcNormFactors(counts, method = "TMM")
  }
  gc()
  metadata$tmm_factors <- tmm

  ##### Bulk data -- create SummarizedExperiment and save #####

  # Bulk data also has two more assays in addition to raw counts
  if (is_bulk(dataset)) {
    normalized <- ReadNormCounts_BulkData(files, "normalized", genes, colnames(counts))
    residuals <- ReadNormCounts_BulkData(files, "residuals", genes, colnames(counts))

    se <- SummarizedExperiment(assays = list(counts = counts,
                                             normalized = normalized,
                                             residuals = residuals),
                               colData = metadata,
                               rowData = genes)

    saveRDS(se, file = file.path(dir_input, str_glue("{dataset}_se.rds")))
  }

  ##### Single cell data -- create SingleCellExperiment object and save #####

  else {
    sce <- SingleCellExperiment(assays = list(counts = counts),
                                colData = metadata,
                                rowData = genes)

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
