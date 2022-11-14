#### Systematic generation of pseudobulks from clusters ###
# TODO turn the two for loops that generate random broad/fine samples into a
# single function. Can make a new metadata DF that contains "cell" and "type"
# that is populated by meta$broadcelltype or meta$subcluster as appropriate,
# and pass it into this function.
# Also pass in: fullmat, ints, numreps, numcells
# Sacrificed readability for speed. Doing matrix multiplication is >10x faster
# than using rowSums().

library(Matrix)
library(SingleCellExperiment)
library(SummarizedExperiment) # TODO or is a DESeq2 object better so we can use normalization functions?
library(stringr)

source("Filenames.R")
source(file.path("functions", "CreatePseudobulk_ByDonor.R"))
source(file.path("functions", "CreatePseudobulk_PureSamplesByDonor.R"))
source(file.path("functions", "CreatePseudobulk_Training.R"))

datasets <- list("mathys") #, "cain", "lau", "leng_SFG", "leng_EC", "lau", "morabito")

for (dataset in datasets) {
  sce <- readRDS(file.path(dir_input, paste0(dataset, "_sce.rds")))
  metadata <- colData(sce)

  CreatePseudobulk_ByDonor(singlecell_counts = counts(sce), metadata = metadata,
                           dataset = dataset, dir_pseudobulk = dir_pseudobulk)

  CreatePseudobulk_PureSamplesByDonor(singlecell_counts = counts(sce),
                                      metadata = metadata, dataset = dataset,
                                      dir_pseudobulk = dir_pseudobulk)

  ### top level - broad cell types ###
  broadtypes <- levels(metadata$broad.cell.type)

  numreps <- 10
  numcells <- 10000 # TODO -- or just use the number of cells available per cell type? Or number of cells in data set?

  ints <- seq(from = 0, to = 1, by = 0.1)

  # Dummy rows/cols so we can rbind/cbind below
  pseudobulk <- Matrix(0, nrow = nrow(counts(sce)), sparse = FALSE)
  rownames(pseudobulk) <- rownames(counts(sce))

  propval <- Matrix(0, ncol = length(broadtypes), sparse = FALSE)
  colnames(propval) <- broadtypes

  # Now add randomly-sampled data sets to fill out pseudobulk data

  for (ii in 1:length(broadtypes)) {

    for (jj in ints) {
      result <- CreatePseudobulk_Training(singlecell_counts = counts(sce),
                                          metadata = metadata,
                                          main_celltype = broadtypes[ii],
                                          proportion = jj, numreps = numreps)
      pseudobulk <- cbind(pseudobulk, result[["counts"]])
      propval <- rbind(propval, result[["propval"]])

      print(c(dataset, broadtypes[ii], jj))
    }
  }

  # Get rid of dummy columns/rows
  pseudobulk <- pseudobulk[,-1]
  propval <- propval[-1,]
  propval <- as.data.frame.matrix(propval)

  se <- SummarizedExperiment(assays = SimpleList(counts = pseudobulk), colData = propval)

  saveRDS(se, file = file.path(dir_pseudobulk, paste0("pseudobulk_", dataset, "_training_broadcelltypes.rds")))
} # Artifically closing for loop to avoid running code below



  ###finer subtypes##
  finetypes <- unique(metadata$fine.cell.type)
  pure_samples_fine <- list()
  for (ii in finetypes) {
    pure_samples_fine[[ii]] <- which(metadata$fine.cell.type == ii)
  }

  numreps <- 10
  numcells <- 1000 # TODO adjust this number since the max fine subtype might be smaller than 20000

  ints <- seq(from = 0, to = 0.3,by = 0.02)

  pseudobulk <- matrix(0, nrow = nrow(counts(sce)),
                       ncol = numreps * length(ints) * length(pure_samples_fine))
  index <- 1
  propval <- matrix(0, nrow = ncol(pseudobulk), ncol = length(pure_samples_fine))
  colnames(propval) <- names(pure_samples_fine)

  for (ii in 1:length(pure_samples_fine)) {
    cellkeep <- pure_samples_fine[[ii]]
    othercells <- setdiff(1:nrow(metadata), cellkeep)

    for (jj in ints) {
      for (kk in 1:numreps) {
        set.seed(ii*10000+jj*100+kk)
        numtokeep <- round(jj*numcells)
        keepinds=sample(cellkeep,numtokeep, replace = T)

        ###remaining cells###
        numtofill <- numcells - numtokeep
        keepinds2 <- sample(othercells, numtofill, replace = T)

        keepinds <- c(keepinds, keepinds2)
        pseudobulk[, index] <- rowSums(counts(sce)[, keepinds])
        colnames(pseudobulk)[index] <- paste(names(pure_samples_broad)[ii], jj, kk, sep = "_")

        tab1 <- table(metadata$fine.cell.type[keepinds])
        propval[index, names(tab1)] = tab1 / sum(tab1)
        index <- index + 1
        print(c(dataset, "fine", index - 1))
      }
    }
  }

  rownames(propval) <- colnames(pseudobulk)

  save(pseudobulk, propval, file = file.path(dir_outpu, paste0("pseudobulk_",dataset,"_finecelltypes_30percentlimit.rda")))

