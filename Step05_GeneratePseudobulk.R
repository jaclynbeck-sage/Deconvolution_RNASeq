# This script generates pseudobulk data sets from single cell data sets. The
# following pseudobulk data sets are created for both broad and fine cell types:
#   1. Pseudobulk of each sample
#   2. Pseudobulk of each cell type + sample combination, to create "pure"
#      pseudobulk samples of each cell type
#   3. A training pseudobulk set with randomly-sampled cells
#
# The training pseudobulk set is generated as follows:
#   1. For each cell type, determine a reasonable range of proportions to
#      test based on that cell type's abundance in the single cell data set
#   2. For each cell type, for each proportion in {0 to max_proportion}, create
#      5 pseudobulk samples where that cell type has the specified proportion:
#      2a. Cells are randomly sampled from the pool of cells of that cell type
#          to make up the specified proportion of the sample.
#      2b. The other cell types fill in the rest of the sample with
#          randomly-generated proportions, sampling from their respective
#          pools of cell types.
#      2c. The counts from all sampled cells are added together.
#   3. For each pseudobulk sample, the ground truth proportion of cells and
#      percent RNA from each cell type is calculated and added as metadata
#
# This ensures that:
#   1. Each cell type is missing from at least 5 samples in the set
#   2. Each cell type is forced to cover a range of proportions across samples
#   3. Randomness in assigning proportions to each non-main cell type will
#      create an even wider range of proportions than you would get from
#      treating the other cell types as a single pool

library(Matrix)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(stringr)
library(dplyr)
library(edgeR)

source(file.path("functions", "FileIO_HelperFunctions.R"))
source(file.path("functions", "General_HelperFunctions.R"))
source(file.path("functions", "Step05a_CreatePseudobulk_BySample.R"))
source(file.path("functions", "Step05b_CreatePseudobulk_PureSamples.R"))
#source(file.path("functions", "Step05c_CreatePseudobulk_Training.R"))

datasets <- c("cain", "lau", "leng", "mathys", "seaRef")

for (dataset in datasets) {
  sce <- Load_SingleCell(dataset, "broad_class", output_type = "counts")
  metadata <- colData(sce)

  message(str_glue("{dataset}: creating pseudobulk by sample..."))
  CreatePseudobulk_BySample(counts(sce), metadata, dataset)

  message(str_glue("{dataset}: creating pseudobulk pure samples..."))
  CreatePseudobulk_PureSamples(counts(sce), metadata, dataset)

  # Skip making training pseudobulk but keep the code just in case
  next

  for (granularity in c("broad_class", "sub_class")) {
    message(str_glue("{dataset}: creating pseudobulk training set for {granularity} cell types..."))

    # Create a generic "celltype" column that is populated with either broad or
    # fine cell types depending on the for-loop
    if (granularity == "broad_class") {
      metadata$celltype <- metadata$broad_class
    } else if (granularity == "sub_class") {
      metadata$celltype <- metadata$sub_class
    }

    celltypes <- levels(metadata$celltype)

    # 5 samples per proportion
    num_samples <- 5
    n_cells <- table(metadata$celltype)

    # The range of percents we use to create the training set is dependent on
    # how abundant each cell type is in this data set
    pcts <- table(metadata$sample, metadata$celltype)
    pcts <- sweep(pcts, 1, rowSums(pcts), "/")

    # The largest proportion of each cell type from any sample, x 2 determines
    # the range of percents we test for that cell type
    maxs <- colMaxs(pcts, useNames = TRUE)
    maxs <- ceiling(20 * maxs) / 10 # Rounds 2*X to the next-highest 10%
    maxs[maxs > 1] <- 1.0

    # Each cell type gets its own range of percents, divided into 20 increments
    ints <- sapply(maxs, function(X) {
      seq(from = 0, to = X, by = (X / 20))
    })

    pseudobulk <- list()
    propCells <- list()
    pctRNA <- list()

    # Create randomly-sampled data sets to fill out pseudobulk data

    for (ct in 1:length(celltypes)) {
      for (prop in ints[, ct]) {
        # Limit the total number of cells in the resample to be proportional
        # to the number of cells for this cell type. We don't want to resample
        # a population of 100 cells 10,000 times, for example
        numcells <- min(10000, n_cells[celltypes[ct]] / prop)

        result <- CreatePseudobulk_Training(singlecell_counts = counts(sce),
                                            cell_assigns = metadata$celltype,
                                            main_celltype = celltypes[ct],
                                            proportion = prop,
                                            num_cells = numcells,
                                            num_samples = num_samples)
        name <- paste(ct, prop)
        pseudobulk[[name]] <- result$counts
        propCells[[name]] <- result$propCells
        pctRNA[[name]] <- result$pctRNA

        print(c(dataset, celltypes[ct], prop))
      }
    }

    pseudobulk <- do.call(cbind, pseudobulk)
    propCells <- do.call(rbind, propCells)
    pctRNA <- do.call(rbind, pctRNA)

    pseudobulk <- as(pseudobulk, "matrix")
    tmm <- edgeR::normLibSizes(pseudobulk, method = "TMMwsp")

    se <- SummarizedExperiment(assays = SimpleList(counts = pseudobulk),
                               colData = data.frame(tmm_factors = tmm),
                               metadata = list("propCells" = propCells,
                                               "pctRNA" = pctRNA))

    Save_Pseudobulk(se, dataset, "training", granularity)
  }
}
