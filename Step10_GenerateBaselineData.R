# This script generates a data set of 'random' estimates for each bulk data set,
# for use as a baseline in benchmarking. The random estimates are generated 3
# different ways:
#   1. Random cell type percentages drawn from a uniform distribution, with
#      all cell types treated equally.
#   2. An 'educated' random guess drawn from a gaussian distribution, centered
#      around the mean percentage of each cell type in one of the single cell
#      data sets. The standard deviation of each distribution is equal to the
#      standard deviation of that cell type in the reference data.
#   3. One cell type is drawn from a gaussian distribution centered at 0.95 with
#      sd = 0.05, and the other cell types are drawn from gaussian distributions
#      with mean = 0.05/num_celltypes and sd = 0.05/num_celltypes.
# All estimates are normalized per sample so they sum to 1.
#
# This script also creates a 'zero' dataset, where all estimates are 0%. This
# gives us a lower bound (correlation) or upper bound (rMSE, mAPE) on how bad
# the error would be if we didn't guess at all.

library(Matrix)
library(SummarizedExperiment)
library(dplyr)
library(stringr)

source(file.path("functions", "General_HelperFunctions.R"))

reference_datasets <- c("cain", "lau", "leng", "mathys", "seaRef")
bulk_datasets <- c("Mayo", "MSBB", "ROSMAP")

set.seed(12345)
for (granularity in c("broad_class", "sub_class")) {

  # Get which celltypes are present by looking at the metadata from the Cain
  # dataset, which was used to map cell type labels on to all other single cell
  # datasets
  signature <- Load_SignatureMatrix("cain", granularity, "cpm")
  celltypes <- colnames(signature)

  # Loop over all 3 bulk datasets
  for (bulk_dataset in bulk_datasets) {
    se <- Load_BulkData(bulk_dataset, "counts", "none")
    samples <- colnames(se)


    ##### Random uniform dataset #####
    # Generate 100 'parameter sets' of random estimates to mimic the structure
    # of algorithm output

    params <- data.frame("test_data_name" = bulk_dataset,
                         "granularity" = granularity,
                         "method" = "random_uniform")
    name_base <- paste(params, collapse = "_")

    dataset_uniform <- lapply(1:100, function(N) {
      ests <- runif(length(celltypes) * length(samples), min = 0, max = 1)
      ests <- matrix(ests, nrow = length(samples),
                     dimnames = list(samples, celltypes))

      # Sum to 1
      ests <- sweep(ests, 1, rowSums(ests), "/")

      return(list("estimates" = ests, "params" = params))
    })

    names(dataset_uniform) <- paste(name_base, 1:100, sep = "_")

    Save_AlgorithmOutputList(dataset_uniform, algorithm = "Baseline",
                             test_dataset = bulk_dataset, name_base = name_base)


    ##### Educated guesses #####
    # 20 per reference data set = 100 total

    params <- data.frame("test_data_name" = bulk_dataset,
                         "granularity" = granularity,
                         "method" = "random_educated")
    name_base <- paste(params, collapse = "_")

    dataset_educated <- list()
    for (reference_dataset in reference_datasets) {
      # Add reference dataset name to params
      params_tmp <- cbind(params, "reference_data_name" = reference_dataset)

      # Percentages are stored in the pseudobulk metadata
      pb <- Load_Pseudobulk(reference_dataset, "sc_samples", granularity, "counts")
      pcts <- metadata(pb)$pctRNA

      means <- colMeans(pcts)
      sds <- colSds(pcts)

      # 20 param sets with random estimates based on mean and sd of reference
      # dataset.
      list_tmp <- lapply(1:20, function(N) {
        ests <- rnorm(length(samples) * length(celltypes), mean = means, sd = sds)
        ests <- matrix(ests, nrow = length(samples), byrow = TRUE,
                       dimnames = list(samples, celltypes))
        ests[ests < 0] <- 0

        # Sum to 1
        ests <- sweep(ests, 1, rowSums(ests), "/")

        return(list("estimates" = ests, "params" = params_tmp))
      })

      names(list_tmp) <- paste(name_base, reference_dataset, 1:20, sep = "_")
      dataset_educated <- append(dataset_educated, list_tmp)
    }

    Save_AlgorithmOutputList(dataset_educated, algorithm = "Baseline",
                             test_dataset = bulk_dataset, name_base = name_base)


    ##### Bias toward one cell type #####
    # 10 per cell type = 70 total for broad cell types, 240 for fine cell types
    params <- data.frame("test_data_name" = bulk_dataset,
                         "granularity" = granularity,
                         "method" = "random_biased")
    name_base <- paste(params, collapse = "_")

    dataset_bias <- list()
    for (ct in 1:length(celltypes)) {
      # Add the bias cell type name to params
      params_tmp <- cbind(params, "bias_celltype" = celltypes[ct])

      list_tmp <- lapply(1:10, function(N) {
        # Generate a matrix with all low values, then replace the target cell type
        # column with high values
        mean_sd <- 0.05 / length(celltypes)
        ests <- rnorm(length(samples) * length(celltypes), mean = mean_sd, sd = mean_sd)
        ests <- matrix(ests, nrow = length(samples),
                       dimnames = list(samples, celltypes))

        bias_vals <- rnorm(length(samples), mean = 0.95, sd = 0.05)
        ests[,ct] <- bias_vals

        ests[ests < 0] <- 0

        # Sum to 1
        ests <- sweep(ests, 1, rowSums(ests), "/")
        return(list("estimates" = ests, "params" = params_tmp))
      })

      names(list_tmp) <- paste(name_base, celltypes[ct], 1:10, sep = "_")
      dataset_bias <- append(dataset_bias, list_tmp)
    }

    Save_AlgorithmOutputList(dataset_bias, algorithm = "Baseline",
                             test_dataset = bulk_dataset, name_base = name_base)


    ##### Data set of all 0 estimates #####

    params <- data.frame("test_data_name" = bulk_dataset,
                         "granularity" = granularity,
                         "method" = "zeros")
    name_base <- paste(params, collapse = "_")

    zeros <- matrix(rep(0, length(samples) * length(celltypes)),
                    nrow = length(samples),
                    dimnames = list(samples, celltypes))

    dataset_zeros <- list(list("estimates" = zeros, "params" = params))

    names(dataset_zeros) <- name_base

    Save_AlgorithmOutputList(dataset_zeros, algorithm = "Baseline",
                             test_dataset = bulk_dataset, name_base = name_base)
  }
}



