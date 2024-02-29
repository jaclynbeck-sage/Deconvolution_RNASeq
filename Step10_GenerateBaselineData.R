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

# For direct comparison with the algorithm errors, we need to calculate errors
# on the randomly-generated estimates for each combination of
# reference dataset + normalization + regression_method. The simplest way to
# do this without duplicating error calculation code or other downstream
# processing is to make copies of each of the estimates, one file for each
# combination, even though the files are all the same.
params_permute <- expand.grid(reference_data_name = c("none"),
                              normalization = c("cpm", "tmm", "tpm"),
                              regression_method = c("none", "edger", "deseq2", "dream"),
                              reference_input_type = c("signature"),
                              stringsAsFactors = FALSE)

##### Main loop #####

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
                         "random_method" = "random_uniform")
    name_base <- paste(params, collapse = "_")

    dataset_uniform <- lapply(1:100, function(N) {
      ests <- runif(length(celltypes) * length(samples), min = 0, max = 1)
      ests <- matrix(ests, nrow = length(samples),
                     dimnames = list(samples, celltypes))

      # Sum to 1
      ests <- sweep(ests, 1, rowSums(ests), "/")

      # Needed so error calculations have unique params for each item
      params_tmp <- cbind(params, "trial" = N)

      return(list("estimates" = ests, "params" = params_tmp))
    })

    names(dataset_uniform) <- paste(name_base, 1:100, sep = "_")


    ##### Educated guesses #####
    # 20 per reference data set = 100 total

    params <- data.frame("test_data_name" = bulk_dataset,
                         "granularity" = granularity,
                         "random_method" = "random_educated")
    name_base <- paste(params, collapse = "_")

    dataset_educated <- list()
    for (reference_dataset in reference_datasets) {
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

        # Needed so error calculations have unique params for each item
        params_tmp <- cbind(params, "trial" = paste(reference_dataset, N, sep = "_"))

        return(list("estimates" = ests, "params" = params_tmp))
      })

      names(list_tmp) <- paste(name_base, reference_dataset, 1:20, sep = "_")
      dataset_educated <- append(dataset_educated, list_tmp)
    }


    ##### Bias toward one cell type #####
    # 10 per cell type = 70 total for broad cell types, 240 for fine cell types
    params <- data.frame("test_data_name" = bulk_dataset,
                         "granularity" = granularity,
                         "random_method" = "random_biased")
    name_base <- paste(params, collapse = "_")

    dataset_bias <- list()
    for (ct in 1:length(celltypes)) {
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

        # Needed so error calculations have unique params for each item
        params_tmp <- cbind(params, "trial" = paste(ct, N, sep = "_"))

        return(list("estimates" = ests, "params" = params_tmp))
      })

      names(list_tmp) <- paste(name_base, celltypes[ct], 1:10, sep = "_")
      dataset_bias <- append(dataset_bias, list_tmp)
    }


    ##### Data set of all 0 estimates #####

    params <- data.frame("test_data_name" = bulk_dataset,
                         "granularity" = granularity,
                         "random_method" = "zeros",
                         "trial" = 1)
    name_base <- paste(params, collapse = "_")

    zeros <- matrix(rep(0, length(samples) * length(celltypes)),
                    nrow = length(samples),
                    dimnames = list(samples, celltypes))

    dataset_zeros <- list(list("estimates" = zeros, "params" = params))

    names(dataset_zeros) <- name_base

    full_dataset <- c(dataset_uniform, dataset_educated,
                      dataset_bias, dataset_zeros)


    ##### Write copies of the estimate list #####
    # One copy per combination of test dataset / normalization / regression

    for (R in 1:nrow(params_permute)) {
      # Add one row of params_permute to all the params in the list
      err_list_copy <- full_dataset
      err_list_copy <- lapply(err_list_copy, function(err_item) {
        err_item$params <- cbind(err_item$params, params_permute[R,])
        return(err_item)
      })

      # Create the same naming scheme as other algorithm output
      params_base <- err_list_copy[[1]]$params %>%
        select(reference_data_name, test_data_name, granularity,
               reference_input_type, normalization, regression_method)
      name_base <- paste(params_base, collapse = "_")

      Save_AlgorithmOutputList(err_list_copy, algorithm = "Baseline",
                               test_dataset = bulk_dataset, name_base = name_base)
    }
  }
}
