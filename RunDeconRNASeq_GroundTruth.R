# Runs DeconRNASeq on a variety of parameters and saves the list of results to a
# file.
#
# NOTE: This script relies on Dtangle markers, so these must be calculated
# before-hand.
#
# This script uses parallel processing to run each parameter set on its own
# core. To run in serial instead, comment out "registerDoParallel(cl)" below and
# change the foreach loop's "%dopar%" to "%do%".

library(dplyr)
library(foreach)
library(doParallel)

##### Parallel execution setup #####

cores <- 8
cl <- makeCluster(cores, type = "PSOCK", outfile = "")
registerDoParallel(cl)

# Libraries that need to be loaded into each parallel environment
required_libraries <- c("DeconRNASeq", "Matrix", "SummarizedExperiment",
                        "stringr", "dplyr")

#### Parameter setup #####

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito",
              "seaRef") #, "seaAD")

params_loop1 <- expand.grid(dataset = datasets,
                            datatype = c("donors", "training"),
                            granularity = c("broad", "fine"),
                            stringsAsFactors = FALSE) %>% arrange(datatype)

params_loop2 <- expand.grid(filter_level = c(0, 1, 2, 3),
                            n_markers = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.75, 1.0,
                                          10, 50, 100, 200, 500),
                            marker_type = c("ratio", "diff", "p.value", "regression"),
                            use_scale = c(TRUE, FALSE),
                            stringsAsFactors = FALSE) %>% arrange(filter_level)

# Some filter_type / n_markers combos are not valid, get rid of them
# (filter levels 1 & 2 don't use n_markers or marker_type arguments)
params_loop2 <- subset(params_loop2, !(filter_level < 3 &
                                         (n_markers != 1 | marker_type != "ratio")))

params_loop2$marker_type[params_loop2$filter_level < 3] <- "None"
params_loop2$n_markers[params_loop2$filter_level < 3] <- -1

#### Iterate through parameters in parallel ####

foreach (P = 1:nrow(params_loop1), .packages = required_libraries) %dopar% {
  # These need to be sourced inside the loop for parallel processing
  source(file.path("functions", "General_HelperFunctions.R"))
  source(file.path("functions", "FileIO_HelperFunctions.R"))

  dataset <- params_loop1$dataset[P]
  datatype <- params_loop1$datatype[P]
  granularity <- params_loop1$granularity[P]

  signature <- Load_SignatureMatrix(dataset, granularity)
  signature <- as.data.frame(signature)

  pseudobulk <- Load_Pseudobulk(dataset, datatype, granularity, output_type = "cpm")
  pseudobulk <- assay(pseudobulk, "counts")

  # Each dataset / datatype / granularity combo gets its own list
  decon_list <- list()

  ##### Filter levels, number of markers, and DeconRNASeq arguments #####
  # NOTE: This set of parameters (params_loop2) are all executed in the same
  # thread because they use the same signature and pseudobulk data

  for (R in 1:nrow(params_loop2)) {
    filter_level <- params_loop2$filter_level[R]
    n_markers <- params_loop2$n_markers[R]
    marker_type <- params_loop2$marker_type[R]
    use_scale <- params_loop2$use_scale[R]

    signature_filt <- FilterSignature(signature, filter_level, dataset,
                                      granularity, n_markers, marker_type)

    keepgene <- intersect(rownames(signature_filt), rownames(pseudobulk))
    pseudobulk_filt <- as.data.frame(as.matrix(pseudobulk[keepgene,]))

    name <- str_glue("{dataset}_{granularity}_{datatype}_{R}")

    res <- DeconRNASeq(pseudobulk_filt, signature_filt, proportions = NULL,
                       known.prop = FALSE, use.scale = use_scale, fig = FALSE)

    rownames(res$out.all) <- colnames(pseudobulk_filt)
    res$params <- cbind(params_loop1[P,], params_loop2[R,])

    decon_list[[name]] <- res
    gc()

    print(paste(res$params, collapse = "  "))
  } # end params loop 2

  # Save the completed list
  print(str_glue("Saving final list for {dataset} {datatype} {granularity}..."))
  Save_AlgorithmOutputList(decon_list, "deconRNASeq", dataset, datatype, granularity)
  rm(decon_list)
  gc()

  return(NULL)
}

stopCluster(cl)
