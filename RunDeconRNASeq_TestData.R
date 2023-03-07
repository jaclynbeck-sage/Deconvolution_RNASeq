# Runs DeconRNASeq on test data with unknown cell proportions, using the best-
# performing parameter sets for each reference data set.
#
# NOTE: This script relies on Dtangle markers, so these must be calculated
# before-hand.
#
# Since the parameter lists are small, we don't bother with parallel execution
# here.

library(DeconRNASeq)
library(Matrix)
library(stringr)
library(dplyr)
library(foreach)

source(file.path("functions", "FileIO_HelperFunctions.R"))
source(file.path("functions", "General_HelperFunctions.R"))
source(file.path("functions", "DeconRNASeq_InnerLoop.R"))

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito",
              "seaRef") #, "seaAD")

params_loop1 <- expand.grid(dataset = datasets,
                            datatype = c("ROSMAP"),
                            granularity = c("broad"),#, "fine"),
                            stringsAsFactors = FALSE) %>% arrange(datatype)

for (P in 1:nrow(params_loop1)) {
  dataset <- params_loop1$dataset[P]
  granularity <- params_loop1$granularity[P]
  datatype <- params_loop1$datatype[P]

  # Reference data: signature matrix and Ensembl ID -> Symbol mapping
  genes <- Load_GeneConversion(dataset)

  signature <- Load_SignatureMatrix(dataset, granularity)
  signature <- as.data.frame(signature)

  # Test data
  # Gene names in bulk data are Ensembl IDs. They will get converted to gene
  # symbols in this function, so the rownames should match between bulk data and
  # signature matrix.
  bulk <- Load_BulkData(datatype, genes, output_type = "cpm")
  bulk <- as.data.frame(bulk)

  best_params <- readRDS(file.path(dir_output,
                                   str_glue("best_params_{dataset}_{granularity}.rds")))

  best_params <- subset(best_params, algorithm == "deconRNASeq")
  params_list <- do.call(rbind, lapply(best_params$params, as.data.frame)) %>%
                  select(-dataset, -granularity)

  ##### Run with the best-performing parameter sets #####

  decon_list <- foreach (R = 1:nrow(params_list)) %do% {
    res <- DeconRNASeq_InnerLoop(signature, bulk,
                                 cbind(params_loop1[P,], params_list[R,]))
    return(res)
  }

  names(decon_list) <- paste0("deconRNASeq_",
                              str_glue("{dataset}_{granularity}_{datatype}_"),
                              1:nrow(params_list))

  # Save the completed list
  print(str_glue("Saving final list for {datatype} / {dataset} {granularity}..."))
  Save_AlgorithmOutputList(decon_list, "deconRNASeq", dataset, datatype, granularity)

  rm(decon_list)
  gc()
}
