# Runs MuSiC on test data with unknown cell proportions, using the best-
# performing parameter sets for each reference data set.
#
# Since the parameter lists are small, we don't bother with parallel execution
# here.

library(MuSiC)
library(SingleCellExperiment)
library(dplyr)

source(file.path("functions", "FileIO_HelperFunctions.R"))
source(file.path("functions", "General_HelperFunctions.R"))
source(file.path("functions", "Music_InnerLoop.R"))

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

  # Reference data: single cell data and Ensembl ID -> Symbol mapping
  genes <- Load_GeneConversion(dataset)

  sce <- Load_SingleCell(dataset, granularity, output_type = "counts")
  A <- Load_AvgLibSize(dataset, granularity)

  # Test data
  # Gene names in bulk data are Ensembl IDs. They will get converted to gene
  # symbols in this function, so the rownames should match between bulk data and
  # signature matrix.
  bulk <- Load_BulkData(datatype, genes, output_type = "counts")

  keepgene <- intersect(rownames(sce), rownames(bulk))

  bulk <- as.matrix(bulk[keepgene, ]) # Needs to be a matrix, not data.frame
  sce <- sce[keepgene, ]

  best_params <- readRDS(file.path(dir_output,
                                   str_glue("best_params_{dataset}_{granularity}.rds")))

  best_params <- subset(best_params, algorithm == "music")
  params_list <- do.call(rbind, lapply(best_params$params, as.data.frame)) %>%
    select(-dataset, -granularity)

  ##### Run with the best-performing parameter sets #####

  music_list <- lapply(1:nrow(params_list), function(R) {
    res <- Music_InnerLoop(sce, bulk, A, cbind(params_loop1[P,], params_list[R,]))
    return(res)
  })

  # It's possible for some items in music_list to be null if there was an error.
  # Filter them out.
  music_list <- music_list[lengths(music_list) > 0]

  names(music_list) <- paste0("music_",
                              str_glue("{dataset}_{granularity}_{datatype}_"),
                              1:length(music_list))

  print(str_glue("Saving final list for {datatype} / {dataset} {granularity}..."))
  Save_AlgorithmOutputList(music_list, "music", dataset, datatype, granularity)

  rm(music_list, bulk, sce)
  gc()
}
