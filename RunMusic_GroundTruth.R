
library(dplyr)
library(foreach)
library(doParallel)

##### Parallel execution setup #####

cores <- 8
cl <- makeCluster(cores, type = "PSOCK", outfile = "")
registerDoParallel(cl)

# Libraries that need to be loaded into each parallel environment
required_libraries <- c("MuSiC", "SummarizedExperiment", "stringr", "dplyr")

#### Parameter setup #####

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito",
              "seaRef") #, "seaAD")

params_loop1 <- expand.grid(dataset = datasets,
                            datatype = c("donors", "training"),
                            granularity = c("broad", "fine"),
                            stringsAsFactors = FALSE) %>% arrange(dataset)

params_loop2 <- expand.grid(ct.cov = c(TRUE, FALSE),
                            centered = c(TRUE, FALSE),
                            normalize = c(TRUE, FALSE))

# HACK to MuSiC to account for the case where some weights are < 0, which leads
# to error-causing NaNs later. seaRef is the only data set this has happened on,
# so I need to figure out what's actually going on here.
#weight.cal.ct.orig <- MuSiC::weight.cal.ct
#weight.cal.ct <- function (...) {
#  weight = weight.cal.ct.orig(...)
#  weight[weight < 0] = 0
#  return(weight)
#}

#R.utils::reassignInPackage("weight.cal.ct", pkgName="MuSiC", weight.cal.ct);

#### Iterate through parameters in parallel ####

foreach (P = 1:nrow(params_loop1), .packages = required_libraries) %dopar% {
  # These need to be sourced inside the loop for parallel processing
  source(file.path("functions", "General_HelperFunctions.R"))
  source(file.path("functions", "FileIO_HelperFunctions.R"))

  dataset <- params_loop1$dataset[P]
  datatype <- params_loop1$datatype[P]
  granularity <- params_loop1$granularity[P]

  sce <- Load_SingleCell(dataset, granularity, output_type = "counts")

  pseudobulk <- Load_Pseudobulk(dataset, datatype, granularity,
                                output_type = "counts")
  pseudobulk <- assay(pseudobulk, "counts")

  A <- Load_AvgLibSize(dataset, granularity)

  # These SHOULD have the same rownames, but just in case.
  keepgene <- intersect(rownames(sce), rownames(pseudobulk))
  pseudobulk <- as.matrix(pseudobulk[keepgene, ])
  sce <- sce[keepgene, ]

  # Each dataset / datatype / granularity combo gets its own list
  music_list <- list()

  ##### Iterate through combinations of MuSiC arguments #####
  # NOTE: This set of parameters (params_loop2) are all executed in the same
  # thread because they use the same single cell and pseudobulk data

  for (R in 1:nrow(params_loop2)) {
    ct_cov <- params_loop2$ct.cov[R]
    centered <- params_loop2$centered[R]
    normalize <- params_loop2$normalize[R]

    name <- str_glue("{dataset}_{granularity}_{datatype}_{R}")

    result <- music_prop(bulk.mtx = pseudobulk, sc.sce = sce,
                         clusters = "celltype",
                         samples = "donor", verbose = TRUE,
                         ct.cov = ct_cov, centered = centered,
                         normalize = normalize)

    # Convert proportion of cells to percent RNA
    result$Est.pctRNA.weighted <- ConvertPropCellsToPctRNA(result$Est.prop.weighted, A)
    result$Est.pctRNA.allgene <- ConvertPropCellsToPctRNA(result$Est.prop.allgene, A)

    result$params <- cbind(params_loop1[P,], params_loop2[R,])

    # Remove "Weight.gene", "r.squared.full", and "Var.prop". "Weight.gene"
    # especially is a very large array and is unneeded, so this reduces
    # output size.
    music_list[[name]] <- result[c("Est.prop.weighted", "Est.prop.allgene",
                                   "Est.pctRNA.weighted", "Est.pctRNA.allgene",
                                   "params")]

    gc()
    print(name)
  } # End params loop

  print("Saving final list...")
  Save_AlgorithmOutputList(music_list, "music", dataset, datatype, granularity)

  rm(music_list, pseudobulk, sce)
  gc()
}

stopCluster(cl)
