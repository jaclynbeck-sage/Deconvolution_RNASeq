library(dplyr)
library(foreach)
library(doParallel)

##### Parallel execution setup #####

cores <- 12
cl <- makeCluster(cores, type = "PSOCK", outfile = "")
registerDoParallel(cl)

# Libraries that need to be loaded into each parallel environment
required_libraries <- c("MuSiC", "SummarizedExperiment", "stringr", "dplyr")

#### Parameter setup #####

datasets <- c("cain", "lau", "lengEC", "lengSFG", "mathys", "morabito",
              "seaRef") #, "seaAD")

# We can't run Music2 on training data unless we create a training set where
# each sample is from a pool of control OR AD donors but not both. We currently
# don't generate our training sets that way.
params_loop1 <- expand.grid(dataset = datasets,
                            datatype = c("donors"),
                            granularity = c("broad", "fine"),
                            stringsAsFactors = FALSE) %>% arrange(dataset)

params_loop2 <- expand.grid(music2_fn_type = c("default", "t_statistics"), #, "toast"), # toast is broken
                            ct.cov = c(TRUE, FALSE),
                            centered = c(TRUE, FALSE),
                            normalize = c(TRUE, FALSE))

foreach (P = 1:nrow(params_loop1), .packages = required_libraries) %dopar% {
  source(file.path("functions", "FileIO_HelperFunctions.R"))
  source(file.path("functions", "General_HelperFunctions.R"))

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

  # MuSiC2 needs to know which samples are control and which are AD
  metadata <- colData(sce)

  controls <- unique(metadata$donor[metadata$diagnosis == "Control"])
  case <- unique(metadata$donor[metadata$diagnosis == "AD"])

  # Each dataset / datatype / granularity combo gets its own list
  music_list <- list()

  ##### Iterate through combinations of MuSiC2 arguments #####
  # NOTE: This set of parameters (params_loop2) are all executed in the same
  # thread because they use the same single cell and pseudobulk data

  for (R in 1:nrow(params_loop2)) {
    music2_fn_type <- params_loop2$music2_fn_type[R]
    ct.cov <- params_loop2$ct.cov[R]
    centered <- params_loop2$centered[R]
    normalize <- params_loop2$normalize[R]

    name <- str_glue("{dataset}_{granularity}_{datatype}_{R}")
    result <- NULL

    # The TOAST function has some problems that need extra handling
    if (music2_fn_type == "toast") {
      # There is an error in music2_prop_toast where it references these five
      # variables without defining them or taking them as an argument to the
      # function. Defining them here counts as a global variable, and the
      # function can use them. Not ideal, but things work.
      markers <- NULL
      cell_size = NULL
      # ct.cov is needed and defined above
      # centered is needed and defined above
      # normalize is needed and defined above

      # Fix for music2 calling this function with the wrong coef name
      # TODO this doesn't fix it. music2_prop_toast is broken for now.
      csTest_orig <- TOAST::csTest
      csTest <- function(fitted_model, coef, verbose, ...) {
        return(csTest_orig(fitted_model, coef = "condition", verbose, ...))
      }

      R.utils::reassignInPackage("csTest", pkgName="TOAST", csTest);

      result <- music2_prop_toast(bulk.control.mtx = pseudobulk[, controls],
                                  bulk.case.mtx = pseudobulk[, case],
                                  sc.sce = sce,
                                  select.ct = levels(metadata$celltype),
                                  clusters = "celltype",
                                  samples = "donor")
    }
    else {
      # Workaround for a bug in music2_prop_t_statistics
      exprs <- function(X) {X}

      music2_fn <- music2_prop   # music2_fn_type == "default"
      if (music2_fn_type == "t_statistics") {
        music2_fn <- music2_prop_t_statistics
      }

      result <- music2_fn(bulk.control.mtx = pseudobulk[, controls],
                          bulk.case.mtx = pseudobulk[, case],
                          sc.sce = sce,
                          clusters = "celltype",
                          samples = "donor",
                          select.ct = levels(metadata$celltype),
                          ct.cov = ct.cov,
                          centered = centered,
                          normalize = normalize)
    } # End if/else on music2_fn_type

    result$Est.pctRNA <- ConvertPropCellsToPctRNA(result$Est.prop, A)
    result$params <- cbind(params_loop1[P,], params_loop2[R,])
    music_list[[name]] <- result

    gc()
    print(paste(result$params, collapse = "  "))
  } # End params_loop2 loop

  print(str_glue("Saving final list for {dataset} {datatype} {granularity}..."))
  Save_AlgorithmOutputList(music_list, "music2", dataset, datatype, granularity)

  rm(music_list, pseudobulk, sce)
  gc()
}

stopCluster(cl)
