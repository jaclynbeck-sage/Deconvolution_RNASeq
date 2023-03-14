# Runs Dtangle on test data with unknown cell proportions, using the best-
# performing parameter sets for each reference data set.
#
# NOTE: This script relies on Dtangle markers, so these must be calculated
# before-hand.
#
# Since the parameter lists are small, we don't bother with parallel execution
# here.
library(Matrix)
library(stringr)
library(dplyr)

source(file.path("functions", "FileIO_HelperFunctions.R"))
source(file.path("functions", "DtangleHSPE_HelperFunctions.R"))
source(file.path("functions", "DtangleHSPE_InnerLoop.R"))

##### Edit this variable to run either dtangle or hspe #####

algorithm <- "dtangle" # "dtangle" or "hspe", all lower case

##### Parameter setup #####

library(str_glue("{algorithm}Sparse"), character.only = TRUE)

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

  best_params <- readRDS(file.path(dir_output,
                                   str_glue("best_params_{dataset}_{granularity}.rds")))

  alg_name <- algorithm
  best_params <- subset(best_params, algorithm == alg_name)
  params_list <- do.call(rbind, lapply(best_params$params, as.data.frame)) %>%
    select(-dataset, -granularity)

  ##### Run with the best-performing parameter sets #####

  dtangle_list <- lapply(1:nrow(params_list), function(R) {
    input_type <- params_list$input_type[R]
    input_list <- Get_DtangleHSPEInput(dataset, datatype, granularity, input_type)

    res <- DtangleHSPE_InnerLoop(Y = input_list[["Y"]],
                                 pure_samples = input_list[["pure_samples"]],
                                 params = cbind(params_loop1[P,],
                                                params_list[R,]),
                                 algorithm = algorithm,
                                 limit_n_markers = FALSE)
    return(res)
  })

  names(dtangle_list) <- paste0("dtangle_",
                                str_glue("{dataset}_{granularity}_{datatype}_"),
                                1:length(dtangle_list))

  print(str_glue("Saving final list for {datatype} / {dataset} {granularity}..."))
  Save_AlgorithmOutputList(dtangle_list, "dtangle", dataset, datatype, granularity)

  rm(dtangle_list)
  gc()
}
