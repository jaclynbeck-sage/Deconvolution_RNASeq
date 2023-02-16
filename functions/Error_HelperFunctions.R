library(tidyr)

extract_dtangle_params <- function(params_list) {
  params <- as.data.frame(params_list) %>%
              extract(params_list, c("input", "marker_meth", "gamma_name", "sum_fn_type", "n_markers", "normtype"),
            "input_([a-z]+)_method_([a-z]+)_gamma_([[:alnum:]]+)_summaryfn_([a-z]+)_nmarkers_([0-9\\.|all]+)_normalization_([[:alnum:]]+)")
  rownames(params) <- params_list
  return(params)
}

extract_deconRNASeq_params <- function(params_list) {
  params <- as.data.frame(params_list) %>%
    extract(params_list, c("filterlvl", "n_markers", "usescale", "normtype"),
            "filterlvl_([0-9]+)_nmarkers_([0-9|\\.]+)_usescale_([[:alnum:]]+)_normalization_([[:alnum:]]+)")
  rownames(params) <- params_list
  return(params)
}
