library(tidyr)

extract_dtangle_params <- function(params_list) {
  params <- as.data.frame(params_list) %>%
              extract(params_list, c("marker_meth", "gamma_name", "sum_fn_type", "n_markers", "normtype"),
            "method_([a-z]+)_gamma_([[:alnum:]]+)_summaryfn_([a-z]+)_nmarkers_([0-9\\.|all]+)_normalization_([[:alnum:]]+)")
  rownames(params) <- params_list
  return(params)
}
