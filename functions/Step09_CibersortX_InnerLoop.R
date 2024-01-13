# Uses the "omnideconv" library to run CibersortX. CibersortX runs in a docker
# container that requires a username and token, which can be requested on the
# CibersortX website here: https://cibersortx.stanford.edu/download.php. This
# code expects environment variables CIBERSORT_EMAIL and CIBERSORT_TOKEN to be
# set prior to running.
#
# Arguments:
#   reference_filename = the filename (plus full path) of the single cell data
#                        file that was created by Save_SingleCellToCibersort().
#   bulk_mat = a matrix (genes x samples) of bulk data. Must be a dense matrix.
#   cibersort_signature = a matrix (genes x celltypes) of the signature generated
#                         by CibersortX in Step 07 of this pipeline.
#   params = a single-row data frame or a named vector/list of parameters
#            containing the following variables: reference_data_name,
#            test_data_name, granularity, filter_level, n_markers, marker_type,
#            marker_subtype, marker_input_type. This variable is unused except
#            to get added to the results object.
#
# Returns:
#   a list containing entries for the celltype percentage estimates ("estimates"),
#   "params", which is the parameter set used for this run, and "markers", which
#   is the list of genes from the CibersortX signature matrix
CibersortX_InnerLoop <- function(reference_filename, bulk_mat, cibersort_signature, params) {
  set_cibersortx_credentials(email = Sys.getenv("CIBERSORT_EMAIL"),
                             token = Sys.getenv("CIBERSORT_TOKEN"))

  # Ensure signature is a matrix. For CibersortX, we do not need to filter the
  # signature down or worry about markers so we ignore those settings from 'params'
  cibersort_signature <- as.matrix(cibersort_signature)

  res_pcts <- deconvolute_cibersortx(bulk_mat, cibersort_signature,
                                     single_cell_object = basename(reference_filename),
                                     cell_type_annotations = colnames(cibersort_signature), # dummy variable, not used if a filename is provided but has to have a non-null value
                                     rmbatch_S_mode = TRUE,
                                     verbose = TRUE, container = "docker",
                                     input_dir = dir_cibersort,
                                     output_dir = dir_cibersort,
                                     qn = FALSE,
                                     absolute = FALSE)

  # Cleanup finished docker container
  system("docker rm $(docker ps -a -q --filter ancestor=cibersortx/fractions)")
  gc()

  res <- list("estimates" = res_pcts,
              "params" = params,
              "markers" = rownames(cibersort_signature))

  return(res)
}
