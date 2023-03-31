install.packages(c("BiocManager", "devtools", "stringr", "Metrics",
                   "reticulate", "dplyr", "ggplot2",
                   "DEoptimR", "nloptr"))

# If synapser fails to install because it can't find "synapseclient", go to
# RStudio options (Tools->Global Options) -> Python, uncheck "Automatically
# activate project-local Python environments" and restart R.
install.packages("synapser", repos = c("http://ran.synapse.org", "http://cran.fhcrc.org"))

BiocManager::install(c("Biobase", "SingleCellExperiment", "scuttle",
                       "DeconRNASeq", "Seurat", "MAST", "GEOquery",
                       "rhdf5", "HDF5Array", "zellkonverter"))

# install the MuSiC package
devtools::install_github('xuranw/MuSiC')

# HSPE doesn't provide support for sparse matrices, so I download the repo and
# make a few small changes. It's then re-compiled into a package and installed
# as "hspeSparse".
# NOTE: This may break if the developer updates their repo
system("bash install_hspe.sh")
install.packages(file.path('~', 'hspe', 'hspeSparse_0.1.tar.gz'),
                 repos = NULL, type = "source")

# Same thing with dtangle
system("bash install_dtangle.sh")
install.packages(file.path('~', 'dtangle', 'dtangleSparse_2.0.9.tar.gz'),
                 repos = NULL, type = "source")

# Extra package needed for pre-processing seaRef
reticulate::virtualenv_install("r-reticulate", packages = c("anndata"))

# virtualenv setup for AutogeneS
reticulate::virtualenv_create("autogenes_env",
                              packages = c("scanpy", "anndata", "autogenes",
                                           "scikit-misc"))
