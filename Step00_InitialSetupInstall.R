install.packages(c("BiocManager", "devtools", "stringr", "Metrics",
                   "readxl", "reticulate"))

reticulate::install_miniconda()

# requires install of several linux packages first
install.packages("synapser", repos = c("http://ran.synapse.org", "http://cran.fhcrc.org"))

BiocManager::install(c("Biobase", "SingleCellExperiment", "TOAST", "scuttle",
                       "DeconRNASeq", "Seurat", "GEOquery", "biomaRt",
                       "rhdf5", "HDF5Array"))

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

# Conda setup for AutogeneS
reticulate::conda_create(envname = "autogenes_env",
                         packages = c("python=3.10", "libffi=3.3"))
reticulate::conda_install(envname = "autogenes_env",
                          packages = c("scanpy", "pandas", "numpy", "scipy",
                                       "anndata", "anndata2ri", "autogenes",
                                       "rpy2", "scikit-misc"),
                          pip = TRUE) # Needs to be pip, not conda
