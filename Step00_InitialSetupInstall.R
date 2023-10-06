install.packages(c("BiocManager", "devtools", "stringr", "Metrics",
                   "dplyr", "ggplot2", "DEoptimR", "nloptr", "anndata",
                   "readxl", "Hmisc", "RMariaDB", "lme4", "corrplot",
                   "Hmisc", "circlize", "gMWT", "DWLS"))

# If synapser fails to install because it can't find "synapseclient", go to
# RStudio options (Tools->Global Options) -> Python, uncheck "Automatically
# activate project-local Python environments" and restart R.
install.packages("synapser", repos = c("http://ran.synapse.org", "http://cran.fhcrc.org"))

BiocManager::install(c("Biobase", "SingleCellExperiment", "TOAST", "scuttle",
                       "DeconRNASeq", "Seurat", "MAST", "GEOquery", "biomaRt",
                       "DESeq2", "edgeR", "GenomicFeatures")) #,
                       #"rhdf5", "HDF5Array"))

BiocManager::install("preprocessCore", configure.args="--disable-threading")

# install the MuSiC package from my fork, which has a few fixes and speedups
devtools::install_github("jaclynbeck-sage/MuSiC")

# install sageseqr -- requires archived package from CRAN
gt_url <- "https://cran.r-project.org/src/contrib/Archive/GenomicTools/GenomicTools_0.2.9.7.tar.gz"
download.file(url = gt_url, destfile = "GenomicTools_0.2.9.7.tar")
install.packages(pkgs="GenomicTools_0.2.9.7.tar", type="source", repos=NULL)
unlink("GenomicTools_0.2.9.7.tar") # deletes the tar file

devtools::install_github("Sage-Bionetworks/sageseqr")

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

#install.packages("pak")
#pak::pkg_install("omnideconv/omnideconv", dependencies = c("Imports", "Remotes"))
