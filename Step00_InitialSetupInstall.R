install.packages(c("BiocManager", "devtools", "stringr", "Metrics",
                   "dplyr", "ggplot2", "DEoptimR", "nloptr", "anndata",
                   "readxl", "Hmisc", "RMariaDB", "lme4", "corrplot",
                   "circlize", "gMWT", "pak", "uuid"),
                 clean = TRUE)

# If synapser fails to install because it can't find "synapseclient", go to
# RStudio options (Tools->Global Options) -> Python, uncheck "Automatically
# activate project-local Python environments" and restart R.
install.packages("synapser", repos = c("http://ran.synapse.org", "http://cran.fhcrc.org"))

BiocManager::install(c("Biobase", "SingleCellExperiment", "TOAST", "scuttle",
                       "DeconRNASeq", "MAST", "GEOquery", "biomaRt",
                       "DESeq2", "edgeR", "GenomicFeatures", "snpStats",
                       "HDF5Array", "glmGamPoi"))

# Install version 4 of Seurat instead of version 5. SeuratObject v5 will be installed
# but it will work with Seurat v4. Installing SeuratObject v4 caused install issues.
remotes::install_version("Seurat", "4.4.0",
                         repos = c("https://satijalab.r-universe.dev", getOption("repos")))

BiocManager::install("preprocessCore", configure.args="--disable-threading")

# These packages require several BiocManager packages to be installed first
install.packages(c("GenomicTools.fileHandler"))

# install the MuSiC package from my fork, which has a few fixes and speedups
devtools::install_github("jaclynbeck-sage/MuSiC")

# install Dtangle and HSPE from my forks, which have a few fixes
devtools::install_github("jaclynbeck-sage/dtangle", subdir = "lib_dtangle")
devtools::install_github("jaclynbeck-sage/hspe", subdir = "lib_hspe")

# install sageseqr -- uses my fork with updated package dependencies
devtools::install_github("jaclynbeck-sage/sageseqr")

# Extra package needed for pre-processing seaRef
reticulate::virtualenv_install("r-reticulate", packages = c("anndata"))

# virtualenv setup for AutogeneS
reticulate::virtualenv_create("autogenes_env",
                              packages = c("scanpy", "anndata", "autogenes",
                                           "scikit-misc"))

# Install omnideconv and its DWLS add-on
pak::pkg_install("omnideconv/omnideconv")
pak::pkg_install("omnideconv/DWLS")
