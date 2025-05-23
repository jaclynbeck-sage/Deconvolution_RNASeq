install.packages(c("BiocManager", "devtools", "stringr", "Metrics", "dplyr",
                   "ggplot2", "anndata", "readxl", "Hmisc", "lme4",
                   "reticulate", "pak", "uuid"),
                 clean = TRUE)

# If synapser fails to install because it can't find "synapseclient", go to
# RStudio options (Tools->Global Options) -> Python, uncheck "Automatically
# activate project-local Python environments" and restart R.
reticulate::virtualenv_create("r-reticulate")
install.packages("synapser", repos = c("http://ran.synapse.org", "https://cloud.r-project.org"))

BiocManager::install(c("SingleCellExperiment", "TOAST", "scuttle", "DESeq2",
                       "edgeR", "DeconRNASeq", "Seurat", "GEOquery",
                       "HDF5Array", "scDblFinder"),
                     clean = TRUE)

#BiocManager::install(c("Biobase", "MAST", "GenomicFeatures", "snpStats", "glmGamPoi"))
#BiocManager::install("preprocessCore", configure.args="--disable-threading")

# install my utility package
devtools::install_github("jaclynbeck-sage/sageRNAUtils", clean = TRUE)

# install the MuSiC package from my fork, which has a few fixes and speedups
devtools::install_github("jaclynbeck-sage/MuSiC", clean = TRUE)

# install Dtangle and HSPE from my forks, which have a few fixes
devtools::install_github("jaclynbeck-sage/dtangle", subdir = "lib_dtangle", clean = TRUE)
devtools::install_github("jaclynbeck-sage/hspe", subdir = "lib_hspe", clean = TRUE)

# Extra package needed for pre-processing seaRef
reticulate::virtualenv_install("r-reticulate", packages = c("anndata"))

# virtualenv setup for AutogeneS
#reticulate::virtualenv_create("autogenes_env",
#                              packages = c("scanpy", "anndata", "autogenes",
#                                           "scikit-misc"))

# Install omnideconv and its DWLS add-on. Uses my fork of omnideconv
pak::pkg_install("jaclynbeck-sage/omnideconv")
pak::pkg_install("omnideconv/DWLS")
omnideconv::install_all_python()
