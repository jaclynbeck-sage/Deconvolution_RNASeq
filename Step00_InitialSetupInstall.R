install.packages(
  c(
    "BiocManager", "remotes", "reticulate", "foreach", "doParallel",  # general
    "stringr", "tidyr", "purrr", "dplyr", "ggplot2", "readxl", "vroom",  # tidyverse
    "anndata",  # h5ad ingest / AutogeneS
    "Metrics", "Hmisc", "lme4",  # analysis
    "uuid",  # omnideconv dependency
    "RColorBrewer", "viridis", "patchwork" # plotting
  ),
  clean = TRUE
)

# Extra setup to get synapser to install -- Downgrade to Python 3.11 and install
# rjson version 0.2.21
reticulate::install_python(version = "3.11")
reticulate::virtualenv_create("r-reticulate")

remotes::install_version("rjson", version = "0.2.21")
install.packages("synapser", repos = c("http://ran.synapse.org", "https://cloud.r-project.org"))

BiocManager::install(
  c(
    "rtracklayer", "BSgenome", # gtf file ingest
    "SingleCellExperiment", "Seurat", "edgeR", "scuttle", "DESeq2",  # general
    "GEOquery", "scDblFinder", "HDF5Array",  # data ingest / QC
    "variancePartition",  # regression
    "TOAST",  # MuSiC dependency
    "MAST",  # omnideconv dependency
    "DeconRNASeq"  # deconvolution
  ),
  clean = TRUE
)

#BiocManager::install(c("Biobase", "GenomicFeatures", "snpStats", "glmGamPoi"))
#BiocManager::install("preprocessCore", configure.args="--disable-threading")

# install presto for faster marker finding in Seurat
remotes::install_github("immunogenomics/presto", clean = TRUE)

# install my utility package
remotes::install_github("jaclynbeck-sage/sageRNAUtils", clean = TRUE)

# install mvIC for regression
remotes::install_github("GabrielHoffman/mvIC")

# install the MuSiC package from my fork, which has a few fixes and speedups
remotes::install_github("jaclynbeck-sage/MuSiC", clean = TRUE)

# install Dtangle and HSPE from my forks, which have a few fixes
remotes::install_github("jaclynbeck-sage/dtangle", subdir = "lib_dtangle", clean = TRUE)
remotes::install_github("jaclynbeck-sage/hspe", subdir = "lib_hspe", clean = TRUE)

# Extra package needed for pre-processing seaRef
reticulate::virtualenv_install("r-reticulate", packages = c("anndata"))

# Install omnideconv and its DWLS add-on. Uses my fork of omnideconv
remotes::install_github("jaclynbeck-sage/omnideconv", ref = "use_virtualenv", clean = TRUE)
remotes::install_github("omnideconv/DWLS", clean = TRUE)
omnideconv::install_all_python()

# virtualenv setup for AutogeneS marker finding
reticulate::virtualenv_install("r-omnideconv", packages = c("scanpy", "scikit-misc"))
