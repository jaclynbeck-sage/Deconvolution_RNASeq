install.packages(
  c(
    "BiocManager", "remotes", "reticulate", "foreach", "doParallel", "config",  # general
    "stringr", "tidyr", "purrr", "dplyr", "ggplot2", "readxl", "vroom",  # tidyverse
    "anndata", "hdf5r",  # h5ad ingest / AutogeneS
    "Metrics", "Hmisc", "lme4", "DescTools",  # analysis
    "uuid",  # omnideconv dependency
    "RColorBrewer", "viridis", "patchwork" # plotting
  ),
  dependencies = c("Depends", "Imports"),
  clean = TRUE
)

# Extra setup to get synapser to install -- Downgrade to Python 3.11 and install
# rjson version 0.2.21. Freezing at specific Python version 3.11.12 which is
# the version used to run data for the paper.
reticulate::install_python(version = "3.11.12")
reticulate::virtualenv_create("r-reticulate")

remotes::install_version("rjson", version = "0.2.21")
install.packages("synapser",
                 repos = c("http://ran.synapse.org", "https://cloud.r-project.org"),
                 dependencies = c("Depends", "Imports"),
                 clean = TRUE)

BiocManager::install(
  c(
    "rtracklayer", "BSgenome", # gtf file ingest
    "SingleCellExperiment", "Seurat", "edgeR", "scuttle", "DESeq2",  # general
    "GEOquery", "scDblFinder", "HDF5Array",  # data ingest / QC
    "variancePartition", "sva",  # regression
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
remotes::install_github("GabrielHoffman/mvIC", clean = TRUE)

# install the MuSiC package from my fork, which has a few fixes and speedups
remotes::install_github("jaclynbeck-sage/MuSiC", clean = TRUE)

# install Dtangle and HSPE from my forks, which have a few fixes
remotes::install_github("jaclynbeck-sage/dtangle", subdir = "lib_dtangle",
                        clean = TRUE)
remotes::install_github("jaclynbeck-sage/hspe", subdir = "lib_hspe",
                        clean = TRUE)

# Extra package needed for pre-processing seaRef
reticulate::virtualenv_install("r-reticulate", packages = c("anndata"))

# Install omnideconv and its DWLS add-on. Uses my fork of omnideconv
remotes::install_github("jaclynbeck-sage/omnideconv", ref = "use_virtualenv",
                        clean = TRUE)
remotes::install_github("omnideconv/DWLS",
                        clean = TRUE)

# Ensure that pip manages all dependencies at once for all packages that need to
# get installed for omnideconv. Calling omnideconv::install_all_python() does
# the installation piecemeal and results in an incompatible environment.
#reticulate::virtualenv_create("r-omnideconv", version = "3.11.12", force = TRUE,
#                              packages = c("numpy", "anndata", "scanpy",
#                                           "scikit-learn", "scikit-misc",
#                                           "git+https://github.com/omnideconv/AutoGeneS.git",
#                                           "git+https://github.com/omnideconv/scaden.git"))

# After running the above virtualenv_create and confirming that the resulting
# environment works, I used `pip freeze > omnideconv_requirements.txt` to get a
# requirements file. This will ensure that the exact package versions used are
# always installed in the docker container.
reticulate::virtualenv_create("r-omnideconv",
                              requirements = file.path("docker", "omnideconv_requirements.txt"))
