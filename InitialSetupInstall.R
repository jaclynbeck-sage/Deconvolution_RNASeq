install.packages(c("BiocManager", "devtools", "stringr", "Metrics",
                   "dtangle", "readxl"))

# requires install of several linux packages first
install.packages("synapser", repos = c("http://ran.synapse.org", "http://cran.fhcrc.org"))

BiocManager::install(c("Biobase", "SingleCellExperiment", "TOAST", "scuttle",
                       "DeconRNASeq", "Seurat", "GEOquery", "edgeR",
                       "rhdf5", "HDF5Array"))

BiocManager::install("preprocessCore", configure.args="--disable-threading")

# install the MuSiC package
devtools::install_github('xuranw/MuSiC')

# install Spacexr / RCTD package
#devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)

# Requires install of libxml on Linux
BiocManager::install(c("biomaRt"))

# HSPE
devtools::install_github("gjhunt/hspe/lib_hspe")
