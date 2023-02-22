install.packages(c("BiocManager", "devtools", "stringr", "Metrics",
                   "readxl"))

# requires install of several linux packages first
install.packages("synapser", repos = c("http://ran.synapse.org", "http://cran.fhcrc.org"))

BiocManager::install(c("Biobase", "SingleCellExperiment", "TOAST", "scuttle",
                       "DeconRNASeq", "Seurat", "GEOquery", "rhdf5", "HDF5Array"))

# install the MuSiC package
devtools::install_github('xuranw/MuSiC')

# install Spacexr / RCTD package
#devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)

# Requires install of libxml on Linux
BiocManager::install(c("biomaRt"))

# HSPE
#devtools::install_github("gjhunt/hspe/lib_hspe")
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

