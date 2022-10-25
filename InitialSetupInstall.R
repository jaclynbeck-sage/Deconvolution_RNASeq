install.packages(c("BiocManager", "devtools", "stringr"))

# requires install of several linux packages first
install.packages("synapser", repos = c("http://ran.synapse.org", "http://cran.fhcrc.org"))

install.packages(c("dtangle"))

BiocManager::install(c("Biobase", "SingleCellExperiment", "TOAST"))
                       #"edgeR", "DelayedArray", "rhdf5", "HDF5Array", "scuttle"))

# install the MuSiC package
devtools::install_github('xuranw/MuSiC')

# install Spacexr / RCTD package
#devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)

# Temporary? Install biomart. Requires install of libxml on Linux
BiocManager::install(c("biomart"))



