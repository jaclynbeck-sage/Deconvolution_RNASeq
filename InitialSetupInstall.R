install.packages(c("BiocManager", "devtools", "stringr", "Metrics"))

# requires install of several linux packages first
install.packages("synapser", repos = c("http://ran.synapse.org", "http://cran.fhcrc.org"))

install.packages(c("dtangle"))

BiocManager::install(c("Biobase", "SingleCellExperiment", "TOAST", "scuttle",
                       "DeconRNASeq"))
                       #"edgeR", "DelayedArray", "rhdf5", "HDF5Array", "scuttle"))

BiocManager::install("preprocessCore", configure.args="--disable-threading")

# install the MuSiC package
devtools::install_github('xuranw/MuSiC')

# install Spacexr / RCTD package
#devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)

# Temporary? Install biomart. Requires install of libxml on Linux
BiocManager::install(c("biomart"))

# HSPE
devtools::install_github("gjhunt/hspe/lib_hspe")


