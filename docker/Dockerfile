FROM rocker/rstudio:4.2.1

RUN apt-get update -y && apt-get upgrade -y
RUN apt-get install -y dpkg-dev zlib1g-dev libssl-dev libffi-dev
RUN apt-get install -y curl libcurl4-openssl-dev
RUN apt-get install -y build-essential nano less libxt-dev libxml2 liblzma-dev
RUN apt-get install -y libigraph0v5 libgeos-dev libhdf5-dev libbz2-dev

RUN R -e "install.packages(c('BiocManager', 'devtools', 'stringr', 'Metrics'))"
RUN R -e "install.packages(c('dtangle', 'readxl'))"
RUN R -e "install.packages('synapser', repos = c('http://ran.synapse.org', 'http://cran.fhcrc.org'))"
RUN R -e "BiocManager::install(c('Biobase', 'SingleCellExperiment', 'TOAST'))"
RUN R -e "BiocManager::install(c('scuttle', 'DeconRNASeq', 'Seurat', 'biomaRt'))"
RUN R -e "BiocManager::install(c('GEOquery', 'edgeR', 'rhdf5', 'HDF5Array'))"
RUN R -e "BiocManager::install('preprocessCore', configure.args='--disable-threading')"

RUN R -e "devtools::install_github('xuranw/MuSiC')"
RUN R -e "devtools::install_github('gjhunt/hspe/lib_hspe')"
