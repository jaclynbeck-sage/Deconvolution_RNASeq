FROM rocker/rstudio:4.2.1

RUN apt-get update -y && apt-get upgrade -y
RUN apt-get install -y dpkg-dev zlib1g-dev libssl-dev libffi-dev ssh
RUN apt-get install -y curl libcurl4-openssl-dev
RUN apt-get install -y build-essential nano less libxt-dev libxml2 liblzma-dev
RUN apt-get install -y libigraph0v5 libgeos-dev libhdf5-dev libbz2-dev

RUN cd ~ && git clone https://github.com/jaclynbeck-sage/Deconvolution_RNASeq.git

RUN R -e "install.packages(c('BiocManager', 'devtools', 'stringr', 'Metrics'))"
RUN R -e "install.packages(c('readxl'))"
RUN R -e "install.packages('synapser', repos = c('http://ran.synapse.org', 'http://cran.fhcrc.org'))"
RUN R -e "BiocManager::install(c('Biobase', 'SingleCellExperiment', 'scuttle'))"
RUN R -e "BiocManager::install(c('TOAST', 'DeconRNASeq', 'Seurat', 'biomaRt'))"
RUN R -e "BiocManager::install(c('GEOquery', 'rhdf5', 'HDF5Array'))"

RUN R -e "devtools::install_github('xuranw/MuSiC')"

RUN cd ~/Deconvolution_RNASeq && bash install_hspe.sh && bash install_dtangle.sh

RUN R -e "install.packages(file.path('~', 'hspe', 'hspeSparse_0.1.tar.gz'), repos = NULL, type = 'source')"
RUN R -e "install.packages(file.path('~', 'dtangle', 'dtangleSparse_2.0.9.tar.gz'), repos = NULL, type = 'source')"

