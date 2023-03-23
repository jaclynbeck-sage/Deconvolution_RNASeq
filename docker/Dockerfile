FROM rocker/rstudio:4.2.1

RUN apt-get update -y && apt-get upgrade -y
RUN apt-get install -y dpkg-dev zlib1g-dev libssl-dev libffi-dev ssh
RUN apt-get install -y curl libcurl4-openssl-dev
RUN apt-get install -y build-essential nano less libxt-dev libxml2 liblzma-dev
RUN apt-get install -y libigraph0v5 libgeos-dev libhdf5-dev libbz2-dev

USER rstudio

RUN R -e "install.packages(c('BiocManager', 'devtools', 'stringr', 'Metrics'))"
RUN R -e "install.packages(c('readxl', 'reticulate', 'DEoptimR', 'nloptr'))"
RUN R -e "reticulate::install_miniconda()"
RUN R -e "install.packages('synapser', repos = c('http://ran.synapse.org', 'http://cloud.r-project.org'))"
RUN R -e "BiocManager::install(c('Biobase', 'SingleCellExperiment', 'scuttle'))"
RUN R -e "BiocManager::install(c('TOAST', 'DeconRNASeq', 'Seurat', 'biomaRt'))"
RUN R -e "BiocManager::install(c('GEOquery', 'rhdf5', 'HDF5Array'))"

RUN R -e "devtools::install_github('xuranw/MuSiC')"

RUN cd /home/rstudio && git clone https://github.com/jaclynbeck-sage/Deconvolution_RNASeq.git
RUN cd /home/rstudio/Deconvolution_RNASeq && bash install_hspe.sh && bash install_dtangle.sh

RUN R -e "install.packages(file.path('~', 'hspe', 'hspeSparse_0.1.tar.gz'), repos = NULL, type = 'source')"
RUN R -e "install.packages(file.path('~', 'dtangle', 'dtangleSparse_2.0.9.tar.gz'), repos = NULL, type = 'source')"

RUN R -e "reticulate::conda_create(envname = 'autogenes_env', packages = c('python=3.10', 'libffi=3.3'))"

pkgs="c('scanpy', 'pandas', 'numpy', 'scipy', 'anndata', 'anndata2ri', 'autogenes', 'rpy2', 'scikit-misc')"
RUN R -e "reticulate::conda_install(envname = 'autogenes_env', packages = $pkgs, pip = TRUE)"

USER root
