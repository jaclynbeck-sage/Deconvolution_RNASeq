FROM rocker/rstudio:4.4

RUN apt-get update -y && apt-get upgrade -y && apt-get clean
RUN apt-get install -y ssh nano less curl \
  libpng-dev libglpk-dev liblzma-dev libbz2-dev libffi-dev libhdf5-dev \
  libcurl4-openssl-dev patch \
  && apt-get clean
# libcurl4-openssl-dev required for scDblFinder
# patch required for hdf5r

# Install docker for CIBERSORTx
RUN cd /home/rstudio && curl -fsSL https://get.docker.com -o install-docker.sh
RUN cd /home/rstudio && sh install-docker.sh

USER rstudio

RUN cd /home/rstudio && git clone https://github.com/jaclynbeck-sage/Deconvolution_RNASeq.git \
  && cd Deconvolution_RNASeq && git checkout sc_qc
RUN cd /home/rstudio/Deconvolution_RNASeq && R -e "source('Step00_InitialSetupInstall.R')" \
  && /home/rstudio/.virtualenvs/r-omnideconv/bin/pip cache purge

USER root
