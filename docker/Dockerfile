FROM rocker/rstudio:4.4

RUN apt-get update -y && apt-get upgrade -y
RUN apt-get install -y ssh nano less curl
RUN apt-get install -y libpng-dev libglpk-dev liblzma-dev libbz2-dev libffi-dev libhdf5-dev
RUN apt-get install -y libcurl4-openssl-dev # for scDblFinder
RUN apt-get install -y patch # For hdf5r
#RUN apt-get install -y libxt-dev
RUN apt-get clean

# Install docker for CIBERSORTx
RUN cd /home/rstudio && curl -fsSL https://get.docker.com -o install-docker.sh
RUN cd /home/rstudio && sh install-docker.sh

USER rstudio

RUN cd /home/rstudio && git clone https://github.com/jaclynbeck-sage/Deconvolution_RNASeq.git
RUN cd /home/rstudio/Deconvolution_RNASeq && git checkout sc_qc
RUN cd /home/rstudio/Deconvolution_RNASeq && R -e "source('Step00_InitialSetupInstall.R')"

USER root
