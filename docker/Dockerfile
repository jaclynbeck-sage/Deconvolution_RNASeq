FROM rocker/rstudio:4.3.2

RUN apt-get update -y && apt-get upgrade -y
RUN apt-get install -y ssh nano less python3 python3-pip python3-venv
RUN apt-get install -y libpng-dev libglpk-dev curl libxt-dev
RUN apt-get install -y liblzma-dev libbz2-dev libhdf5-dev
RUN apt-get install -y libmysqlclient21 # for GenomicFeatures
RUN apt-get clean

# Install docker for CIBERSORTx
RUN cd /home/rstudio && curl -fsSL https://get.docker.com -o install-docker.sh
RUN cd /home/rstudio && sh install-docker.sh

USER rstudio

RUN cd /home/rstudio && git clone https://github.com/jaclynbeck-sage/Deconvolution_RNASeq.git
RUN cd /home/rstudio/Deconvolution_RNASeq && R -e "source('Step00_InitialSetupInstall.R')"

USER root
