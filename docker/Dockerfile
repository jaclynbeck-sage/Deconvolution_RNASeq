FROM rocker/rstudio

RUN apt-get update -y && apt-get upgrade -y
RUN apt-get install -y dpkg-dev zlib1g-dev libssl-dev libffi-dev ssh
RUN apt-get install -y curl libcurl4-openssl-dev
RUN apt-get install -y build-essential nano less libxt-dev libxml2 liblzma-dev
RUN apt-get install -y libglpk-dev libgeos-dev libhdf5-dev libbz2-dev

USER rstudio

RUN cd /home/rstudio && git clone https://github.com/jaclynbeck-sage/Deconvolution_RNASeq.git
RUN cd /home/rstudio/Deconvolution_RNASeq && R -e "source('Step00_InitialSetupInstall.R')"

USER root
