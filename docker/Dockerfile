FROM rocker/rstudio:4.2.3

RUN apt-get update -y && apt-get upgrade -y
RUN apt-get install -y ssh nano less python3 python3-pip python3-venv
RUN apt-get install -y libpng-dev libglpk-dev curl libxt-dev
RUN apt-get install -y liblzma-dev libbz2-dev libmysqlclient21 # for GenomicFeatures

# TODO Remove libopenblas-dev, libopenblas-pthread-dev, replace with libblas-dev and liblapack-dev?
# Unclear if these packages are still needed
#RUN apt-get install -y libffi-dev libgeos-dev libhdf5-dev

USER rstudio

RUN cd /home/rstudio && git clone https://github.com/jaclynbeck-sage/Deconvolution_RNASeq.git
RUN cd /home/rstudio/Deconvolution_RNASeq && R -e "source('Step00_InitialSetupInstall.R')"

USER root
