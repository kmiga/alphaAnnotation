FROM rocker/r-base:4.0.4

MAINTAINER Julian Lucas, juklucas@ucsc.edu


RUN apt-get update && \
    apt-get install -y git


## Install NTRprism
WORKDIR /opt/

RUN git clone https://github.com/altemose/NTRprism.git && \
    cd NTRprism && \
    git reset --hard 438fb135cc7dca4d19ab636ba97d94fbe544feaf && \
	chmod +x NTRprism_PlotHeatmap.r NTRprism_PlotSpectrum.r NTRprism_ProcessFasta_v0.22.pl

ENV PATH=/opt:$PATH


RUN R -e "install.packages(c('ggplot2', 'reshape2'))"

WORKDIR /data
