FROM rocker/r-base:4.0.4

MAINTAINER Julian Lucas, juklucas@ucsc.edu



## Install NTRprism
WORKDIR /opt/

RUN wget https://github.com/altemose/NTRprism/archive/refs/tags/v1.0.0.tar.gz \
	&& tar -xvf v1.0.0.tar.gz \
	&& cd NTRprism-1.0.0 \
	&& chmod +x NTRprism_PlotHeatmap.r NTRprism_PlotSpectrum.r NTRprism_ProcessFasta_v0.21.pl

ENV PATH=/opt/NTRprism-1.0.0:$PATH

RUN R -e "install.packages(c('ggplot2', 'reshape2'))"

WORKDIR /data
