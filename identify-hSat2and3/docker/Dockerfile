FROM ubuntu:18.04

MAINTAINER Julian Lucas, juklucas@ucsc.edu


RUN apt-get update && apt-get install -y \
    git


## Install chm13_hsat
WORKDIR /opt/

RUN git clone https://github.com/altemose/chm13_hsat.git \
	&& cd chm13_hsat \
	&& git reset --hard 2e8c31a93cb73400037e45151cf4023a0a3706d7 \
	&& chmod +x Assembly_HSat2and3_v3.pl

ENV PATH=/opt/chm13_hsat:$PATH

WORKDIR /data
