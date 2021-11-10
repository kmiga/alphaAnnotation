FROM ubuntu:18.04

MAINTAINER Julian Lucas, juklucas@ucsc.edu


RUN apt-get update && apt-get install -y \
    git


## Install chm13_hsat
WORKDIR /opt/

RUN git clone https://github.com/altemose/chm13_hsat.git \
	&& cd chm13_hsat \
	&& git reset --hard 1096d62a929fecc0ede419fb4141d17441b4283f \
	&& chmod +x Assembly_HSat2and3_v2.pl

ENV PATH=/opt/chm13_hsat:$PATH

WORKDIR /data
