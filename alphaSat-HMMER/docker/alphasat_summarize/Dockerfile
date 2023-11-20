FROM ubuntu:18.04

MAINTAINER Julian Lucas, juklucas@ucsc.edu


RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    wget \
    git \
    pigz \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev


RUN pip3 install --no-cache-dir --upgrade pip \
    && pip3 install pandas pybedtools

## Install BEDTOOLS
WORKDIR /opt/

RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary \
    && mv bedtools.static.binary bedtools \
    && chmod a+x bedtools

ENV PATH=/opt:$PATH


## Copy in scripts that aggregate HMM monomers from HumAS-HMMER into HORs
RUN mkdir -p /opt/scripts
COPY ./scripts /opt/scripts/ 


WORKDIR /data
