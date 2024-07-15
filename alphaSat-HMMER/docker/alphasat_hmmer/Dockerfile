FROM ubuntu:18.04

MAINTAINER Julian Lucas, juklucas@ucsc.edu


RUN apt-get update && apt-get install -y \
    python-pip \
    python \
    python3.5 \
    wget \
    git \
    pigz


## Install HMMER (need nhmmer)
WORKDIR /opt/hmmer

RUN wget http://eddylab.org/software/hmmer/hmmer-3.2.1.tar.gz \
	&& tar zxf hmmer-3.2.1.tar.gz \
	&& cd hmmer-3.2.1 \
	&& ./configure --prefix /usr/local \
	&& make \
	&& make check \
	&& make install \
	&& (cd easel; make install) \
	&& (cd ..; rm hmmer-3.2.1.tar.gz)


## Install BEDOPS (need bedmap)
WORKDIR /opt/bedops

RUN wget https://github.com/bedops/bedops/releases/download/v2.4.40/bedops_linux_x86_64-v2.4.40.tar.bz2 \
	&& tar jxvf bedops_linux_x86_64-v2.4.40.tar.bz2 \
	&& cp bin/* /usr/local/bin \
	&& rm -rf bin \
	&& rm bedops_linux_x86_64-v2.4.40.tar.bz2


## Install HumAS-HMMER for AnVIL
WORKDIR /opt/

RUN git clone https://github.com/fedorrik/HumAS-HMMER_for_AnVIL.git \
	&& cd HumAS-HMMER_for_AnVIL \
	&& git reset --hard f7e1e13b0d63e2007b53d90b274875a079750ee2 \
	&& chmod +x hmmer-run.sh hmmer-run_SF.sh hmmertblout2bed.awk overlap_filter.py

ENV PATH=/opt/HumAS-HMMER_for_AnVIL:$PATH

WORKDIR /data