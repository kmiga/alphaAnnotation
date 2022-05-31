FROM ubuntu:18.04

MAINTAINER Julian Lucas, juklucas@ucsc.edu


RUN apt-get update && apt-get install -y \
    python3 \
    git \
    bedtools
    

## Install HOR StV (Structural Variant)
WORKDIR /opt/

RUN git clone https://github.com/fedorrik/stv.git \
	&& cd stv \
	&& git reset --hard c11f823cb74e286c2f1bccd69bd34ec14f5cfb9a \
	&& chmod +x scripts/* \
	&& sed -i 's/\/home\/fedor\/Programs\/my/\/opt/' stv.sh

ENV PATH=/opt/stv:$PATH

WORKDIR /data