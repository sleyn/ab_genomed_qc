FROM ubuntu:20.04

LABEL author "Semion Leyn"
LABEL email "semen.leyn@gmail.com"
LABEL description "Setup Python 3 environment for A.baumannii genomes analysis"

WORKDIR /home

RUN apt-get update \
    && apt-get -y install python3 \
                          python3-pip \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN pip3 --no-cache-dir install biopython pandas
