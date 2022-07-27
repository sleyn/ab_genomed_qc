FROM ubuntu:20.04 AS build

LABEL author "Semion Leyn"
LABEL email "semen.leyn@gmail.com"
LABEL description "MAFFT 5.05 aligner"

WORKDIR /home

RUN apt-get update \
    && apt-get -y install wget \
                          make \
                          gcc \
                          g++ \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN wget https://mafft.cbrc.jp/alignment/software/mafft-7.505-with-extensions-src.tgz && \
    tar -xvzf mafft-7.505-with-extensions-src.tgz && \
    cd mafft-7.505-with-extensions/core && \
    make clean && \
    make && \
    make install && \
    cd ../extensions && \
    make clean && \
    make && \
    make install

#############################################################################
# Configure main release
FROM ubuntu:20.04
WORKDIR /root/
COPY --from=build /usr/local/libexec/mafft /usr/local/libexec/mafft
COPY --from=build /usr/local/bin/ /usr/local/bin/
COPY --from=build /usr/local/share/man/man1/mafft* /usr/local/share/man/man1
