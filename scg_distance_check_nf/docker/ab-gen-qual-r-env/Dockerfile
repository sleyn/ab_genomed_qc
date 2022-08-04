FROM fedora:36

LABEL author "Semion Leyn"
LABEL email "semen.leyn@gmail.com"
LABEL description "Setup R environment for A.baumannii genomes analysis"

WORKDIR /home

RUN dnf -y install R

RUN R -e "install.packages(c('seqinr'), dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages(c('tidyverse'), dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages(c('hexbin'), dependencies=TRUE, repos='http://cran.rstudio.com/')"
