# Setup all tools for A.baumannii genomes analysis
FROM mambaorg/micromamba:0.24.0
COPY --chown=$MAMBA_USER:$MAMBA_USER ./scripts/environment.yml /apps/environment.yml

RUN micromamba install -y -n base -f /apps/environment.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

WORKDIR /apps/
COPY scripts/patric2fasta.py /apps
ENV PATH=/apps:${PATH}

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]