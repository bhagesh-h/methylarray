#
# methylarray developer toolbox.
#
# Provides Nextflow, nf-core tools, nf-test, nf-metro, prek and a Docker client
# so that the pipeline can be developed and tested on a machine where only
# Docker is installed. End users do not need this image: they run the pipeline
# with their own Nextflow installation.

FROM docker:29.9.3-cli AS docker-cli

FROM mambaorg/micromamba:2.10.0-ubuntu24.04

ARG NXF_VER=26.10.6
ARG NXF_SHA256=80e5332f37cd9a3b328fb843271209e090d91a7813f238e92f4efe3b5ee5d8d9

LABEL org.opencontainers.image.title="methylarray-devkit" \
      org.opencontainers.image.description="Developer toolbox for the methylarray Nextflow pipeline" \
      org.opencontainers.image.source="https://github.com/bhagesh-h/methylarray" \
      org.opencontainers.image.licenses="MIT"

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install --yes --name base --file /tmp/environment.yml \
 && micromamba clean --all --yes \
 && rm /tmp/environment.yml

USER root
ENV PATH=/opt/conda/bin:${PATH}

# Static Docker client, used to drive the host daemon (Docker outside of Docker)
# or a Docker-in-Docker sidecar.
COPY --from=docker-cli /usr/local/bin/docker /usr/local/bin/docker

# Self-contained Nextflow distribution, verified against its SHA-256 checksum.
# The conda launcher pulled in by nf-core and nf-test is replaced by a link to it,
# because the micromamba entrypoint puts /opt/conda/bin first on PATH.
# Downloaded with the environment's Python so that the classic builder works
# (no BuildKit or buildx needed on the host).
RUN python -c "import urllib.request as u; u.urlretrieve('https://github.com/nextflow-io/nextflow/releases/download/v${NXF_VER}/nextflow-${NXF_VER}-dist', '/usr/local/bin/nextflow')" \
 && echo "${NXF_SHA256}  /usr/local/bin/nextflow" | sha256sum --check --strict \
 && chmod 0755 /usr/local/bin/nextflow \
 && ln -sf /usr/local/bin/nextflow /opt/conda/bin/nextflow

ENV NXF_VER=${NXF_VER} \
    NXF_ANSI_LOG=false

# Fail the build if any tool is missing or cannot start.
RUN test "$(readlink -f "$(command -v nextflow)")" = /usr/local/bin/nextflow \
 && NXF_HOME=/tmp/nxf-check HOME=/tmp nextflow -version \
 && nf-core --version \
 && nf-test version \
 && nf-metro --version \
 && prek --version \
 && docker --version \
 && java -version \
 && rm -rf /tmp/nxf-check

USER $MAMBA_USER
WORKDIR /tmp
