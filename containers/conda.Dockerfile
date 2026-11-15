#
# Task image built from one module's environment.yml, for modules whose image
# cannot come from Seqera Containers (Wave). The software comes only from the
# environment.yml, so the Docker, Singularity/Apptainer and conda profiles run
# the same package versions. Build from the repository root, for example:
#   dev/nfdev image modules/local/exploration/pca bhagesh/methylarray-exploration:1.0.0
# which runs:
#   docker build -f containers/conda.Dockerfile \
#       --build-arg ENVIRONMENT=modules/local/exploration/pca/environment.yml \
#       -t bhagesh/methylarray-exploration:1.0.0 .
#
FROM mambaorg/micromamba:2.10.0-ubuntu24.04

ARG ENVIRONMENT
COPY --chown=$MAMBA_USER:$MAMBA_USER ${ENVIRONMENT} /tmp/environment.yml
RUN micromamba install --yes --name base --file /tmp/environment.yml \
 && micromamba clean --all --yes \
 && rm /tmp/environment.yml

ENV PATH=/opt/conda/bin:$PATH \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1

LABEL org.opencontainers.image.source="https://github.com/bhagesh-h/methylarray" \
      org.opencontainers.image.licenses="MIT"
