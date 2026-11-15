#
# Image for the SeSAMe modules (modules/local/sesame/*).
#
# The software comes only from the modules' environment.yml, so the Docker,
# Singularity/Apptainer and conda profiles run the same package versions.
# Build from the repository root:
#   docker build -f containers/sesame/Dockerfile -t bhagesh/methylarray-sesame:1.28.1 .
#
FROM mambaorg/micromamba:2.10.0-ubuntu24.04

COPY --chown=$MAMBA_USER:$MAMBA_USER modules/local/sesame/preprocess/environment.yml /tmp/environment.yml
RUN micromamba install --yes --name base --file /tmp/environment.yml \
 && micromamba clean --all --yes \
 && rm /tmp/environment.yml

ENV PATH=/opt/conda/bin:$PATH \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1

# Fail the build if any package cannot be loaded
RUN Rscript -e 'for (p in c("sesame", "sesameData", "ExperimentHub", "HDF5Array", "SummarizedExperiment")) suppressPackageStartupMessages(library(p, character.only = TRUE)); cat("R", as.character(getRversion()), "sesame", as.character(packageVersion("sesame")), "\n")'

LABEL org.opencontainers.image.source="https://github.com/bhagesh-h/methylarray" \
      org.opencontainers.image.description="SeSAMe 1.28.1 with sesameData, ExperimentHub, HDF5Array and SummarizedExperiment (R 4.5.3)" \
      org.opencontainers.image.licenses="MIT"
