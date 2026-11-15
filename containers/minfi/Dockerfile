#
# Image for the minfi and ENmix module (modules/local/minfi/preprocess).
#
# The software comes only from the modules' environment.yml, so the Docker,
# Singularity/Apptainer and conda profiles run the same package versions.
# Build from the repository root:
#   docker build -f containers/minfi/Dockerfile -t bhagesh/methylarray-minfi:1.56.0 .
#
FROM mambaorg/micromamba:2.10.0-ubuntu24.04

COPY --chown=$MAMBA_USER:$MAMBA_USER modules/local/minfi/preprocess/environment.yml /tmp/environment.yml
RUN micromamba install --yes --name base --file /tmp/environment.yml \
 && micromamba clean --all --yes \
 && rm /tmp/environment.yml

ENV PATH=/opt/conda/bin:$PATH \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1

# Fail the build if any package cannot be loaded
RUN Rscript -e 'for (p in c("minfi", "ENmix", "IlluminaHumanMethylationEPICv2manifest", "IlluminaHumanMethylationEPICv2anno.20a1.hg38", "HDF5Array", "SummarizedExperiment")) suppressPackageStartupMessages(library(p, character.only = TRUE)); cat("R", as.character(getRversion()), "minfi", as.character(packageVersion("minfi")), "ENmix", as.character(packageVersion("ENmix")), "\n")'

LABEL org.opencontainers.image.source="https://github.com/bhagesh-h/methylarray" \
      org.opencontainers.image.description="minfi 1.56.0 and ENmix 1.46.2 with the EPIC v2 manifest and annotation, HDF5Array and SummarizedExperiment (R 4.5.3)" \
      org.opencontainers.image.licenses="MIT"
