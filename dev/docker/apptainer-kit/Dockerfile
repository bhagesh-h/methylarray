#
# methylarray developer toolbox with Apptainer.
#
# Extends the devkit image so that the pipeline's singularity and apptainer
# profiles can be tested on a machine where only Docker is installed.
# The official Apptainer packages are used, including the setuid starter,
# because the host restricts unprivileged user namespaces.

ARG DEVKIT_IMAGE=methylarray-devkit:dev
FROM ${DEVKIT_IMAGE}

ARG APPTAINER_VER=1.5.5
ARG APPTAINER_SHA256=25095d3711504b3e8657dbb2e172ef8db9db43c4503d940c9c36b105db8dce77
ARG APPTAINER_SUID_SHA256=d9dd6b26540f200f394465bf870dcf7c3a392effa57d5786f8682cde2dafa493

LABEL org.opencontainers.image.title="methylarray-apptainer-kit" \
      org.opencontainers.image.description="Developer toolbox for methylarray with Apptainer" \
      org.opencontainers.image.source="https://github.com/bhagesh-h/methylarray" \
      org.opencontainers.image.licenses="MIT"

USER root

RUN base="https://github.com/apptainer/apptainer/releases/download/v${APPTAINER_VER}" \
 && python -c "import urllib.request as u, sys; [u.urlretrieve(f'${base}/{n}', f'/tmp/{n}') for n in sys.argv[1:]]" \
        "apptainer_${APPTAINER_VER}_amd64.deb" "apptainer-suid_${APPTAINER_VER}_amd64.deb" \
 && printf '%s  %s\n' "${APPTAINER_SHA256}" "/tmp/apptainer_${APPTAINER_VER}_amd64.deb" \
                      "${APPTAINER_SUID_SHA256}" "/tmp/apptainer-suid_${APPTAINER_VER}_amd64.deb" \
    | sha256sum --check --strict \
 && apt-get update \
 && DEBIAN_FRONTEND=noninteractive apt-get install --yes --no-install-recommends \
        "/tmp/apptainer_${APPTAINER_VER}_amd64.deb" "/tmp/apptainer-suid_${APPTAINER_VER}_amd64.deb" \
        squashfs-tools tzdata \
 && rm -rf /var/lib/apt/lists/* "/tmp/apptainer_${APPTAINER_VER}_amd64.deb" "/tmp/apptainer-suid_${APPTAINER_VER}_amd64.deb" \
 && ln -sf /usr/share/zoneinfo/Etc/UTC /etc/localtime \
 && apptainer --version \
 && singularity --version

USER $MAMBA_USER
