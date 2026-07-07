FROM mambaorg/micromamba:ubuntu
LABEL maintainer="Danilo Tatoni <danilotatoni@gmail.com>"
# netmhcpan dependencies
USER root
RUN apt-get update -qq && \
    apt-get install -qq -y --no-install-recommends \
        bzip2 gcc g++ make zlib1g-dev wget tcsh gawk vim \
    && rm -rf /var/lib/apt/lists/*
RUN mkdir -p /opt/iedb && chown mambauser:mambauser /opt/iedb
USER mambauser
WORKDIR /opt/iedb
# IEDB MHC-I tools (netmhcpan). Pulled from a GitHub release mirror because
# downloads.iedb.org blocks GitHub Actions runner IPs. The mirror asset is
# version-less, so updating IEDB is a re-run of setup/mirror_iedb.sh with no
# change here. Override with --build-arg IEDB_MHCI_URL=... to pull elsewhere.
ARG IEDB_MHCI_URL=https://github.com/ctglab/ENEO/releases/download/iedb-tools/IEDB_MHC_I.tar.gz
RUN IEDB_MHCI_ARCHIVE="$(basename "$IEDB_MHCI_URL")" && \
    echo "Downloading IEDB MHC-I tools from $IEDB_MHCI_URL ..." && \
    if ! wget --tries=5 --waitretry=10 --timeout=30 "$IEDB_MHCI_URL" -O "$IEDB_MHCI_ARCHIVE"; then \
        echo "ERROR: failed to download IEDB MHC-I tools from $IEDB_MHCI_URL" >&2; \
        exit 1; \
    fi && \
    tar -xzf "$IEDB_MHCI_ARCHIVE" && \
    rm "$IEDB_MHCI_ARCHIVE"
# conda dependencies
RUN micromamba install -n base -y \
    -c bioconda -c conda-forge \
    python=3.10 \
    bedtools bcftools fastp tabix samtools pip scipy pandas bionumpy cyvcf2 numpy toml pyyaml \
    && micromamba clean --all --yes
# bind netmhcpan
ENV PATH="$MAMBA_ROOT_PREFIX/bin:$PATH:/opt/iedb/mhc_i/method/netmhcpan-4.1-executable/netmhcpan_4_1_executable/"
# Explicitly ensure the ARG is set for any subsequent RUN commands in this build stage
ARG MAMBA_DOCKERFILE_ACTIVATE=1
WORKDIR /opt
