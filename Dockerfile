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
RUN IEDB_MHCI_URL="https://downloads.iedb.org/tools/mhci/LATEST" && \
    WGET_OPTS="--tries=5 --waitretry=10 --timeout=30 --user-agent=Mozilla/5.0" && \
    echo "Fetching IEDB MHC-I directory listing from $IEDB_MHCI_URL/ ..." && \
    if ! IEDB_MHCI_LISTING=$(wget $WGET_OPTS -qO- "$IEDB_MHCI_URL/"); then \
        echo "ERROR: failed to fetch IEDB listing from $IEDB_MHCI_URL/" >&2; \
        exit 1; \
    fi && \
    IEDB_MHCI_ARCHIVE=$(printf '%s\n' "$IEDB_MHCI_LISTING" | grep -oE 'IEDB_MHC_I-[0-9.]+\.tar\.gz' | sort -V | tail -n1 || true) && \
    if [ -z "$IEDB_MHCI_ARCHIVE" ]; then \
        echo "ERROR: could not resolve the latest IEDB MHC-I archive name; listing returned was:" >&2; \
        printf '%s\n' "$IEDB_MHCI_LISTING" >&2; \
        exit 1; \
    fi && \
    echo "Resolved IEDB MHC-I archive: $IEDB_MHCI_ARCHIVE" && \
    wget $WGET_OPTS "$IEDB_MHCI_URL/$IEDB_MHCI_ARCHIVE" && \
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
