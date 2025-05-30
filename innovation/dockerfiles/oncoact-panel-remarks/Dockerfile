# Use Debian Buster slim as base image
FROM debian:buster-slim


# Install OS-level dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-dev \
    build-essential \
    libhts-dev \
    libbz2-dev \
    liblzma-dev \
    wget \
    unzip \
    ca-certificates \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip to the latest version
RUN python3 -m pip install --upgrade pip

# Install required Python packages
RUN pip install --no-cache-dir \
    pandas \
    cyvcf2 \
    openpyxl

SHELL ["bash", "-c"]

# Copy workflow scripts into the container.
RUN mkdir -p /scripts
COPY oncoact_panel_remarks.py /scripts
# Find a clever way to retrieve these automatically, now update manually!!!
COPY excludedExons.tsv /data/resources/ops/panel/

ENTRYPOINT ["/scripts/oncoact_panel_remarks.py"]