# Use Debian Buster slim as base image
FROM debian:buster-slim


# Install OS-level dependencies, including jq for JSON parsing
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-dev \
    curl \
    build-essential \
    libhts-dev \
    libbz2-dev \
    liblzma-dev \
    openjdk-11-jre-headless \
    wget \
    unzip \
    ca-certificates \
    jq \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip to the latest version
RUN python3 -m pip install --upgrade pip

# Install required Python packages
RUN pip install --no-cache-dir \
    pandas \
    cyvcf2 \
    openpyxl

SHELL ["bash", "-c"]

# Get COMPAR jar
RUN wget -O compar.jar https://github.com/hartwigmedical/hmftools/releases/download/compar-v1.3.4/compar_v1.3.4.jar

# Copy workflow scripts into the container.
RUN mkdir -p /scripts
COPY create_excel_compar.py /scripts
COPY DriverGenePanel.37.tsv /data/resources/gene_panel/
COPY DriverGenePanel.38.tsv /data/resources/gene_panel/