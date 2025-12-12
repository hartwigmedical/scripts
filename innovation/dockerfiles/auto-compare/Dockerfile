# Use Debian Buster slim as base image
FROM debian:trixie-slim


# Install OS-level dependencies, including jq for JSON parsing
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3-full \
    python3-dev \
    build-essential \
    libhts-dev \
    libbz2-dev \
    liblzma-dev \
    openjdk-25-jre-headless \
    wget \
    unzip \
    ca-certificates \
    jq \
    && rm -rf /var/lib/apt/lists/*

# Create a venv for all Python libraries (avoids PEP 668 issues)
RUN python3 -m venv /opt/venv

# Add venv to PATH to make this the default python
ENV PATH="/opt/venv/bin:${PATH}"

# Upgrade pip to the latest version
RUN python3 -m pip install --upgrade pip

# Install required Python packages
RUN pip install --no-cache-dir \
    pandas \
    cyvcf2 \
    openpyxl

SHELL ["bash", "-c"]

# Get COMPAR jar
RUN wget -O compar.jar https://github.com/hartwigmedical/hmftools/releases/download/compar-v1.4.0/compar_v1.4.0.jar

# Copy workflow scripts into the container.
RUN mkdir -p /scripts
COPY create_excel_compar.py /scripts
COPY DriverGenePanel.37.tsv /data/resources/gene_panel/
COPY DriverGenePanel.38.tsv /data/resources/gene_panel/