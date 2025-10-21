# syntax=docker/dockerfile:1
FROM debian:bookworm-slim

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies: Java runtime + unzip + certificates
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        openjdk-17-jre-headless \
        unzip \
        wget \
        ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# Set FastQC version
ENV FASTQC_VERSION=0.12.1
ENV FASTQC_DIR=/opt/FastQC

# Download and install FastQC
RUN wget -q https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip -O /tmp/fastqc.zip && \
    unzip /tmp/fastqc.zip -d /opt && \
    rm /tmp/fastqc.zip && \
    chmod +x ${FASTQC_DIR}/fastqc && \
    ln -s ${FASTQC_DIR}/fastqc /usr/local/bin/fastqc

# Default working directory
WORKDIR /data

# Default entrypoint
ENTRYPOINT ["fastqc"]
CMD ["--help"]