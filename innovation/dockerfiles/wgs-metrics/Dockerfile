FROM debian:13-slim

# Install Java 21
RUN apt-get update && apt-get install -y --no-install-recommends \
    openjdk-21-jre-headless \
    && rm -rf /var/lib/apt/lists/*

ARG BAMTOOLS_VERSION=1.5
ARG GATK_VERSION=4.6.1.0
ARG SAMTOOLS_VERSION=1.21

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    python3 \
    r-base \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libncurses5-dev \
    autoconf \
    unzip \
    ca-certificates \
    && ln -s /usr/bin/python3 /usr/bin/python \
    && rm -rf /var/lib/apt/lists/*

# Install bam-tools
RUN mkdir -p /opt/tools && \
    wget -q "https://github.com/hartwigmedical/hmftools/releases/download/bam-tools-v${BAMTOOLS_VERSION}/bam-tools_v${BAMTOOLS_VERSION}.jar" \
    -O /opt/tools/bam-tools_v${BAMTOOLS_VERSION}.jar

# Install GATK (for CollectGcBiasMetrics)
RUN wget -q "https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip" \
    -O /tmp/gatk.zip && \
    unzip -q /tmp/gatk.zip -d /opt/ && \
    ln -s /opt/gatk-${GATK_VERSION}/gatk /usr/local/bin/gatk && \
    rm /tmp/gatk.zip

# Install samtools
RUN wget -q "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" \
    -O /tmp/samtools.tar.bz2 && \
    cd /tmp && \
    tar xjf samtools.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    rm -rf /tmp/samtools*

# Verify installations
RUN java -version && \
    samtools --version | head -1 && \
    gatk --version | head -3 && \
    java -cp /opt/tools/bam-tools_v${BAMTOOLS_VERSION}.jar \
    com.hartwig.hmftools.bamtools.metrics.BamMetrics --version 2>&1 || true

CMD ["/bin/bash"]