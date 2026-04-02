FROM debian:13-slim

ARG QSEE_VERSION=1.0-beta.9

# Install Java 21, R, and build dependencies for R packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    openjdk-21-jre-headless \
    wget \
    ca-certificates \
    r-base \
    r-base-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install R packages required by qsee plot_qc.R
RUN Rscript -e 'install.packages(c("dplyr", "ggplot2", "patchwork", "scales"), repos="https://cloud.r-project.org")' \
    && Rscript -e 'library(dplyr); library(ggplot2); library(patchwork); library(scales)'

# Install qsee jar
RUN mkdir -p /opt/tools && \
    wget -q "https://github.com/hartwigmedical/hmftools/releases/download/qsee-v${QSEE_VERSION}/qsee_v${QSEE_VERSION}.jar" \
    -O /opt/tools/qsee.jar

# Verify installation
RUN java -version && \
    java -cp /opt/tools/qsee.jar \
    com.hartwig.hmftools.qsee.cohort.CohortPercentilesTrainer --version 2>&1 || true

CMD ["/bin/bash"]
