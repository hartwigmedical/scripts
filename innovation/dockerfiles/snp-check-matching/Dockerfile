# syntax=docker/dockerfile:1
FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

# OS deps (includes bash so RUN/CMD can use it)
RUN apt-get update && apt-get install -y --no-install-recommends \
      bash \
      curl \
      openjdk-21-jre-headless \
      ca-certificates \
      unzip \
      jq \
    && rm -rf /var/lib/apt/lists/*

SHELL ["bash","-o","pipefail","-c"]

# Python deps
RUN pip install --upgrade pip \
 && pip install pandas pysam mysql-connector-python

# AMBER jar
ARG AMBER_VERSION=4.2
RUN install -d /opt/amber /data/resources /data/scripts \
 && curl -fsSL -o /opt/amber/amber.jar \
    "https://github.com/hartwigmedical/hmftools/releases/download/amber-v${AMBER_VERSION}/amber_v${AMBER_VERSION}.jar"
ENV AMBER_JAR=/opt/amber/amber.jar

WORKDIR /data
COPY Amber.snpcheck.38.vcf /data/resources/
COPY match_amber_sample.py /data/scripts/

# Default to an interactive bash
CMD ["bash"]