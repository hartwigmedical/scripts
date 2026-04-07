# syntax=docker/dockerfile:1
FROM debian:bookworm-slim

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        samtools \
        bcftools \
        bedtools \
        gawk \
        ca-certificates && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /data
