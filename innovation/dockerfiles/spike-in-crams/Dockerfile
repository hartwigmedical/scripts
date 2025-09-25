FROM debian:bullseye-slim

ENV DEBIAN_FRONTEND=noninteractive


RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        samtools \
        bc \
        python3 \
        python3-pip \
        vcftools \
        ca-certificates \
    && apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir \
    pysam==0.22.1

WORKDIR /work

CMD ["bash"]