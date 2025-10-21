# syntax=docker/dockerfile:1
FROM debian:bookworm-slim

# Avoid interactive prompts and reduce image size
ENV DEBIAN_FRONTEND=noninteractive

# Install samtools and cleanup
RUN apt-get update && \
    apt-get install -y --no-install-recommends samtools ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# Default workdir
WORKDIR /data

# Default entrypoint to samtools
ENTRYPOINT ["samtools"]
CMD ["--help"]