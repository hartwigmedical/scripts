################## BASE IMAGE ######################

FROM ubuntu:18.04

################## METADATA ######################
LABEL base.image="biocontainers:latest"
LABEL version="1"
LABEL software="chromaticHMF-GRIDSS2"
LABEL software.version="0.0.1"
LABEL about.summary="ChromaticHMF-GRIDSS2 somatic genomic rearrangement toolkit"
LABEL about.home="https://github.com/PapenfussLab/gridss"
LABEL about.tags="Genomics"

#RUN apt-get clean all && apt-get update
# samtools and bwa  dependencies
RUN apt-get install -y zlib1g libncurses5 libbz2 liblzma libcurl4
# samtools and bwa compilation dependencies
RUN apt-get install -y build-essential zlib1g-dev libncurses5-dev libbz2-dev liblzma-dev libcurl4-openssl-dev
# circos dependencies
# See https://hub.docker.com/r/dennishazelett/circos/dockerfile for ubuntu 14 version
sudo apt-get -y install libgd-dev
# java runtime required for GRIDSS and Chromatic HMF
RUN apt-get install -y openjdk-8-jre-headless r-base
# install R packages