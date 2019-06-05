################## BASE IMAGE ######################

FROM ubuntu:18.04

################## METADATA ######################
LABEL base.image="biocontainers:latest"
LABEL version="1"
LABEL software="GRIDSS PURPLE LINX"
LABEL software.version="0.0.1"
LABEL about.summary="GRIDSS PURPLE LINX somatic genomic rearrangement toolkit"
LABEL about.home="https://github.com/PapenfussLab/gridss"
LABEL about.tags="Genomics"

#RUN apt-get clean all && apt-get update
RUN apt-get install -y openjdk-8-jre-headless r-base samtools sambamba bwa circos

RUN mkdir /jar /scripts /data
COPY ../../../dev/gridss/target/gridss-2.2.3-gridss-jar-with-dependencies.jar /jar
COPY ../../hmftools/amber/target/amber-2.3-jar-with-dependencies.jar /jar
COPY ../../hmftools/bachelor/target/bachelor-1.5-jar-with-dependencies.jar /jar
COPY ../../hmftools/purity-ploidy-estimator/target/purity-ploidy-estimator-2.29-jar-with-dependencies.jar /jar
COPY ../../hmftools/sage/target/sage-1.1-jar-with-dependencies.jar /jar
COPY ../../hmftools/sv-analyser/target/sv-analyser-1.0-jar-with-dependencies.jar /jar
COPY ../../hmftools/sv-graphs/target/sv-graphs-1.0-jar-with-dependencies.jar /jar
COPY ../../hmftools/sv-visualiser/target/sv-visualiser-1.0-jar-with-dependencies.jar /jar
COPY ../../hmftools/variant-annotator/target/variant-annotator-2.15-jar-with-dependencies.jar /jar
COPY standalone-gridss-purple-linx.sh /scripts


