#!/usr/bin/env bash

# ./extractRefGenomeChromosomes.sh /data/common/refgenomes/Homo_sapiens.GRCh38.no.alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna /data/common/refgenomes/Homo_sapiens.GRCh38.no.alt/chr/

refGenome=$1 && shift
outputDir=$1 && shift
samtools=/data/common/tools/samtools_v1.9/samtools

${samtools} faidx ${refGenome} chr1 > ${outputDir}/chr1.fasta 
${samtools} faidx ${refGenome} chr2> ${outputDir}/chr2.fasta
${samtools} faidx ${refGenome} chr3 > ${outputDir}/chr3.fasta
${samtools} faidx ${refGenome} chr4 > ${outputDir}/chr4.fasta
${samtools} faidx ${refGenome} chr5 > ${outputDir}/chr5.fasta
${samtools} faidx ${refGenome} chr6 > ${outputDir}/chr6.fasta
${samtools} faidx ${refGenome} chr7 > ${outputDir}/chr7.fasta
${samtools} faidx ${refGenome} chr8 > ${outputDir}/chr8.fasta
${samtools} faidx ${refGenome} chr9 > ${outputDir}/chr9.fasta
${samtools} faidx ${refGenome} chr10 > ${outputDir}/chr10.fasta
${samtools} faidx ${refGenome} chr11 > ${outputDir}/chr11.fasta
${samtools} faidx ${refGenome} chr12 > ${outputDir}/chr12.fasta
${samtools} faidx ${refGenome} chr13 > ${outputDir}/chr13.fasta
${samtools} faidx ${refGenome} chr14 > ${outputDir}/chr14.fasta
${samtools} faidx ${refGenome} chr15 > ${outputDir}/chr15.fasta
${samtools} faidx ${refGenome} chr16 > ${outputDir}/chr16.fasta
${samtools} faidx ${refGenome} chr17 > ${outputDir}/chr17.fasta
${samtools} faidx ${refGenome} chr18 > ${outputDir}/chr18.fasta
${samtools} faidx ${refGenome} chr19 > ${outputDir}/chr19.fasta
${samtools} faidx ${refGenome} chr20 > ${outputDir}/chr20.fasta
${samtools} faidx ${refGenome} chr21 > ${outputDir}/chr21.fasta
${samtools} faidx ${refGenome} chr22 > ${outputDir}/chr22.fasta
${samtools} faidx ${refGenome} chrY > ${outputDir}/chrY.fasta
${samtools} faidx ${refGenome} chrX > ${outputDir}/chrX.fasta 
