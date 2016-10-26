#!/bin/bash

GATK_JAR='/data/common/tools/gatk_v3.4.46/GenomeAnalysisTK.jar'
REF_FASTA='/data/common/refgenomes/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta'
DESIGN_DIR=$( dirname $0 )

designs=( 
    '32SNPtaqman_design.vcf'
    '59SNPtaqman_design.vcf'
)

bamPath=$1

if [[ ! -f ${bamPath} ]]; then echo "Provide path to bam"; exit 1; fi
if [[ ! "${bamPath}" =~ \.bam$ ]]; then echo "Only BAM files as input"; exit 1; fi

echo "[INFO] SNPcheck Genotying Starting"

fileName=$( basename $bamPath )
sampleName=$( echo $fileName | sed 's/\.bam//' )

for designFile in "${designs[@]}"; 
do 
    designName=$( echo $designFile | sed 's/_design\.vcf//' )
    designPath=${DESIGN_DIR}'/'${designFile}
    outFile="${sampleName}_${designName}.vcf"
    echo "[INFO] Genotying by ${designName} for ${sampleName} to ${outFile}"
    java -Xms2g -Xmx9g -jar ${GATK_JAR} -T UnifiedGenotyper -R ${REF_FASTA} -L ${designPath} --output_mode EMIT_ALL_SITES -I ${bamPath} -o ${outFile}
done

echo "[INFO] SNPcheck Genotying Done"
