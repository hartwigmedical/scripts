#!/bin/bash

MNV_DETECTOR=/data/common/tools/mnvdetector_v1.0/mnv-detector.jar
MNV_VALIDATOR=/data/common/tools/mnvdetector_v1.0/mnv-validator.jar
BAM_SLICER_SCRIPT=/data/common/repos/scripts/hmftools/bamslicer/bam-slicer.sh
SNP_EFF=/data/common/tools/snpEff_v4.1h/snpEff.jar

VCF=$1
SAMPLE=$2

if [ -z "$VCF" ] || [ -z "$SAMPLE" ];
  then
    echo "Usage: $ ./mnv_detector.sh vcf sample"
    echo "   vcf        vcf input file"
    echo "   sample     sample to search for. e.g CPCT11111111T"
    exit 1
fi

SLICED_BAM="${SAMPLE}_dedup.realigned.sliced.bam"
VCF_FILE_NAME=`basename $VCF`
VCF_NAME="${VCF_FILE_NAME%.*}"
MNV_BED="${VCF_NAME}.mnvs.bed"
MNV_VCF="${VCF_NAME}.potential_mnvs.vcf"
FINAL_VCF="${VCF_NAME}.mnvs.vcf"
ANNOTATED_VCF="${VCF_NAME}.mnvs.annotated.vcf"

java -jar $MNV_DETECTOR \
    -v $VCF \
    -bed_out $MNV_BED \
    -vcf_out $MNV_VCF

rm "${MNV_VCF}.idx"

$BAM_SLICER_SCRIPT $SAMPLE $MNV_BED

java -jar $MNV_VALIDATOR \
    -v $VCF \
    -b $SLICED_BAM \
    -o $FINAL_VCF

rm "${FINAL_VCF}.idx"

java -jar $SNP_EFF \
    -c /data/common/tools/snpEff_v4.1h/snpEff.config "GRCh37.74" \
    -v $FINAL_VCF \
    -hgvs -lof -no-downstream -no-upstream -no-intergenic \
    > $ANNOTATED_VCF

rm $FINAL_VCF
rm snpEff_genes.txt
rm snpEff_summary.html
