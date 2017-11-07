#!/bin/bash

MNV_DETECTOR=/data/common/tools/mnvdetector_v1.0/mnv-detector.jar
MNV_VALIDATOR=/data/common/tools/mnvdetector_v1.0/mnv-validator.jar
BAM_SLICER_SCRIPT=/data/common/repos/scripts/hmftools/bamslicer/bam-slicer.sh

SNPEFF_VRSN="v4.3s"
SNPEFF_FLAG=" -hgvs -lof -no-downstream -no-upstream -no-intergenic -noShiftHgvs"
SNPEFF_ROOT=/data/common/tools/snpEff_${SNPEFF_VRSN}/
SNPEFF_DB="GRCh37.75"


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

if ! java -jar $MNV_DETECTOR \
    -v $VCF \
    -bed_out $MNV_BED \
    -vcf_out $MNV_VCF ;
then exit 1
fi

rm "${MNV_VCF}.idx"

if ! $BAM_SLICER_SCRIPT $SAMPLE $MNV_BED ;
then exit 1
fi

if ! java -jar $MNV_VALIDATOR \
    -v $VCF \
    -b $SLICED_BAM \
    -o $FINAL_VCF ;
then exit 1
fi

rm "${FINAL_VCF}.idx"

if ! java -jar ${SNPEFF_ROOT}/snpEff.jar \
    -c "${SNPEFF_ROOT}/snpEff.config" "${SNPEFF_DB}" \
    -v $FINAL_VCF \
    ${SNPEFF_FLAG} \
    > $ANNOTATED_VCF ;
then exit 1
fi

rm $FINAL_VCF
rm snpEff_genes.txt
rm snpEff_summary.html
