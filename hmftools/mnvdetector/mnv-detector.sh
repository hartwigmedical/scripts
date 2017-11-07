#!/bin/bash

MNV_DETECTOR=/data/common/tools/mnvdetector_v1.0/mnv-detector.jar
MNV_VALIDATOR=/data/common/tools/mnvdetector_v1.0/mnv-validator.jar
BAM_SLICER_SCRIPT=/data/common/repos/scripts/hmftools/bamslicer/bam-slicer.sh

SNPEFF_VRSN="v4.3s"
SNPEFF_FLAG=" -hgvs -lof -no-downstream -no-upstream -no-intergenic -noShiftHgvs"
SNPEFF_ROOT=/data/common/tools/snpEff_${SNPEFF_VRSN}/
SNPEFF_DB="GRCh37.75"

RUN=$1
#VCF=$1
#SAMPLE=$2

#if [ -z "$VCF" ] || [ -z "$SAMPLE" ];
if [ -z "$RUN" ];
  then
    echo "Usage: $ ./mnv_detector.sh <runpath>"
    echo "   runpath = Path to pipeline run dir"
    echo "Notes:"
    echo "  - should have a ./metadata file in rundir"
    echo "  - somatic vcf should meet format *_post_processed.vcf"
    exit 1
fi

VCF=$( find $RUN -wholename "*/somaticVariants/*_post_processed.vcf" )
SAMPLE=$( cat $RUN/metadata | jq -r '.tumor_sample' )

#echo "SAMPLE: $SAMPLE"

## some sanity checks
if [[ ! -d "${RUN}" ]]; then echo "[EXIT] Run ($RUN) not found" && exit 1; fi
if [[ ! -f "${VCF}" ]]; then echo "[EXIT] Vcf not found" && exit 1; fi
if [[ -z "${SAMPLE}" ]]; then echo "[EXIT] Sample not found" && exit 1; fi

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
