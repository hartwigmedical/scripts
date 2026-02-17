#!/bin/bash

# Usage check, only permit 1 or 3 arguments
if [ "$#" -lt 1 ] || [ "$#" -gt 3 ] || [ "$#" -eq 2 ]; then
  echo "Usage: <setname> <input-bucket-name> <output-bucket-name>"
  echo "Ex: 123456_HMFregCORE_FS12345678_CORE0100000 research-pipeline-output-prod-1 example-output-bucket"
  echo "Input and output bucket-names are optional, defaults to diagnostic-pipeline-output-prod-1 and wgs-combined-snps-vcfs (used in production process)"
  echo "NOTE: When specifying a bucket, both input and output buckets must be provided."
  exit 1
else
  # Redirect output to log only after usage validated
  INPUT_DIR="$HOME/inputfiles"
  mkdir INPUT_DIR
  exec > "$INPUT_DIR/prepare_inputs.log" 2>&1
  set -euo pipefail
fi

setname=$1
DEFAULT_BUCKET="diagnostic-pipeline-output-prod-1"
BUCKET_NAME=${2:-$DEFAULT_BUCKET}

DEFAULT_OUTPUT_BUCKET="wgs-combined-snps-vcfs"
OUTPUT_BUCKET_NAME=${3:-$DEFAULT_OUTPUT_BUCKET}

MOUNT_POINT_BAM="$HOME/testdir/"
mkdir -p "$INPUT_DIR"

# Get barcodes
ISOLATION_BARCODE=$(echo "$setname" | cut -d'_' -f4)
SAMPLE_BARCODE=$(lama_get_patient_reporter_data "${ISOLATION_BARCODE}" | jq .tumorSampleBarcode | tr -d '"')
HOSPITAL_SAMPLE_LABEL=$(lama_get_patient_reporter_data "${ISOLATION_BARCODE}" | jq -r .hospitalSampleLabel)
REPORTING_ID=$(lama_get_patient_reporter_data "${ISOLATION_BARCODE}" | jq -r .reportingId)

if [[ -n "$HOSPITAL_SAMPLE_LABEL" && "$HOSPITAL_SAMPLE_LABEL" != "null" ]]; then
    CONVERTED_REPORTING_ID="${REPORTING_ID}-${HOSPITAL_SAMPLE_LABEL}"
else
    CONVERTED_REPORTING_ID="${REPORTING_ID}"
fi

# Mount buckets
mkdir -p "${MOUNT_POINT_BAM}"

fusermount -u "${MOUNT_POINT_BAM}" 2>/dev/null || true

gcsfuse --implicit-dirs "${BUCKET_NAME}" "${MOUNT_POINT_BAM}"

# SNP intervals
cp "/data/resources/reporting-resources/snps/id_snps_intervals.hg37.bed" "$INPUT_DIR/id_snps_intervals.hg37.bed"

# Copy hartwig SNP VCF
ALL_SNP_VCFS=($(find "${MOUNT_POINT_BAM}${setname}" -type f -path "*/snp_genotype/*output.vcf"))

if [ ${#ALL_SNP_VCFS[@]} -eq 0 ]; then
  echo "Warning: Geen SNP VCF bestanden gevonden in ${MOUNT_POINT_BAM}${setname}" >&2
fi

for SNP_VCF_PATH in "${ALL_SNP_VCFS[@]}"; do
    if [[ "$SNP_VCF_PATH" != *-ref* ]]; then
        cp "$SNP_VCF_PATH" "$INPUT_DIR/hartwig_snpfile_tum.vcf"
    fi
done



# Purple VCFs
cp "${MOUNT_POINT_BAM}${setname}/purple/"*.purple.germline.vcf.gz "$INPUT_DIR/"
cp "${MOUNT_POINT_BAM}${setname}/purple/"*.purple.somatic.vcf.gz "$INPUT_DIR/"

# Save metadata
echo "$SAMPLE_BARCODE" > "$INPUT_DIR/sample_barcode.txt"
echo "$CONVERTED_REPORTING_ID" > "$INPUT_DIR/converted_reporting_id.txt"
echo "$OUTPUT_BUCKET_NAME" > "$INPUT_DIR/output_bucket.txt"
cp "/data/tools/gatk/3.8.0/GenomeAnalysisTK.jar" "$INPUT_DIR/GenomeAnalysisTK.jar"

gsutil cp -r "${INPUT_DIR}" "gs://${OUTPUT_BUCKET_NAME}/${setname}/"
gsutil cp  "$INPUT_DIR/prepare_inputs.log" "gs://${OUTPUT_BUCKET_NAME}/${setname}/"

# Unmount & cleanup
fusermount -u "${MOUNT_POINT_BAM}"
rm -r "${MOUNT_POINT_BAM}"
rm -r "${INPUT_DIR}"

echo "Done. Files prepared in: gs://${OUTPUT_BUCKET_NAME}/${setname}/"

