#!/bin/bash
exec > "$HOME/script.log" 2>&1
set -euo pipefail

# Usage check
if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
  echo "Usage: <setname> <bucket-name>"
  echo "Ex: 123456_HMFregCORE_FS12345678_CORE0100000 research-pipeline-output-prod-1"
  echo "Bucket-name is optional, defaults to diagnostic-pipeline-output-prod-1"
  exit 1
fi

setname=$1
DEFAULT_BUCKET="diagnostic-pipeline-output-prod-1"
BUCKET_NAME=${2:-$DEFAULT_BUCKET}
MOUNT_POINT_BAM="$HOME/testdir/"
MOUNT_POINT_REFGENOME="$HOME/refgenometmp/"
ISOLATION_BARCODE=$(echo "$setname" | cut -d'_' -f4)
SAMPLE_BARCODE=$(lama_get_patient_reporter_data ${ISOLATION_BARCODE} | jq .tumorSampleBarcode | tr -d '"')

# Paths relative to bucket root (adjust as needed)
REFERENCE="${MOUNT_POINT_REFGENOME}/reference_genome/37/Homo_sapiens.GRCh37.GATK.illumina.fasta"
INTERVALS="/data/resources/reporting-resources/snps/id_snps_intervals.hg37.bed"

# Temp file names
SNP_OUTPUT_VCF="${SAMPLE_BARCODE}_snp_genotype_output.vcf"
MERGED_OUTPUT_VCF="${SAMPLE_BARCODE}_merged.vcf"
FINAL_OUTPUT_VCF="${SAMPLE_BARCODE}_merged_final.vcf"

echo ${REFERENCE}
echo ${INTERVALS}
echo ${SNP_OUTPUT_VCF}
echo ${MERGED_OUTPUT_VCF}

GATK="/data/tools/gatk/3.8.0/GenomeAnalysisTK.jar"  # Adjust if needed

echo "Mounting gs://${BUCKET_NAME} at ${MOUNT_POINT_BAM}..."
mkdir -p "${MOUNT_POINT_BAM}"
if mountpoint -q "${MOUNT_POINT_BAM}"; then
  echo "Mount point ${MOUNT_POINT_BAM} already mounted. Unmounting first..."
  fusermount -u "${MOUNT_POINT_BAM}"
fi


mkdir -p "${MOUNT_POINT_REFGENOME}"
if mountpoint -q "${MOUNT_POINT_REFGENOME}"; then
  echo "Unmounting first"
  fusermount -u "${MOUNT_POINT_REFGENOME}"
fi

/usr/bin/gcsfuse --implicit-dirs "${BUCKET_NAME}" "${MOUNT_POINT_BAM}"
/usr/bin/gcsfuse --implicit-dirs "common-resources" "${MOUNT_POINT_REFGENOME}"

# Find all .bam files under 'aligner' subfolders
ALL_BAMS=($(find "${MOUNT_POINT_BAM}${setname}" -type f -path "*/aligner/*.bam"))

# Separate REF and TUM
for BAM_PATH in "${ALL_BAMS[@]}"; do
    if [[ "$BAM_PATH" == *-ref.bam ]]; then
        BAM_REF="$BAM_PATH"
    else
        BAM_TUM="$BAM_PATH"
    fi
done

echo "Tum BAM: ${BAM_TUM}"
echo "Running GATK HaplotypeCaller..."
/usr/lib/jvm/adoptopenjdk-8-hotspot-amd64/jre/bin/java -Xmx20G \
  -jar "$GATK" \
  -T UnifiedGenotyper \
  -nct $(grep -c '^processor' /proc/cpuinfo) \
  --input_file $BAM_TUM \
  -o $HOME/${SNP_OUTPUT_VCF} \
  -L ${INTERVALS} \
  --reference_sequence $REFERENCE \
  --output_mode EMIT_ALL_SITES \

somatic_vcf=$(ls ${MOUNT_POINT_BAM}${setname}/purple/*.purple.somatic.vcf.gz)
germline_vcf=$(ls ${MOUNT_POINT_BAM}${setname}/purple/*.purple.germline.vcf.gz)
zcat ${MOUNT_POINT_BAM}${setname}/purple/*.purple.somatic.vcf.gz | grep -E '^#|REPORTED' > $HOME/${SAMPLE_BARCODE}_somatic_filtered.vcf
zcat ${MOUNT_POINT_BAM}${setname}/purple/*.purple.germline.vcf.gz | grep -E '^#|REPORTED' > $HOME/${SAMPLE_BARCODE}_germline_filtered.vcf

/usr/lib/jvm/adoptopenjdk-8-hotspot-amd64/jre/bin/java -jar "$GATK" \
  -T CombineVariants \
  -R $REFERENCE \
  --variant $HOME/${SAMPLE_BARCODE}_somatic_filtered.vcf \
  --variant $HOME/${SAMPLE_BARCODE}_germline_filtered.vcf \
  --variant $HOME/$SNP_OUTPUT_VCF  \
  -o $HOME/${MERGED_OUTPUT_VCF}\
  -genotypeMergeOptions UNSORTED


echo "Unmounting ${MOUNT_POINT_REFGENOME}..."
fusermount -u "${MOUNT_POINT_REFGENOME}"

echo "Unmounting ${MOUNT_POINT_BAM}..."
fusermount -u "${MOUNT_POINT_BAM}"

awk '
BEGIN { OFS="\t" }

# Meta-header lines: print unchanged
/^##/ { print; next }

# Header line
/^#CHROM/ {
    total_cols = NF
    second_last_col = NF - 1
    last_col = NF
    # Print all columns except the last
    for (i = 1; i < NF; i++) {
        printf "%s%s", $i, (i < NF - 1 ? OFS : ORS)
    }
    next
}

# Data lines
{
    if ($(NF-2) == "./.") {
        $(NF-2) = $NF
    }
    # Print all columns except the last
    for (i = 1; i < NF; i++) {
        printf "%s%s", $i, (i < NF - 1 ? OFS : ORS)
    }
}
' $HOME/${MERGED_OUTPUT_VCF} > $HOME/${FINAL_OUTPUT_VCF}

gsutil cp $HOME/${FINAL_OUTPUT_VCF} gs://wgs-combined-snps-vcfs/${setname}/${FINAL_OUTPUT_VCF}
gsutil cp $HOME/script.log gs://wgs-combined-snps-vcfs/${setname}/${SAMPLE_BARCODE}_merge_script.log

rm -r "${MOUNT_POINT_BAM}"
rm -r "${MOUNT_POINT_REFGENOME}"
rm $HOME/${SAMPLE_BARCODE}_somatic_filtered.vcf
rm $HOME/${SAMPLE_BARCODE}_somatic_filtered.vcf.idx
rm $HOME/${SAMPLE_BARCODE}_germline_filtered.vcf
rm $HOME/${SAMPLE_BARCODE}_germline_filtered.vcf.idx
rm $HOME/${MERGED_OUTPUT_VCF}
rm $HOME/${FINAL_OUTPUT_VCF}
rm $HOME/${SNP_OUTPUT_VCF}
rm $HOME/script.log
