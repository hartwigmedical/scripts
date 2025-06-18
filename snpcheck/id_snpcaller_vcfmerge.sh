#!/bin/bash
exec > "$HOME/script.log" 2>&1
set -euo pipefail

# Usage check, only permit 1 or 3 arguments
if [ "$#" -lt 1 ] || [ "$#" -gt 3 ] || [ "$#" -eq 2 ]; then
  echo "Usage: <setname> <input-bucket-name> <output-bucket-name>"
  echo "Ex: 123456_HMFregCORE_FS12345678_CORE0100000 research-pipeline-output-prod-1 example-output-bucket"
  echo "Input and output bucket-names are optional, defaults to diagnostic-pipeline-output-prod-1 and wgs-combined-snps-vcfs (used in production process)"
  echo "NOTE: When specifying a bucket, both input and output buckets must be provided."
  exit 1
fi

setname=$1
DEFAULT_BUCKET="diagnostic-pipeline-output-prod-1"
BUCKET_NAME=${2:-$DEFAULT_BUCKET}

DEFAULT_OUTPUT_BUCKET="wgs-combined-snps-vcfs"
OUTPUT_BUCKET_NAME=${3:-$DEFAULT_OUTPUT_BUCKET}

# Mount dirs
MOUNT_POINT_BAM="$HOME/testdir/"
MOUNT_POINT_REFGENOME="$HOME/refgenometmp/"
echo "Mount point for BAM files: $MOUNT_POINT_BAM"
echo "Mount point for reference genome: $MOUNT_POINT_REFGENOME"

# Get barcodes
ISOLATION_BARCODE=$(echo "$setname" | cut -d'_' -f4)
SAMPLE_BARCODE=$(lama_get_patient_reporter_data ${ISOLATION_BARCODE} | jq .tumorSampleBarcode | tr -d '"')
echo "Extracted ISOLATION_BARCODE: $ISOLATION_BARCODE"
echo "Retrieved SAMPLE_BARCODE from LAMA: $SAMPLE_BARCODE"

# Paths relative to bucket root
REFERENCE="${MOUNT_POINT_REFGENOME}/reference_genome/37/Homo_sapiens.GRCh37.GATK.illumina.fasta"
INTERVALS="/data/resources/reporting-resources/snps/id_snps_intervals.hg37.bed"
echo "Reference genome path: $REFERENCE"
echo "SNP intervals path: $INTERVALS"

# Temp file names
SNP_OUTPUT_VCF="${SAMPLE_BARCODE}_snp_genotype_output.vcf"
MERGED_OUTPUT_VCF="${SAMPLE_BARCODE}_merged.vcf"
FINAL_OUTPUT_VCF="${SAMPLE_BARCODE}_reported_variants_and_snps_combined.vcf"

# GATK jar location
GATK="/data/tools/gatk/3.8.0/GenomeAnalysisTK.jar"

# Mount required buckets
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

/usr/bin/gcsfuse --implicit-dirs "${BUCKET_NAME}" "${MOUNT_POINT_BAM}" || {
  echo "Error: Failed to mount BAM bucket for tumor BAM"
  exit 1
}

/usr/bin/gcsfuse --implicit-dirs "common-resources" "${MOUNT_POINT_REFGENOME}"|| {
  echo "Error: Failed to mount common-resources bucket for reference genome"
  exit 1
}

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

if [[ -z "$BAM_TUM" ]]; then
  echo "Error: Could not identify tumor BAM file."
  exit 1
fi

echo "Tum BAM: ${BAM_TUM}"

# Calling id snps
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

# Get snpfile molecular pipeline
ALL_SNP_VCFS=($(find "${MOUNT_POINT_BAM}${setname}" -type f -path "*/snp_genotype/*output.vcf"))

for SNP_VCF_PATH in "${ALL_SNP_VCFS[@]}"; do
    if [[ "$SNP_VCF_PATH" == *-ref* ]]; then
        SNP_VCF_REF="$SNP_VCF_PATH"
    else
        SNP_VCF_TUM="$SNP_VCF_PATH"
    fi
done

if [[ -z "$SNP_VCF_TUM" ]]; then
  echo "Error: Could not locate hartwig tumor SNP VCF."
  exit 1
fi

echo "TUM snp vcf: ${SNP_VCF_TUM}"

cp ${SNP_VCF_PATH} $HOME/hartwig_snpfile_tum.vcf

# Get and filter germline + somatic vcf molecular pipeline
somatic_vcf=$(ls ${MOUNT_POINT_BAM}${setname}/purple/*.purple.somatic.vcf.gz)
germline_vcf=$(ls ${MOUNT_POINT_BAM}${setname}/purple/*.purple.germline.vcf.gz)

zcat "${somatic_vcf}" | grep -E '^#|REPORTED' > "$HOME/${SAMPLE_BARCODE}_somatic_filtered.vcf" || {
  echo "Error: Failed to filter somatic VCF"
  exit 1
}

zcat "${germline_vcf}" | grep -E '^#|REPORTED' > "$HOME/${SAMPLE_BARCODE}_germline_filtered.vcf" || {
  echo "Error: Failed to filter germline VCF"
  exit 1
}


# Combine hartwig+id snp vcf files and purple somatic+germline vcf files
/usr/lib/jvm/adoptopenjdk-8-hotspot-amd64/jre/bin/java -jar "$GATK" \
  -T CombineVariants \
  -R $REFERENCE \
  --variant $HOME/${SAMPLE_BARCODE}_somatic_filtered.vcf \
  --variant $HOME/${SAMPLE_BARCODE}_germline_filtered.vcf \
  --variant $HOME/$SNP_OUTPUT_VCF  \
  --variant $HOME/hartwig_snpfile_tum.vcf \
  -o $HOME/${MERGED_OUTPUT_VCF}\
  -genotypeMergeOptions UNSORTED

# Snp call column is seperate with "NA" as header, copy snp calls to tumor sample column
header_line=$(grep -m 1 '^#CHROM' "$HOME/$MERGED_OUTPUT_VCF")
IFS=$'\t' read -r -a headers <<< "$header_line"

tumor_col=-1
na_col=-1
ref_col=-1

## Identify Tumor, NA and ref column indices
for i in "${!headers[@]}"; do
  col="${headers[$i]}"
  if [[ "$col" == "NA" ]]; then
    na_col=$i
  elif [[ "$col" == *"-ref" ]]; then
    ref_col=$i
  elif [[ "$col" != "FORMAT" && "$col" != "#CHROM" && "$col" != "POS" && "$col" != "ID" && "$col" != "REF" && "$col" != "ALT" && "$col" != "QUAL" && "$col" != "FILTER" && "$col" != "INFO" ]]; then
    tumor_col=$i
  fi
done

## Check if NA and tumor column were found
if [[ $tumor_col -eq -1 || $na_col -eq -1 ]]; then
  echo "Error: Could not find both tumor and NA columns in the merged VCF."
  exit 1
fi

## Convert 0-based to 1-based indexing for awk
tumor_col_awk=$((tumor_col + 1))
na_col_awk=$((na_col + 1))
ref_col_awk=$((ref_col + 1))

## Step 1: Copy snp genotype (na) col to tumor if tumor col no variant call
awk -v tum="$tumor_col_awk" -v na="$na_col_awk" 'BEGIN { OFS="\t" }
  /^##/ { print; next }
  /^#CHROM/ { print; next }
  {
    if ($tum == "./." && $na != "./." && $na != "") {
      $tum = $na
    }
    print
  }
' "$HOME/${MERGED_OUTPUT_VCF}" > "$HOME/temp_with_na_corrected.vcf"

## Step 2: Remove ref and NA cols
awk -v ref="$ref_col_awk" -v na="$na_col_awk" 'BEGIN { OFS="\t" }
  /^##/ { print; next }
  /^#CHROM/ {
    for (i = 1; i <= NF; i++) {
      if (i != ref && i != na) {
        printf "%s%s", $i, (i == NF || (i+1 == ref || i+1 == na) ? ORS : OFS)
      }
    }
    next
  }
  {
    for (i = 1; i <= NF; i++) {
      if (i != ref && i != na) {
        printf "%s%s", $i, (i == NF || (i+1 == ref || i+1 == na) ? ORS : OFS)
      }
    }
  }
' "$HOME/temp_with_na_corrected.vcf" > "$HOME/${FINAL_OUTPUT_VCF}"

# Copy output to bucket used in share script production
gsutil cp "$HOME/${FINAL_OUTPUT_VCF}" "gs://${OUTPUT_BUCKET_NAME}/${setname}/${FINAL_OUTPUT_VCF}" || {
  echo "Error: Failed to upload VCF to output bucket"
  exit 1
}

# Unmount buckets
echo "Unmounting ${MOUNT_POINT_REFGENOME}..."
fusermount -u "${MOUNT_POINT_REFGENOME}"

echo "Unmounting ${MOUNT_POINT_BAM}..."
fusermount -u "${MOUNT_POINT_BAM}"

# Delete all temp files
rm -r "${MOUNT_POINT_BAM}"
rm -r "${MOUNT_POINT_REFGENOME}"
rm $HOME/${SAMPLE_BARCODE}_somatic_filtered.vcf
rm $HOME/${SAMPLE_BARCODE}_somatic_filtered.vcf.idx
rm $HOME/${SAMPLE_BARCODE}_germline_filtered.vcf
rm $HOME/${SAMPLE_BARCODE}_germline_filtered.vcf.idx
rm $HOME/${MERGED_OUTPUT_VCF}
rm $HOME/${FINAL_OUTPUT_VCF}
rm $HOME/${SNP_OUTPUT_VCF}
rm $HOME/hartwig_snpfile_tum.vcf
rm $HOME/temp_with_na_corrected.vcf

gsutil cp $HOME/script.log gs://${OUTPUT_BUCKET_NAME}/${setname}/${SAMPLE_BARCODE}_merge_script.log
rm $HOME/script.log