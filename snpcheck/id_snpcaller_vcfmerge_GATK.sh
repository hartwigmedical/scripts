#!/bin/bash

set -euo pipefail

LOCAL_INPUT_DIR="$HOME/input_files"
mkdir -p "$LOCAL_INPUT_DIR"
exec > "$LOCAL_INPUT_DIR/snpcaller_merge.log" 2>&1

setname="$1"
INPUT_DIR="gs://wgs-combined-snps-vcfs/${setname}/inputfiles"
SAMPLE_BARCODE=$(gsutil cat "$INPUT_DIR/sample_barcode.txt")
CONVERTED_REPORTING_ID=$(gsutil cat "$INPUT_DIR/converted_reporting_id.txt")
OUTPUT_BUCKET_NAME=$(gsutil cat "$INPUT_DIR/output_bucket.txt")

# define BAM_TUM path
BAM_TUM=#add location
REFERENCE=#add location

# Copy input files to local
gsutil cp "$INPUT_DIR/id_snps_intervals.hg37.bed" "$LOCAL_INPUT_DIR/id_snps_intervals.hg37.bed"
gsutil cp "$INPUT_DIR/hartwig_snpfile_tum.vcf" "$LOCAL_INPUT_DIR/hartwig_snpfile_tum.vcf"
gsutil cp "$INPUT_DIR/*.purple.germline.vcf.gz" "$LOCAL_INPUT_DIR/"
gsutil cp "$INPUT_DIR/*.purple.somatic.vcf.gz" "$LOCAL_INPUT_DIR/"
gsutil cp "$INPUT_DIR/GenomeAnalysisTK.jar" "$LOCAL_INPUT_DIR/GenomeAnalysisTK.jar"
gsutil cp "$INPUT_DIR/java" "$LOCAL_INPUT_DIR/java"

# Create input file vars
INTERVALS="$LOCAL_INPUT_DIR/id_snps_intervals.hg37.bed"
SNP_VCF_TUM="$LOCAL_INPUT_DIR/hartwig_snpfile_tum.vcf"
GERMLINE_VCF=$(ls $LOCAL_INPUT_DIR/*.purple.germline.vcf.gz | head -n 1)
SOMATIC_VCF=$(ls $LOCAL_INPUT_DIR/*.purple.somatic.vcf.gz | head -n 1)
GATK="$LOCAL_INPUT_DIR/GenomeAnalysisTK.jar"
JAVA="$LOCAL_INPUT_DIR/java"

# Create output file name vars
SNP_OUTPUT_VCF="${LOCAL_INPUT_DIR}/${SAMPLE_BARCODE}_snp_genotype_output.vcf"
MERGED_OUTPUT_VCF="${LOCAL_INPUT_DIR}/${SAMPLE_BARCODE}_merged.vcf"
FINAL_OUTPUT_VCF="${LOCAL_INPUT_DIR}/${CONVERTED_REPORTING_ID}.reported.variants.and.snps.vcf"

# GATK HaplotypeCaller (UnifiedGenotyper)
$JAVA -Xmx20G -jar "$GATK" \
  -T UnifiedGenotyper \
  -nct $(nproc) \
  --input_file "$BAM_TUM" \
  -o "$SNP_OUTPUT_VCF" \
  -L "$INTERVALS" \
  --reference_sequence "$REFERENCE" \
  --output_mode EMIT_ALL_SITES

# Filter somatic/germline VCFs
zcat "$SOMATIC_VCF" | grep -E '^#|REPORTED' > "${LOCAL_INPUT_DIR}/${SAMPLE_BARCODE}_somatic_filtered.vcf"
zcat "$GERMLINE_VCF" | grep -E '^#|REPORTED' > "${LOCAL_INPUT_DIR}/${SAMPLE_BARCODE}_germline_filtered.vcf"

# CombineVariants
$JAVA -jar "$GATK" \
  -T CombineVariants \
  -R "$REFERENCE" \
  --variant "${LOCAL_INPUT_DIR}/${SAMPLE_BARCODE}_somatic_filtered.vcf" \
  --variant "${LOCAL_INPUT_DIR}/${SAMPLE_BARCODE}_germline_filtered.vcf" \
  --variant "$SNP_OUTPUT_VCF" \
  --variant "$SNP_VCF_TUM" \
  -o "$MERGED_OUTPUT_VCF" \
  -genotypeMergeOptions UNSORTED

# Copy NA column to tumor if missing
header_line=$(grep -m 1 '^#CHROM' "$MERGED_OUTPUT_VCF")
IFS=$'\t' read -r -a headers <<< "$header_line"

tumor_col=-1
na_col=-1
ref_col=-1
for i in "${!headers[@]}"; do
  col="${headers[$i]}"
  [[ "$col" == "NA" ]] && na_col=$i
  [[ "$col" == *"-ref" ]] && ref_col=$i
  [[ "$col" != "FORMAT" && ! "$col" =~ ^#?(CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO)$ ]] && tumor_col=$i
done

(( tumor_col_awk = tumor_col + 1 ))
(( na_col_awk = na_col + 1 ))
(( ref_col_awk = ref_col + 1 ))

awk -v tum="$tumor_col_awk" -v na="$na_col_awk" 'BEGIN { OFS="\t" }
  /^##/ { print; next }
  /^#CHROM/ { print; next }
  {
    if ($tum == "./." && $na != "./." && $na != "") {
      $tum = $na
    }
    print
  }
' "$MERGED_OUTPUT_VCF" > "${LOCAL_INPUT_DIR}/temp_with_na_corrected.vcf"

awk -v ref="$ref_col_awk" -v na="$na_col_awk" 'BEGIN { OFS="\t" }
  /^##/ { print; next }
  /^#CHROM/ {
    for (i = 1; i <= NF; i++) if (i != ref && i != na) printf "%s%s", $i, (i == NF ? ORS : OFS)
    next
  }
  {
    for (i = 1; i <= NF; i++) if (i != ref && i != na) printf "%s%s", $i, (i == NF ? ORS : OFS)
  }
' "${LOCAL_INPUT_DIR}/temp_with_na_corrected.vcf" > "$FINAL_OUTPUT_VCF"

gsutil cp "${FINAL_OUTPUT_VCF}" "gs://${OUTPUT_BUCKET_NAME}/${setname}/${FINAL_OUTPUT_VCF}" || {
  echo "Error: Failed to upload VCF to output bucket"
  exit 1
}

gsutil cp "$LOCAL_INPUT_DIR/snpcaller_merge.log" "gs://${OUTPUT_BUCKET_NAME}/${setname}/"

gsutil rm -r "${INPUT_DIR}"
gsutil rm -r "${LOCAL_INPUT_DIR}"

echo "Final VCF copied to output bucket: gs://${OUTPUT_BUCKET_NAME}/${setname}/${FINAL_OUTPUT_VCF} "

