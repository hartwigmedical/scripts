#!/usr/bin/env bash

source message_functions || exit 1

set -e

fastq_dir=$1
sample_identifier=$2
output_dir=$3

SCRIPT=$(basename "$0")
MIP_ALIGN_VERSION="1.0.0-beta.1"
MIP_COUNT_VERSION="1.0.0-beta.1"
IDENTIFIER_REGEX="^[0-9a-zA-Z-]+$"

print_usage(){
    echo ""
    echo "Description: Runs the MIP analysis for a FASTQ dataset (align + count)"
    echo "Usage: $SCRIPT <path-to-dir-with-fastq> <output-sample_identifier-identifier> <output-dir>"
    echo "Exmpl: $SCRIPT /path/to/fastq/dir FR12345678 /path/to/output/dir"
    echo ""
    exit 1
}
[[ -n "${fastq_dir}" && -n "${sample_identifier}" && -n "${output_dir}" ]] || print_usage
[[ $sample_identifier =~ $IDENTIFIER_REGEX ]] || die "Exit: unsupported chars in sample_identifier identifier [regex:${IDENTIFIER_REGEX}]"

FASTQ_DIR=$(realpath "${fastq_dir}")
OUTPUT_DIR=$(realpath "${output_dir}")
SAMPLE_IDENTIFIER="${sample_identifier}"

BAM_FILE="${SAMPLE_IDENTIFIER}_mip-align.bam"
BAM_PATH="${OUTPUT_DIR}/${BAM_FILE}"
TSV_FILE="${SAMPLE_IDENTIFIER}_mip-counts.tsv"
TSV_PATH="${OUTPUT_DIR}/${TSV_FILE}"
VCF_FILE_37="${SAMPLE_IDENTIFIER}_mip-calls-hg37.vcf"
VCF_PATH_37="${OUTPUT_DIR}/${VCF_FILE_37}"
VCF_FILE_38="${SAMPLE_IDENTIFIER}_mip-calls-hg38.vcf"
VCF_PATH_38="${OUTPUT_DIR}/${VCF_FILE_38}"

main () {
    info "Starting with ${SCRIPT}"

    info "Running mip-align [${MIP_ALIGN_VERSION}] on fastq [${FASTQ_DIR}]"
    run_mip_align

    info "Running mip-count [${MIP_COUNT_VERSION}] on bam [${BAM_FILE}]"
    run_mip_count

    info "Outputs are expected to be present at:"
    info "  ${BAM_PATH}"
    info "  ${TSV_PATH}"
    info "  ${VCF_PATH_37}"
    info "  ${VCF_PATH_38}"

    info "Finished with ${SCRIPT}"
}

run_mip_align () {
    if [[ -f "${BAM_PATH}" ]]; then
        warn "SKIPPING mip-align because BAM exists [${BAM_PATH}]"
    else
        # TODO: make output owned by user instead of root
        #docker run --user "$(id -u):$(id -g)" \
        docker run -v "${FASTQ_DIR}:/data/input" -v "${OUTPUT_DIR}:/data/output" \
          europe-west4-docker.pkg.dev/hmf-build/hmf-docker/mip-align:${MIP_ALIGN_VERSION} \
          -i "/data/input" -o "/data/output/${BAM_FILE}" -p 4 || die "Exit: mip-align has non-zero exit status"
    fi
}

run_mip_count () {
    if [[ -f "${TSV_PATH}" ]]; then
        warn "SKIPPING mip-count because TSV exists [${TSV_PATH}]"
    else
        # TODO: make output owned by user instead of root
        #docker run --user "$(id -u):$(id -g)" \
        docker run -v "${OUTPUT_DIR}:/data/input" -v "${OUTPUT_DIR}:/data/output" \
            europe-west4-docker.pkg.dev/hmf-build/hmf-docker/mip-count:${MIP_COUNT_VERSION} \
            -i "/data/input/${BAM_FILE}" -o "/data/output" -t "${SAMPLE_IDENTIFIER}" || die "Exit: mip-count has non-zero exit status"
    fi
}

main