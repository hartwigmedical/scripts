#!/usr/bin/env bash

source message_functions || exit 1

set -e

fastq_dir=$1
sample=$2
output_dir=$3

MIP_ALIGN_VERSION="1.0.0-beta.1"
MIP_COUNT_VERSION="1.0.0-beta.1"
IDENTIFIER_REGEX="^[0-9a-zA-Z-]+$"

BAM_FILE="${sample}_mip-align.bam"
BAM_PATH="${output_dir}/${BAM_FILE}"
TSV_FILE="${sample}_mip-counts.tsv"
TSV_PATH="${output_dir}/${TSV_FILE}"
SCRIPT=$(basename "$0")

print_usage(){
    echo ""
    echo "Description: Runs the MIP analysis for a FASTQ dataset (align + count)"
    echo "Usage: $SCRIPT <path-to-dir-with-fastq> <output-sample-identifier> <output-dir>"
    echo "Exmpl: $SCRIPT /path/to/fastq/dir FR12345678 /path/to/output/dir"
    echo ""
    exit 1
}
[[ -n "${fastq_dir}" && -n "${sample}" && -n "${output_dir}" ]] || print_usage
[[ $sample =~ $IDENTIFIER_REGEX ]] || die "Exit: unsupported chars in sample identifier [regex:${IDENTIFIER_REGEX}]"

main () {
    info "Starting with ${SCRIPT}"
    run_mip_align "${BAM_PATH}" "${fastq_dir}" "${output_dir}" "${BAM_FILE}"
    run_mip_count "${TSV_PATH}" "${BAM_FILE}" "${sample}"
    info "Finished with ${SCRIPT}"
}

run_mip_align () {
    local bam_path=$1 && shift
    local fastq_dir=$1 && shift
    local output_dir=$1 && shift
    local bam_file=$1 && shift
    info "RUNNING mip-align [${MIP_ALIGN_VERSION}] for fastq in [${fastq_dir}] to bam [${bam_file}]"
    if [[ -f "${bam_path}" ]]; then
        warn "SKIPPING mip-align because BAM exists [${bam_path}]"
    else
        #docker run --user "$(id -u):$(id -g)" \
        docker run -v "${fastq_dir}:/data/input" -v "${output_dir}:/data/output" \
          europe-west4-docker.pkg.dev/hmf-build/hmf-docker/mip-align:${MIP_ALIGN_VERSION} \
          -i "/data/input" -o "/data/output/${bam_file}" -p 4 || die "Exit: mip-align has non-zero exit status"
    fi
}

run_mip_count () {
    local tsv_path=$1 && shift
    local bam_file=$1 && shift
    local sample=$1 && shift
    info "RUNNING mip-count [${MIP_COUNT_VERSION}] for bam [${bam_file}] to counts [${tsv_path}]"
    if [[ -f "${tsv_path}" ]]; then
        warn "SKIPPING mip-count because TSV exists [${tsv_path}]"
    else
        #docker run --user "$(id -u):$(id -g)" \
        docker run -v "${output_dir}:/data/input" -v "${output_dir}:/data/output" \
            europe-west4-docker.pkg.dev/hmf-build/hmf-docker/mip-count:${MIP_COUNT_VERSION} \
            -i "/data/input/${bam_file}" -o "/data/output" -t "${sample}" || die "Exit: mip-count has non-zero exit status"
    fi
}

main