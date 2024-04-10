#!/usr/bin/env bash

source message_functions || exit 1

fastq_dir_old=${1:-./bcl2fastq/old/conversion}
fastq_dir_new=${2:-./bcl2fastq/new/conversion}

info "Verifying bcl2fastq FASTQ output"
info "  FASTQ dir OLD: ${fastq_dir_old}"
info "  FASTQ dir NEW: ${fastq_dir_new}"

for sample_index in {1..8}; do
    for read_index in {1..2}; do
        # Note: Only works with ISEQ testdata with one lane!
        file_regex="S${sample_index}_L001_R${read_index}_001.fastq.gz"
        old_file=$(find "${fastq_dir_old}/" -name "*_${file_regex}")
        new_file=$(find "${fastq_dir_new}/" -name "*_${file_regex}")
        [[ -f "${old_file}" ]] || die "Unable to find file [${old_file}]"
        [[ -f "${new_file}" ]] || die "Unable to find file [${new_file}]"
        info "Diffing file [${file_regex}]"
        info "  OLD: ${old_file}"
        info "  NEW: ${new_file}"
        diff_count=$(diff <(zcat "${old_file}" | grep -v ^@) <(zcat "${new_file}" | grep -v ^@) | wc -l)
        if [[ "$diff_count" -ne 0 ]]; then
            warn "Diff was not empty for file [${file_regex}]!!!"
        fi
    done
done
