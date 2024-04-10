#!/usr/bin/env bash

source message_functions || exit 1

old_path="./bcl2fastq/old/forensics"
new_path="./bcl2fastq/new/forensics"

info "Comparing forensics output"

forensics_files=(Adapter_Metrics.csv Demultiplex_Stats.csv fastq_list.csv Index_Hopping_Counts.csv Quality_Metrics.csv RunInfo.xml SampleSheet.csv Top_Unknown_Barcodes.csv)

for forensics_file in "${forensics_files[@]}"; do
    old_file="${old_path}/output/results/Reports/${forensics_file}"
    new_file="${new_path}/output/results/Reports/${forensics_file}"
    info "Diffing file [${forensics_file}]"
    info "  OLD: ${old_file}"
    info "  NEW: ${new_file}"
    diff_count=$(diff "${old_file}" "${new_file}" | wc -l)
    if [[ "$diff_count" -ne 0 ]]; then
        warn "Diff was not empty for file [${forensics_file}]:"
    diff "${old_file}" "${new_file}"
    fi
done
