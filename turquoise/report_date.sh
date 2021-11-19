#!/usr/bin/env bash

input_file="/data/ops/lims/prod/sharing_db_generated.tsv"
wd="$HOME/.turquoise/data"
out_dir="$HOME/.turquoise"

[[ $# -eq 1 ]] && input_file="${1}"
[[ ! -e ${wd} ]] && echo "Work dir $wd does not exist, creating" && mkdir ${wd}

timestamp="$(date +%Y%m%d_%H%M)"

process_events() {
    event_desc="reported"
    curr_events="${wd}/${event_desc}_${timestamp}.events"

    cat ${1} | awk '$4 ~ "-20" {print $2, $4, $5}' | while read sample date report_type; do
        echo "${sample} $(date -d "${date}" "+%Y-%m-%d")T00:00:00.000000Z[UTC]" ${report_type} >> ${curr_events}
    done
            
    first=1
    echo "["
    cat ${curr_events} | while read sample date report_type; do
        [[ ${first} -ne 1 ]] && echo ","
        [[ ${first} -eq 1 ]] && first=0
        final_date="$(awk -v sample="${sample}" '$1 == sample {print $2}' < ${curr_events} | tail -n1)"
        report_type="$(awk -v sample="${sample}" '$1 == sample {print $3}' < ${curr_events} | tail -n1)"
        labels=$(printf "{\"name\": \"report_type\",\n\"value\": \"%s\"\n}" "${report_type}")
        subject="$(printf "{\"name\": \"${sample}\", \"type\": \"sample\"}")"
        printf "{\"timestamp\": \"%s\", \"type\": \"reported\", \"labels\": [%s], \"subjects\": [%s]}" "${final_date}" "$labels" "${subject}" 
    done
    echo "]"
}

process_events ${input_file} | jq '.' > ${out_dir}/reported.json
