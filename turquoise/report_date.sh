#!/usr/bin/env bash

input_file="/data/ops/lims/prod/reporting_db.tsv"
wd="data"

[[ $# -eq 1 ]] && input_file="${1}"
[[ ! -e ${wd} ]] && echo "Work dir $wd does not exist, creating" && mkdir ${wd}

timestamp="$(date +%Y%m%d_%H%M)"

process_events() {
    event_desc="reported"
    prev_events="$(ls ${wd}/${event_desc}_*.events | tail -n1) 2>/dev/null"
    [[ -e "$prev_events" ]] || touch "${wd}/${event_desc}_0.events"
    prev_events="$(ls ${wd}/${event_desc}_*.events | tail -n1)"
    curr_events="${wd}/${event_desc}_${timestamp}.events"

    cat $1 | awk -F '\t' '$4 ~ '-202' {print $2, $4}' | sort -u > ${curr_events}

    first=1
    echo "["
    diff ${prev_events} ${curr_events} | grep '^>' | sort -u | sed 's/^> //' | while read sample date; do
        [[ ${first} -ne 1 ]] && echo ","
        [[ ${first} -eq 1 ]] && first=0
        ts="$(date -d "${date}" "+%Y-%m-%d")T00:00:00:000000Z[UTC]"
        subject="$(printf "{\"name\": \"${sample}\", \"type\": \"sample\"}")"
        printf "{\"timestamp\": \"%s\", \"type\": \"reported\", \"subjects\": [%s]}" "${ts}" "${subject}" 
    done
    echo "]"
}

process_events ${input_file} | jq '.' > reported.json
