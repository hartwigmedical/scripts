#!/usr/bin/env bash

input_file="/data/ops/lims/prod/lims.json"
wd="$HOME/.turquoise/data"
tmp_file="$(mktemp -p ${wd} arrived.XXXXX)"
ref_samples="${wd}/reference_samples.tmp"

out_dir="$HOME/.turquoise"
ref_out="${out_dir}/reference_arrived.json"
tum_out="${out_dir}/arrived.json"

[[ $# -eq 1 ]] && input_file="$1"
[[ ! -e $wd ]] && echo "Work dir ${wd} does not exist, creating" && mkdir -p ${wd}

timestamp="$(date +%Y%m%d_%H%M)"

process_events() {
    field_name="$1"
    event_desc="$2"
  
    prev_events="$(ls ${wd}/${event_desc}_*.events | tail -n1) 2>/dev/null"
    [[ -e "$prev_events" ]] || touch "${wd}/${event_desc}_0.events"
    prev_events="$(ls ${wd}/${event_desc}_*.events | tail -n1)"
    curr_events="${wd}/${event_desc}_${timestamp}.events"

    rm -f "${ref_samples}"
 
    jq -cr ".samples | .[] | select(.sample_source != null) | select(.ref_sample_id != null) | select(.ref_sample_id | test(\".+\")) \
        | select(.sample_source | test (\"DNA\")) | select(.${field_name} | test(\"^2020.*\")) \
        | select(.sample_name | test(\"^(DRUP|WIDE|CPCT)\")) | {sample_name, ${field_name}, \"label\", ref_sample_id}" \
        ${input_file} > ${curr_events}
    first=1
    echo "[" > "${tum_out}"
    diff ${prev_events} ${curr_events} | grep '^>' | sort -u | sed 's/^> //' | while read event; do
        echo "${event}" | jq -rc "[.sample_name, .${field_name}, .label, .ref_sample_id] | @tsv"
    done | while read sample event_date study reference_sample; do
        arr_ts="$(echo ${event_date})T00:00:00.000000Z[UTC]"
        labels="$(printf "{\"name\":\"study\", \"value\":\"%s\"}" "$study")"
        subject=$(printf "{\"name\": \"%s\",\"type\": \"sample\", \"labels\": [%s]}" "${sample}" "${labels}")
        payload=$(printf "{\"timestamp\": \"%s\", \"type\": \"tumor.arrived\", \"subjects\": [%s]}" \
            "${arr_ts}" "${subject}")
        [[ ${first} -ne 1 ]] && echo ","
        [[ ${first} -eq 1 ]] && first=0
        echo "${payload}"
        echo "${sample} ${reference_sample} ${study}" >> "${ref_samples}"
    done >> "${tum_out}"
    echo "]" >> "${tum_out}"
  
    ref_complete="${wd}/reference_samples_complete.tmp" 
    rm -f "${ref_complete}"
    cat "${ref_samples}" | while read tumour_sample reference_sample_key study; do
      echo "${tumour_sample} $(jq -cr ".samples.${reference_sample_key} | [.arrival_date] | @tsv" ${input_file}) ${study}" >> "${ref_complete}"
    done

    echo "[" > "${ref_out}"
    first=1 
    cat "${ref_complete}" | while read sample_name event_date study; do
        arr_ts="$(echo ${event_date})T00:00:00.000000Z[UTC]"
        labels="$(printf "{\"name\":\"study\", \"value\":\"%s\"}" "$study")"
        subject=$(printf "{\"name\": \"%s\",\"type\": \"sample\", \"labels\": [%s]}" "${sample_name}" "${labels}")
        payload=$(printf "{\"timestamp\": \"%s\", \"type\": \"reference.arrived\", \"subjects\": [%s]}" "${arr_ts}" "${subject}")
        [[ ${first} -ne 1 ]] && echo "," >> "${ref_out}"
        [[ ${first} -eq 1 ]] && first=0
        echo "${payload}" >> "${ref_out}"
    done
    echo "]" >> "${ref_out}"
}

process_events arrival_date arrived 

jq '.' "${tum_out}" > arrived.pp && mv arrived.pp "${tum_out}"
jq '.' "${ref_out}" > ref.pp && mv ref.pp "${ref_out}"
