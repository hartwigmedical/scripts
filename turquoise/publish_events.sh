#!/usr/bin/env bash

input_file="/data/ops/lims/prod/lims.json"
wd="data"

[[ $# -eq 1 ]] && input_file="$1"
[[ ! -e $wd ]] && echo "Work dir $wd does not exist, creating" && mkdir $wd

timestamp="$(date +%Y%m%d_%H%M)"

process_events() {
  field_name="$1"
  event_desc="$2"
  
  prev_events="$(ls ${wd}/${event_desc}_*.events | tail -n1) 2>/dev/null"
  [[ -e "$prev_events" ]] || touch "${wd}/${event_desc}_0.events"
  prev_events="$(ls ${wd}/${event_desc}_*.events | tail -n1)"
  curr_events="${wd}/${event_desc}_${timestamp}.events"
  
  # Could use the analysis type rather than the sample source to cut out the shallowseq for instance
  jq -cr ".samples | .[] | select(.sample_source != null) | select(.sample_source | test (\"DNA\")) | select(.${field_name} \
    | test(\"^20.*\")) | select(.sample_name | test(\"^(DRUP|WIDE|CPCT)\")) | {sample_name, ${field_name}, \"label\"}" \
    $input_file > $curr_events
  first=1
  echo "["
  diff $prev_events $curr_events | grep '^>' | sort -u | sed 's/^> //' | while read event; do
    echo "${event}" | jq -rc "[.sample_name, .${field_name}, .label] | @tsv"
  done | while read sample event_date study; do
    arr_ts="$(echo $event_date)T00:00:00.000000Z[UTC]"
    labels="$(printf "{\"name\":\"study\", \"value\":\"%s\"}" "$study")"
    subject=$(printf "{\"name\": \"%s\",\"type\": \"sample\", \"labels\": [%s]}" "$sample" "$labels")
    payload=$(printf "{\"timestamp\": \"%s\", \"type\": \"%s\", \"subjects\": [%s]}" \
      "$arr_ts" "$event_desc" "$subject")
    [[ $first -ne 1 ]] && echo ","
    [[ $first -eq 1 ]] && first=0
    echo "$payload"
  done 
  echo "]"
}

process_events arrival_date arrived | jq '.' > arrived.json
process_events report_date reported | jq '.' > reported.json
