#!/usr/bin/env bash

input_file="/data/ops/lims/prod/lims.json"
wd="$HOME/.turquoise/data"
out_file="$HOME/.turquoise/lims.json"

[[ $# -eq 1 ]] && input_file="$1"
[[ ! -e $wd ]] && echo "Work dir ${wd} does not exist, creating" && mkdir -p ${wd}

timestamp="$(date +%Y%m%d_%H%M)"

prev="$(ls ${wd}/lims_*.json 2>/dev/null | tail -n1)"
if [[ ! -e "$prev" ]]; then
    cp ${input_file} ${out_file}
else
    diff <(jq '.samples' ${prev}) <(jq '.samples' ${input_file}) | grep '>' | sed 's/^>//' | sed -z 's/^/\{/' \
        | sed -z 's/,\([^,]*\)$/\}\1/' | jq '.' > ${out_file}
fi
cp ${input_file} ${wd}/lims_${timestamp}.json

