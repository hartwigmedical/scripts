#!/usr/bin/env bash

source message_functions || exit 1

sample=$1

if [[ -z ${sample} ]]; then
    error "No sample name provided"
fi

run_id=$( hmf_api_get "reports/created" | jq '.[] | select(.sample_name == "'${sample}'") | .run_id' | tail -1 )
set_name=$( hmf_api_get 'runs/'${run_id}'' | jq -r '.set | .name' )

mkdir temp
mkdir temp_cup


echo ""
echo "Do you want to upload the final OncoAct report and json? y or n"
read answer_report
if [[ ${answer_report} == "y" ]]; then
    gsutil cp gs://patient-reporter-final-prod-1/${sample}*report.pdf ~/temp
    gsutil cp gs://patient-reporter-final-prod-1/${sample}*oncoact.json ~/temp
    upload_file_to_nc_for_viewing ~/temp/${sample}*report.pdf
    upload_file_to_nc_for_viewing ~/temp/${sample}*.json
    echo "The OncoAct report + json for ${sample} are uploaded to nextcloud STAGING/New-Reports-Viewing"
fi

echo ""
echo "Do you want to upload the ORANGE report? y or n"
read answer_orange
if [[ ${answer_orange} == "y" ]]; then
    gsutil cp gs://diagnostic-pipeline-output-prod-1/${set_name}/orange/${sample}*.orange.pdf ~/temp
    upload_file_to_nc_for_viewing ~/temp/${sample}*orange.pdf
    echo "The ORANGE report for ${sample} is uploaded to nextcloud STAGING/New-Reports-Viewing"
fi

echo ""
echo "Do you want to upload the CUPPA RUO report? y or n"
read answer_cuppa_ruo
if [[ ${answer_cuppa_ruo} == "y" ]]; then
    gsutil cp gs://diagnostic-pipeline-output-prod-1/${set_name}/cuppa/${sample}*_cup_report.pdf ~/temp_cup
    upload_file_to_nc_for_viewing ~/temp_cup/${sample}*_cup_report.pdf
    echo "The CUPPA RUO report for ${sample} is uploaded to nextcloud STAGING/New-Reports-Viewing"
fi

rm -r ~/temp/ 2>&1
rm -r ~/temp_cup/ 2>&1