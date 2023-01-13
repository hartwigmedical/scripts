#!/usr/bin/env bash

source message_functions || exit 1
source locate_reporting_api || exit 1

sampleId=$1

if [[ -z ${sampleId} ]]; then
    error "No sample name provided"
fi

echo "--TAT information for sample ${sampleId}--"


#### Last material arrival date ###

temp_nr=$( lama patients/tumorsamples $sampleId | head -n1 | grep -o '^.*arrivalHmf' | wc -w )
tumor_arrival=$( lama patients/tumorsamples $sampleId | awk '{ print $'${temp_nr}' }' | grep -v arrivalHmf )

sampleId_ref=$( echo $(echo ${sampleId} | cut -c1-12)R )
temp_nr=$( lama patients/bloodsamples $sampleId_ref | head -n1 | grep -o '^.*arrivalHmf' | wc -w )
ref_arrival=$( lama patients/bloodsamples $sampleId_ref | awk '{ print $'${temp_nr}' }' | grep -v arrivalHmf )

if [ $(echo $tumor_arrival | wc -w) == 0 ]; then
	echo "Tumor_material_not_yet_arrived"
	exit
fi

if [ $(echo $ref_arrival | wc -w) == 0 ]; then
	echo "Reference_material_not_yet_arrived"
	exit
fi


tumor_arrival_sec=$( date --date=$tumor_arrival +%s )
ref_arrival_sec=$( date --date=$ref_arrival +%s )

if [[ $tumor_arrival_sec -gt $ref_arrival_sec ]]; then
	start_sec=$tumor_arrival_sec
else
    start_sec=$ref_arrival_sec
fi



##### Determine whether sample was already reported

barcode=$( hmf_api_get samples?name=${sampleId} | jq -r .[].barcode )
report_created_id=$(extract_most_recent_reporting_id_on_barcode $barcode )
reported=$( hmf_api_get reports/shared?report_created_id=${report_created_id}  | jq .[] | jq -r '.share_time' | tr 'T' ' ' | sed 's/\s.*$//')

if [ $(echo $reported | wc -w ) == 0 ]; then
    ##### If not show the current TAT
    now_sec=$(date +%s)
    echo "This sample is still in process!"
    echo "The current TAT is:"
    echo $(( ($now_sec - $start_sec )/(60*60*24) ))
else
    ##### If yes show reporting info
    echo "This OncoAct has been reported on:"
    echo $reported
    echo "The TAT of reporting was:"
    reported_sec=$( date --date=$reported +%s )
    echo $(( ($reported_sec - $start_sec )/(60*60*24) ))
    exit
fi




