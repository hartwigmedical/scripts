# Get dates per sample

source message_functions || exit 1
source locate_reporting_api || exit 1
source lims_functions || exit 1

barcode=$1

sampleId=$( find_name_for_barcode ${barcode} )

temp_nr=$( lama patients/tumorsamples ${sampleId} | head -n1 | grep -o '^.*arrivalHmf' | wc -w )
tumor_arrival=$( lama patients/tumorsamples ${sampleId} | awk '{ print $'${temp_nr}' }' | grep -v arrivalHmf )

sampleId_ref=$( echo $(echo ${sampleId} | cut -c1-12)R )
temp_nr=$( lama patients/bloodsamples ${sampleId_ref} | head -n1 | grep -o '^.*arrivalHmf' | wc -w )
ref_arrival=$( lama patients/bloodsamples ${sampleId_ref} | awk '{ print $'${temp_nr}' }' | grep -v arrivalHmf )


if [ $(echo ${tumor_arrival} | wc -w) == 0 ]; then
    tumor_arrival="NA"
fi

if [ $(echo ${ref_arrival} | wc -w) == 0 ]; then
    ref_arrival="NA"
fi

report_created_id=$( extract_first_time_reporting_id_on_barcode $barcode )
reported=$( hmf_api_get reports/shared?report_created_id=${report_created_id}  | jq .[] | jq -r '.share_time' | tr 'T' ' ' | sed 's/\s.*$//' )

if [ $(echo $reported | wc -w ) == 0 ]; then
    reported="NA"
fi

echo ${sampleId}
echo ${tumor_arrival}
echo ${ref_arrival}
echo ${reported}