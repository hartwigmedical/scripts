#!/usr/bin/env bash

echo "--1. Reading barcodes from silo  --"

token=$(gcloud auth print-identity-token)

rm -f available_barcodes.tsv
curl --oauth2-bearer $token -s http://linkage-silo.prod-1/personal-data/known | jq -r '.[].isolationBarcode' > available_barcodes.tsv
wc -l available_barcodes.tsv

echo "--1. Done --"


rm -f linkage_file.tsv
echo -e "isolationBarcode\tcohort\tallowInternalUse\tsampleId\tsampleBarcode\tisolationBarcode1\thartwigNumber\ttype\thartwigSpecimenId\tpathologyNumer\thartwigDonorId\tisolationBarcode2\tinitials\tbirthSurname\tcurrentSurname\tbirthDate\tgender\tpostalCode\thospitalName\thospitalPatientId\tcreatedAt" > linkage_file.tsv


echo "--2. Getting cohort/consent status from LAMA and getting information from Linkage silo --"

while read barcode; do

cohort=$( lama_api_get queries/patient-reporter/isolation-barcode/${barcode} | jq -r '.cohort')
consent=$( lama_api_get queries/consent/isolation-barcode/${barcode} | jq -r '.allowInternalUse')

echo ${barcode} ${cohort} ${consent}
if [[ ${cohort} == "GENAYA" ]]; then
	echo " .. Relevant cohort thus getting information from linkage silo .. "
	hartwig=$( curl --oauth2-bearer $token -s http://linkage-silo.prod-1/tumor-sample/isolation-barcode/${barcode} | jq -r 'to_entries|map(.value)|@tsv' )
	personal=$( curl --oauth2-bearer $token -s http://linkage-silo.prod-1/personal-data/isolation-barcode/${barcode} | jq -r 'to_entries|map(.value)|@tsv' )
	printf "\n ${barcode} \t ${cohort} \t ${consent} \t ${hartwig} \t ${personal}" >> linkage_file.tsv
fi
echo ""

done < available_barcodes.tsv

echo "--2. Done --"
wc -l linkage_file.tsv