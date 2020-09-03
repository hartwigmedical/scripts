#!/usr/bin/env bash

source api_functions
source gcp_functions

set=$1 && shift

if [[ -z "${set}" ]]; then
    echo "[ERROR] No set provided to $(basename $0). Exiting"
    exit 1
fi

run_bucket=$(load_intial_run_bucket_for_set ${set})
if [[ -z "${run_bucket}" ]]; then
    echo "[ERROR] No initial run bucket found for set '${set}'. Exiting"
    exit 1
fi

echo "[INFO] Copying set ${set} from ${run_bucket} to hmf-crunch"

switch_to_hmf_ops_service_account
echo "gsutil -u hmf-database -m cp gs://${run_bucket}/${set}/sage/* gs://hmf-sage/$set/"
echo "gsutil -u hmf-database -m cp gs://${run_bucket}/${set}/amber/* gs://hmf-amber/$set/"
echo "gsutil -u hmf-database -m cp gs://${run_bucket}/${set}/cobalt/* gs://hmf-cobalt/$set/"
echo "gsutil -u hmf-database -m cp gs://${run_bucket}/${set}/gridss/*unfiltered.vcf.gz* gs://hmf-gridss/unfiltered/$set/"
echo "gsutil -u hmf-database -m cp gs://${run_bucket}/${set}/gridss/*somatic*.vcf.gz* gs://hmf-gripss/$set/"