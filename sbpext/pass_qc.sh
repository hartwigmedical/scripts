#!/bin/bash


if [ "$#" -ne 3 ]; then
  echo "Useage: pass_qc api_url api_dir flowcell_id"
  exit 1
fi

api_url=$1
api_dir=$2
flowcell_id=$3

api_crt="${api_dir}/api.crt"
api_key="${api_dir}/api.key"

## Generic function to use api
hmfapi () {
  echo "$@" 1>&2
  http --ignore-stdin --cert="${api_crt}" --cert-key="${api_key}" "$@"
}

flowcell="$(hmfapi GET ${api_url}/hmf/v1/flowcells?flowcell_id=${flowcell_id})"

if [ $(echo "${flowcell}" | jq length) -eq 1 ]; then
    api_id=$(echo "${flowcell}" | jq -r .[0].id);
    hmfapi PATCH ${api_url}/hmf/v1/flowcells/${api_id} undet_rds_p_pass=True > /dev/null;
    hmfapi POST ${api_url}/hmf/v1/flowcells/${api_id}/recalculate
fi