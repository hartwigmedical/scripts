#!/bin/bash

## Desc: Uploads flowcell BCL level data to the storage bucket

## TODO move slack code to all servers that will need it
#source /data/common/dbs/slack/slack_functions.sh

## server specific flowcell location
flowcells_dir="/data1/illumina_data"

bucket="bcl-input-prod"
gcp_cred="/home/sbpext/bcl-upload-prod/bcl-input-prod.json"

if [ "$#" -ne 1 ]; then
  echo "Useage: upload_gcp api_url api_dir"
  exit 1
fi

api_url=$1
api_dir=$2

api_crt="${api_dir}/api.crt"
api_key="${api_dir}/api.key"
logs_dir="/data/sbpuploadlogsgcp"
hostname=$(hostname)

## tmp file to avoid multiple uploads
## keep this the same as in SBP script to avoid sync to s3 and gcp at same time
tmp_upload_file="/tmp/bcl_uploading"

## Generic function to use api
hmfapi () {
  echo $@ 1>&2
  http --ignore-stdin --cert="${api_crt}" --cert-key="${api_key}" $@
}

## Generic function that posts a message to slack channel
slack () {
  local msg=$1 && shift
  hmf_slack_default "${msg}" > /dev/null
}

## Run once at a time
if [ -f ${tmp_upload_file} ]; then
  echo "Uploading in progress, exiting" && exit 0
else
  touch "${tmp_upload_file}"
fi

## Find flowcells
for flowcell_path in $(find ${flowcells_dir} -mindepth 1 -maxdepth 1 -type d)
do
  flowcell_name=$( basename ${flowcell_path} )
  samplesheet_path="${flowcell_path}/SampleSheet.csv"
  rtacomplete_path="${flowcell_path}/RTAComplete.txt"
  runinfoxml_path="${flowcell_path}/RunInfo.xml"
  log_file="${logs_dir}/${flowcell_name}_SBP_Uploaded.done"

  if [ -f "$log_file" ]; then
    echo "[INFO] Flowcell ${flowcell_name} already uploaded: skipping"
    continue
  fi

  if [ ! -f "${samplesheet_path}" ]; then
    echo "[WARN] Flowcell ${flowcell_name} on $hostname is missing SampleSheet (${samplesheet_path}): skipping"
      continue
  fi

  flowcell_id=$(xmllint --xpath "/RunInfo/Run/Flowcell/text()" ${runinfoxml_path})
  sequencer=$(echo ${flowcell_name} | cut -d '_' -f2)
  index=$(echo ${flowcell_name} | cut -d '_' -f3)
  experiment_name=$(grep ExperimentName ${samplesheet_path} | cut -d "," -f2)

  ## -----
  ## TESTING TODO: remove test stuff
  ## -----
  test_string="_Test200130"
  flowcell_id+="$test_string"
  experiment_name+="$test_string"

  echo "[INFO] Working on path \"${flowcell_path}\""
  echo "[INFO]   Hostname: ${hostname}"
  echo "[INFO]   Samplesheet: ${samplesheet_path}"
  echo "[INFO]   Flowcell ID: ${flowcell_id}"
  echo "[INFO]   Flowcell Name: ${flowcell_name}"
  echo "[INFO]   Experiment Name: ${experiment_name}"

  if [ -z $flowcell_id ]; then
    echo "[WARN] Flowcell ${flowcell_name} on $hostname is missing flowcell_id: skipping"
    continue
  fi

  if [ ! -f ${rtacomplete_path} ]; then
    # Sequencing not finished yet but do pre-register
    echo "[INFO] Flowcell ${flowcell_name} has no run complete file yet (${rtacomplete_path}): skipping"
    flowcell="$(hmfapi GET ${api_url}/hmf/v1/flowcells?flowcell_id=${flowcell_id})"
    if [ $(echo "${flowcell}" | jq length) -eq 0 ]; then
      hmfapi POST ${api_url}/hmf/v1/flowcells name="${experiment_name}" sequencer="${sequencer}" index="${index}" flowcell_id="${flowcell_id}" status=Sequencing > /dev/null;
       fi
    continue
  fi

  echo "[INFO] Flowcell ${flowcell_name} is ready for uploading BCL"
  find ${flowcell_path} -name '*.bcl.gz' -or -name '*.cbcl' -or -name '*.bcl.bgzf' | grep bcl > /dev/null
  if [[ $? -ne 0 ]]; then
    #slack "${flowcell_name} on ${hostname} is ready but has no BCL files, this probably means that the sequencer has failed"
     echo "${flowcell_name} on ${hostname} is ready but has no BCL files, this probably means that the sequencer has failed"
    continue
  fi

  gcloud auth activate-service-account --key-file ${gcp_cred}
  cmd="gsutil -m rsync -r ${flowcell_path} gs://${bucket}/${flowcell_name}"
  echo "[INFO] Executing rsync ($cmd)"
  $cmd

  if [[ $? -ne 0 ]]; then
    #slack "[WARN] Something wrong with BCL upload \"${flowcell_path}\" on \"${hostname}\""
    echo "[WARN] Something wrong with BCL upload \"${flowcell_path}\" on \"${hostname}\""
    continue
  fi

  flowcell="$(hmfapi GET ${api_url}/hmf/v1/flowcells?flowcell_id=${flowcell_id})"

  if [ $(echo "${flowcell}" | jq length) -eq 1 ]; then
    api_id=$(echo "${flowcell}" | jq -r .[0].id);
    status=$(echo "${flowcell}" | jq -r .[0].status);

    echo "[INFO] Found existing flowcell at ${api_url} (api_id:${api_id}|status=${status})"

    if [ "${status}" == "Testing" ]; then
      echo "[INFO] Flowcell status is Testing"
    elif [ "${status}" == "Sequencing" ]; then
      echo "[INFO] PATCH flowcell status to Pending (fcid:${flowcell_id}|api_id:${api_id}|host:${hostname})"
      hmfapi PATCH ${api_url}/hmf/v1/flowcells/${api_id} status=Pending index="${index}" sequencer="${sequencer}" bucket="${bucket}" > /dev/null;
    else
      #slack "[WARN] Would try invalid flowcell transition from ${status} to Pending, check flowcell (fcid:${flowcell_id}|api_id:${api_id}|host:${hostname})";
       echo "[WARN] SKIPPING. Would do invalid flowcell transition from ${status} to Pending, check flowcell (fcid:${flowcell_id}|api_id:${api_id}|host:${hostname})";
      continue
    fi
  else
    echo "[INFO] POST flowcell (${flowcell_id}) to ${api_url} and set pending because BCL has been uploaded"
    hmfapi POST ${api_url}/hmf/v1/flowcells name="${experiment_name}" sequencer="${sequencer}" index="${index}" flowcell_id="${flowcell_id}" bucket="${bucket}" status=Pending
  fi

  echo "[INFO] BCL of flowcell $flowcell_name is uploaded from $hostname and is ready for bcl2fastq conversion"
  echo `date` > /data/sbpuploadlogs/${flowcell_name}_SBP_Uploaded.done

done

echo "[INFO] Removing $tmp_upload_file"
rm -f ${tmp_upload_file}
echo "[INFO] Finished with $0"