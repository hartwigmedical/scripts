#!/bin/bash
# Upload data to the storage bucket
# usage: ./upload.sh
# assuming this runs on a crunch
# Author: Andy Repton <arepton@schubergphilis.com>
##

outputBucket="hmf_bcl_storage"
apiUrl="https://api.hartwigmedicalfoundation.nl"
HMF_DWN_ID="f39de0aec3c8b5bb9d78a22ad88428ad"
profile="default"

#int-hartwig-bots
source ~/.slack
if [ -f /tmp/bcl_uploading ];
then
  echo "Uploading in progress, exiting"
  exit 0
fi

DIRECTORY="/data1/illumina_data/"

for bclName in `ls ${DIRECTORY}`
do
  if [ -f /data/sbpuploadlogs/${bclName}_SBP_Uploaded.done ];
  then
    echo "${bclName} already uploaded, skipping";
    continue;
  fi

  if [ -f ${DIRECTORY}${bclName}/RTAComplete.txt ];
  then
    find ${DIRECTORY}${bclName} -name '*.bcl.gz' | grep bcl > /dev/null
    if [[ $? -ne 0 ]];
    then
      echo "${bclName} has no BCL files, skipping";
      continue;
    fi

    sequencer=`echo ${bclName} | awk -F '_' '{print $2}'`
    index=`echo ${bclName} | awk -F '_' '{print $3}'`
    flowcell_id_tmp=`echo ${bclName} | awk -F '_' '{print $4}'`
    flowcell_id=${flowcell_id_tmp:1}
    experiment_name=`grep ExperimentName ${DIRECTORY}${bclName}/SampleSheet.csv  | awk -F , '{print $2}'`
    if [[ $? -eq 1 ]]; then
      echo "Warning! BCL is missing SampleSheet.csv"
      curl -X POST --data-urlencode 'payload={"text":"Warning! BCL is missing SampleSheet.csv, failing upload"}' ${slackChannel}
      exit 2
    fi
    echo "${bclName} is complete"

    # One upload at a time please
    touch /tmp/bcl_uploading

    # Chopping the first 7 characters off the name for the upload, as I don't know the date for the auto scheduler
    aws s3 sync ${DIRECTORY}${bclName} s3://${outputBucket}/${bclName:7} --profile ${profile} --endpoint-url https://s3.object02.schubergphilis.com --grants read=id=${HMF_DWN_ID} readacl=id=${HMF_DWN_ID} --no-follow-symlinks --endpoint-url=https://s3.object02.schubergphilis.com --exclude *fastq* --exclude *Thumbnail_Images* --exclude *Logs* --exclude *Config* --exclude *PeriodicSaveRates*
    if [[ $? -ne 0 ]]; then
      echo "something went wrong uploading the data"
      curl -X POST --data-urlencode 'payload={"text": "Something went wrong uploading BCL '${bclName}' on '$(hostname)' to object storage. I am going to try a second time"}' ${slackChannel}
      sleep 600
      aws s3 sync ${DIRECTORY}${bclName} s3://${outputBucket}/${bclName:7} --profile ${profile} --endpoint-url https://s3.object02.schubergphilis.com --grants read=id=${HMF_DWN_ID} readacl=id=${HMF_DWN_ID} --no-follow-symlinks --endpoint-url=https://s3.object02.schubergphilis.com --exclude *fastq* --exclude *Thumbnail_Images* --exclude *Logs* --exclude *Config* --exclude *PeriodicSaveRates*
      if [[ $? -ne 0 ]]; then
        curl -X POST --data-urlencode 'payload={"text": "Second attempt at uploading BCL '${bclName}' on '$(hostname)' failed. Waiting for a human to take a look"}' ${slackChannel}
        rm -f /tmp/bcl_uploading
        exit 1
      fi
    fi
    rm -f /tmp/bcl_uploading

    /usr/bin/curl -s -v \
      --cert-type pem \
      --cert /home/sbpext/bcl-upload/api.crt \
      --key /home/sbpext/bcl-upload/api.key \
      https://api.acc.hartwigmedicalfoundation.nl/hmf/v1/flowcells \
      -XPOST \
      -H "Content-Type: application/json" \
      -d '{"name": "'${experiment_name}'", "sequencer": "'${sequencer}'", "index": "'${index}'", "flowcell_id": "'${flowcell_id}'", "status": "Pending"}'

    curl -X POST --data-urlencode 'payload={"text":"BCL '${bclName}' is uploaded from '$(hostname)' to https://s3.object02.schubergphilis.com/'${outputBucket}'/'${bclName:7}' and is ready for bcl conversion"}' ${slackChannel}

    echo `date` > /data/sbpuploadlogs/${bclName}_SBP_Uploaded.done
  fi
done
