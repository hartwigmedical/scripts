#!/bin/bash

function main {

  echo ''
  echo '--- STORAGE OVERVIEW ---'
  diskUsage "/data1"
  diskUsage "/data2"
  echo ''

  echo '--- ILLUMINA_DATA runs ---'
  find_subdirs_in_directory "/data1/illumina_data/"
  echo ''

  echo '--- PROCESSED runs ---'; 
  find_subdirs_in_directory "/data2/processed/"
  echo ''

  echo '--- PIPELINE CHECK files ---'
  find_file_in_pipeline_run "*logs/PipelineCheck.log"
  echo ''

  echo '--- WGS METRICS files ---'
  find_file_in_pipeline_run "*QCStats/WGSMetrics*transposed.txt"
  echo ''

  echo '--- QSTAT ---'
  qstat -u "*" | head -10
  totalJobCount=`qstat -u "*" | awk '$1 !~ /^job|---/' | wc -l`
  echo "[HMF_INFO] Total number of jobs in qstat: "$totalJobCount
  echo ""
}

find_subdirs_in_directory() { 
  find "$1" -mindepth 1 -maxdepth 1 -type d
}

find_file_in_pipeline_run() { 
  for run in /data2/processed/*; do 
    find ${run} -wholename "$1" && found=1
  done 
}

function diskUsage {
    local mount=$1
    local available=$( df -h "${mount}" | tail -1 | tr -s ' ' | cut -d" " -f 4 )
    local percString=$( ${cmdPrefix} df -h "${mount}" | tail -1 | tr -s ' ' | cut -d" " -f 5 )
    local percentage=$( echo ${percString} | sed 's/\%//g' )
    echo "${percString} used with ${available} space left for mount ${mount}"
}

main
