#!/bin/bash

## start
echo ''
echo '======== STATUS ======='
echo ''

## show storage status
echo '--- STORAGE OVERVIEW ---'
df -h | grep -P '\/data\d+'
echo ''

## show content of data1/2 dirs
echo '--- ILLUMINA_DATA runs ---'; 
ls -ld /data1/illumina_data/*; 
echo ''

echo '--- PROCESSED runs ---'; 
ls -ld /data2/processed/*
echo ''

echo '--- PIPELINE CHECK files ---'
ls -lh /data2/processed/*/logs/PipelineCheck.log
echo ''

echo '--- WGS METRICS files ---'
ls -lh /data2/processed/*/QCStats/WGSMetrics*transposed.txt
echo ''

## show cluster activity
echo '--- QSTAT ---'
qstat -u "*" | head -15
totalJobCount=`qstat -u "*" | awk '$1 !~ /^job|---/' | wc -l`
echo "[HMF_INFO] Total number of jobs in qstat: "$totalJobCount
echo ""

