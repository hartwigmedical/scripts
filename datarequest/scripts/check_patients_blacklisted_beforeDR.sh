#!/usr/bin/env bash

#put SQL selectie query file path at $1 (for all use: /data/common/repos/scripts/datarequest/selection/all_patients_eligible_datarequests.sql)

execute_sql_on_prod $1 | sort | uniq > temp_datarequest.tsv
awk -F"\t" 'FNR > 1 {print $1}' temp_datarequest.tsv | sort | uniq > testfile.tmp && mv testfile.tmp temp_datarequest.tsv

wc -l /data/lims/patient_blacklist.tsv
wc -l temp_datarequest.tsv
echo 'Overlap between two files:'
comm -1 -2 <(sort /data/lims/patient_blacklist.tsv) <(sort temp_datarequest.tsv) | wc -l

rm temp_datarequest.tsv