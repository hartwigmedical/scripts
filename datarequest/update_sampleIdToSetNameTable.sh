#!/usr/bin/env bash
source /data/common/repos/scripts/hmftools/metadata/load_metadata

runsDir="/data/cpct/runs"
outFile="/data/data_archive/datarequests/procedure/sampleid2setname.tsv"
errFile="/data/data_archive/datarequests/procedure/sampleid2setname.err"

echo "[INFO] Reading ${runsDir}"
echo "[INFO]   writing output table to ${outFile}"
echo "[INFO]   writing error log to ${errFile}"

for setPath in $( find ${runsDir} -mindepth 1 -maxdepth 1 -type d ); do 
  if [[ -f ${setPath}/metadata || -f ${setPath}/metadata.json ]]; then
    setName=$( basename ${setPath} )
    sampleId=$( load_tumor_sample_from_metadata ${setPath} )
    echo -e "${sampleId}\d${setName}"
  fi
done > ${outFile} 2>${errFile}
