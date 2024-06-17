#!/usr/bin/python

import sys
import subprocess

batch = sys.argv[1]

listRelevantBuckets = 'gsutil ls gs://targeted-pipeline-output-prod-1/ | grep '+batch

cmd ="for i in  `"+listRelevantBuckets+ "`; do echo -n \"$i \"; gsutil cat ${i}metadata.json | jq '.tumor.sampleName'   ;done"

capture_output=subprocess.run(cmd,shell=True,capture_output=True,text=True)

res=capture_output.stdout
#containts space seperated list of buckets and sampleIds in this batch

res=res.rstrip('\n')


lines_arr = res.split('\n')



body="""workflow: "panel-plot"
version: "1.0.1"
params:\n"""

outputYaml=""

cnt=0

for lines in lines_arr:
  
  arr = lines.split(' ')
  inputBucket = arr[0]
  sampleId = arr[1]

  outputYaml += 'name: \"'  + batch+"_"+str(cnt) +'\"\n'
  outputYaml += body
  outputYaml+=  '  run_uri: ' + '\"' + inputBucket +'"\n'
  outputYaml += '  sample_id: '  + sampleId +'\n'
  outputYaml += "---\n"
  cnt+=1
outputYaml = outputYaml.rstrip('\n')
outputYaml = outputYaml.rstrip('-')
print(outputYaml)
