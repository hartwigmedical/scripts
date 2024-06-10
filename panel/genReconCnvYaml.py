#!/usr/bin/python

import sys
import subprocess

batch = sys.argv[1]

listRelevantBuckets = 'gsutil ls gs://targeted-pipeline-output-prod-1/ | grep '+batch

cmd ="for i in  `"+listRelevantBuckets+ "`; do echo -n \"$i \"; gsutil cat ${i}metadata.json | jq '.tumor.sampleName'   ;done"

capture_output=subprocess.run(cmd,shell=True,capture_output=True,text=True)

res=capture_output.stdout
res=res.rstrip('\n')
arr = res.split('\n')


body="""workflow: "panel-plot"
version: "1.0.1"
params:\n"""

outputYaml=""

cnt=0

for i in arr:
  
  arr2 = i.split(' ')
  outputYaml += 'name: \"'  + batch+"_"+str(cnt) +'\"\n'
  outputYaml += body
  outputYaml+=  '  run_uri: ' + '\"' + arr2[0] +'"\n'
  outputYaml += '  sample_id: '  + arr2[1] +'\n'
  outputYaml += "---\n"
  cnt+=1
outputYaml = outputYaml.rstrip('\n')
outputYaml = outputYaml.rstrip('-')
print(outputYaml)
