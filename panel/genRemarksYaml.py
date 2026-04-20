#!/usr/bin/env python3
import sys

input_bucket_uri = "gs://targeted-pipeline-output-prod-1"

run_id = sys.argv[1]

if not run_id:
    print("Usage: genRemarksYaml.py <run_id>")
    exit(1)

body="""workflow: "oncoact-panel-remarks"
version: "1.1.0"
params:\n"""

outputYaml=""

outputYaml += 'name: "' + run_id.replace('_Panel', '', 1) +'"\n'
outputYaml += body
outputYaml+=  '  input_bucket_uri: "' + input_bucket_uri + "/" + run_id +'"\n'

print(outputYaml)
