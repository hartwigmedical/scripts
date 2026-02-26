#!/usr/bin/env python3
import sys

input_bucket_uri = "gs://targeted-pipeline-output-prod-1/"

batch = sys.argv[1]

if not batch:
    print("Usage: genRemarksYaml.py <batch>")
    exit(1)

body="""workflow: "oncoact-panel-remarks"
version: "1.0.1"
params:\n"""

outputYaml=""

outputYaml += 'name: "remarks-' + batch +'"\n'
outputYaml += body
outputYaml+=  '  input_bucket_uri: "' + input_bucket_uri +'"\n'
outputYaml += '  run_names: "' + batch +'*"\n'

print(outputYaml)
