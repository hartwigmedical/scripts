#!/usr/bin/env python3
import json
import sys
import subprocess

def run_bash_command(command: str) -> str:
    return subprocess.run(command,shell=True,capture_output=True,text=True).stdout

batch = sys.argv[1]

body="""workflow: "sage-visualisation"
version: "0.1.6"
params:\n"""

relevant_buckets = run_bash_command('gsutil ls gs://targeted-pipeline-output-prod-1/ | grep ' + batch).split("\n")

outputYaml=""
for input_bucket in relevant_buckets:
    if input_bucket:
        metadata_json_text = run_bash_command("gsutil cat " + input_bucket + "metadata.json")
        metadata_json = json.loads(metadata_json_text)
        sample_id = metadata_json["tumor"]["sampleName"]
        barcode = metadata_json["tumor"]["barcode"]

        outputYaml += 'name: "' + barcode +'"\n'
        outputYaml += body
        outputYaml+=  '  run_uri: "' + input_bucket +'"\n'
        outputYaml += '  sample_id: "' + sample_id +'"\n'
        outputYaml += "---\n"

outputYaml = outputYaml.rstrip('\n')
outputYaml = outputYaml.rstrip('-')
print(outputYaml)
