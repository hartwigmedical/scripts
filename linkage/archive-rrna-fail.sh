#!/bin/bash

# Copies the FASTQ files of an RNA sample from the production bucket to gs://hartwig-reads/archive/<sample-id>
# when the RNA pipeline failed because the sample contains too much ribosomal RNA.

set -eo pipefail

if [ $# -ne 1 ]; then
  echo "Usage: $0 <h-number>"
  exit 1
fi

hartwig_number=$1

readme=$(mktemp)
trap 'rm -f $readme' EXIT

sample_id=$(map_sample_ids -t research -- "$hartwig_number")


sample_id_lower=$(echo "$sample_id" | tr '[:upper:]' '[:lower:]')
hartwig_number_lower=$(echo "$hartwig_number" | tr '[:upper:]' '[:lower:]')

error_message=$(gcloud storage cat "gs://rna-pipeline-output-prod-1/runtime/$hartwig_number_lower/rrna-check.log" | sed -n '4,4p;5q' | sed -e 's/.*ERROR: //')
if [ -z "$error_message" ]; then
  echo "Sample $sample_id/$hartwig_number did not fail rRNA check?"
  gcloud storage cat "gs://rna-pipeline-output-prod-1/runtime/$hartwig_number_lower/rrna-check.log"
  exit 1
fi

echo "# $sample_id" > "$readme"
echo "" >> "$readme"
echo "$error_message" >> "$readme"

cat "$readme"

api_sample_id=$(api -j samples "name=$sample_id&type=tumor-rna" | jq '.[0].id')
fastqs=$(api -j fastq "$api_sample_id" sample_id | jq -r '.[] | "gs://\(.bucket)/\(.name_r1)", "gs://\(.bucket)/\(.name_r2)"')

gcloud storage cp --billing-project hmf-crunch "$readme" "gs://hartwig-reads/archive/$sample_id_lower/README.md"

for fastq in $fastqs; do
  archive_path=${fastq//gs:\/\//}
  gcloud storage cp --billing-project hmf-crunch "$fastq" "gs://hartwig-reads/archive/$sample_id_lower/$archive_path"
done
