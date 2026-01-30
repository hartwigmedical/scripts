#!/usr/bin/env bash

TARGET_BUCKET="actin-shared-trial-hospital-curation-submission"

set -eo pipefail

[[ $# -ne 1 ]] && echo "Provide the path to the trial-curator request text file" && exit 1
[[ ! -f $1 ]] && echo "Cannot find specified input file [$1]" && exit 1

ns="$(date -u "+%N")"
targetFile="gs://${TARGET_BUCKET}/requested/$(date -u "+%Y%m%d-%H%M%S")$(echo $ns | cut -c1-3).txt"
gsutil cp "$1" "${targetFile}"
echo "Copied [$1] to [$targetFile]"
