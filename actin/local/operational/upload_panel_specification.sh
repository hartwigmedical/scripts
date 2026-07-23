#!/usr/bin/env bash

[[ $# -ne 1 ]] && echo "Provide the TSV containing the panel specification you want to upload" && exit 1

echo "Ensure the tunnel has been started to the resource API in another terminal"
curl --data-binary @$1 http://localhost:8088/panel_specifications -H "Content-Type: text/tsv"
