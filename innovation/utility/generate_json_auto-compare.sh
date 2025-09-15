#!/usr/bin/env bash
# generate_sample_map.sh
# This script takes a TSV file with two columns (truth_id <TAB> target_id)
# and outputs a JSON mapping where each truth sample ID is a key with a list of target sample IDs as its value.
# Usage: ./generate_sample_map.sh input.tsv

set -euo pipefail

# Check if the input file was provided.
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 input.tsv"
    exit 1
fi

tsv_file="$1"

# Check if the file exists.
if [ ! -f "$tsv_file" ]; then
    echo "Error: File '$tsv_file' not found!"
    exit 1
fi

# Declare an associative array to store the mapping.
declare -A sample_map

# Read each line from the TSV file.
while IFS=$'\t' read -r truth target || [ -n "$truth" ]; do
    # Skip empty lines.
    if [[ -z "$truth" ]] || [[ -z "$target" ]]; then
        continue
    fi
    # Append the target sample ID to the list for this truth sample ID.
    if [[ -v sample_map["$truth"] ]]; then
        sample_map["$truth"]+=",${target}"
    else
        sample_map["$truth"]="$target"
    fi
done < "$tsv_file"

# Build the JSON output.
json="{"
first_entry=true
for truth in "${!sample_map[@]}"; do
    if [ "$first_entry" = true ]; then
        first_entry=false
    else
        json+=", "
    fi
    # Split the comma-separated list of target sample IDs into a JSON array.
    IFS=',' read -ra targets <<< "${sample_map[$truth]}"
    json_array="["
    first_target=true
    for target in "${targets[@]}"; do
        # Trim any surrounding whitespace.
        target="$(echo "$target" | sed 's/^ *//; s/ *$//')"
        if [ "$first_target" = true ]; then
            first_target=false
        else
            json_array+=", "
        fi
        json_array+="\"$target\""
    done
    json_array+="]"
    
    json+="\"$truth\": $json_array"
done
json+="}"

# Print the resulting JSON mapping.
echo "$json"