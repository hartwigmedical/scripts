#!/bin/bash

set -e

input_directory="YOUR_INPUT_DOWNLOAD_FOLDER"
output_directory="YOUR_OUTPUT_DOWNLOAD_FOLDER"
trial_curator_repo_directory="PATH_TO_TRIAL_CURATOR_REPO"

curation_script="$trial_curator_repo_directory/actin_curator/actin_eligibility_curator.py"


# LLM_provider="Google"
LLM_provider="OpenAI"

nct_ids=(
    'NCT04626635'
)

source /opt/anaconda3/etc/profile.d/conda.sh
conda activate trial_curator

export OPENAI_API_KEY="YOUR_API_KEY"

export GEMINI_PROJECT_ID="hmf-aus"
export GEMINI_LOCATION="europe-west4"


for trial_id in "${nct_ids[@]}"; do

	echo "=== PROCESSING ${trial_id} ... ==="

	trial_log_file="${output_directory}/log_${trial_id}_$(date +%Y%m%d_%H%M%S).log"

	export PYTHONPATH=$trial_curator_repo_directory

	python "$curation_script" \
	--LLM_provider $LLM_provider \
	--trial_json "$input_directory/${trial_id}.json" \
	--out_trial_file "$output_directory/${trial_id}.json" \
	--ACTIN_path "$trial_curator_repo_directory/actin_curator/data/ACTIN_rules/ACTIN_rules_w_categories_13062025.csv" \
	--log_level DEBUG \
	2>&1 | tee -a "$trial_log_file"

done
