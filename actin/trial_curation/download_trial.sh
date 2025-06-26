#!/bin/bash

set -e

download_directory="YOUR_DOWNLOAD_FOLDER"

download_script="${PATH_TO_TRIAL_CURATOR_REPO}/trialcurator/download_trial.py"

nct_ids=(
NCT04626635
NCT06644755
)

source /opt/anaconda3/etc/profile.d/conda.sh
conda activate trial_curator

for trial_id in "${nct_ids[@]}"; do
    
    echo "Downloading trial: $trial_id"

    python "$download_script" \
        "$trial_id" \
        -o "$download_directory/${trial_id}.json"

done
