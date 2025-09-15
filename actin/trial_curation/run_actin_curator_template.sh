#!/bin/bash
set -e

# reminder 1: `chmod +x ./run_actin_curator.sh` before using this script for the first time
# reminder 2: Pull the latest changes from whatever git branch you are using
# Once the following set-up is complete (the list of TODOs), Trial ID is the only argument required to run the ACTIN curator
# e.g. ./run_actin_curator.sh NCT06818643


trial_id=$1

if [ -z $trial_id ]; then
	echo "Trial ID is empty."
	exit 1
else
	echo "Trial ID $trial_id entered."
fi


# TODO: Specify where you have cloned https://github.com/hartwigmedical/trial-curator
repo_dir="..."

if [ ! -d "$repo_dir" ]; then
	echo "Specify the trial curator repository path in the script."
	exit 1
fi

cd $repo_dir
export PYTHONPATH="$repo_dir"


# TODO: In production, this would be always 'master'. We are currently testing version 2 ('version2_experimental')
branch_ver="version2_experimental"

git checkout $branch_ver
echo "You are on branch: $branch_ver"


# TODO: In time the curator will use other LLMs. Right now it's designed specifically to OpenAI's GPT-4o (version: gpt-4o-2024-08-06)
llm_provider="OpenAI"

echo "Using LLM provider: $llm_provider"


# TODO: Copy your OpenAI API key below
export OPENAI_API_KEY="..."


# TODO: Specify the Conda environment created for trial curation
conda_env="trial_curator"

source /opt/miniconda3/etc/profile.d/conda.sh
conda activate $conda_env
echo "Activated conda env: $conda_env"


# TODO: Change the input and output folders as needed
input_dir="..."
output_dir="..."

mkdir -p "$input_dir" "$output_dir"


download_app="${repo_dir}/trialcurator/download_trial.py"
if [ ! -f "${input_dir}/${trial_id}.json" ]; then
	echo -e "Downloading ${trial_id}.\n"
else
	echo -e "${trial_id} is already available.\n"
fi

python "$download_app" \
	"$trial_id" \
	--output "${input_dir}/${trial_id}.json"


# TODO: This should be one-time setup, unless new rules are added
actin_filepath="${repo_dir}/actin_curator/data/ACTIN_rules/ACTIN_rules_w_categories_13062025.csv"

actin_curation_engine="actin_curator.actin_curator"
log="${output_dir}/log_${trial_id}_$(date +%Y%m%d_%H%M%S).log"

python -m "$actin_curation_engine" \
	--llm_provider $llm_provider \
	--input_file "$input_dir/${trial_id}.json" \
	--output_file_complete "$output_dir/${trial_id}.json" \
	--output_file_concise "$output_dir/${trial_id}_concise.tsv" \
	--actin_filepath $actin_filepath \
	--log_level DEBUG \
	2>&1 | tee -a "$log"

