#!/usr/bin/env bash

source message_functions || exit 1
source database_functions || exit 1

if [ -z "$1" ]; then
  echo "Usage: $0 <experiment_dir>"
  exit 1
fi

experiment_dir=$1

database="hmfpatients_pilot"
credentials=$(prod_writer_sql_credentials)
patient_db_jar=$(locate_pilot_patient_db)

processed_ids=()

for file in "${experiment_dir}"/*; do
  filename=$(basename "$file")
  patient_id=$(echo "$filename" | cut -d'.' -f1)

  if [[ " ${processed_ids[@]} " =~ " ${patient_id} " ]]; then
        echo "[INFO] Skipping already processed patient ${patient_id}"
      else
        echo "[INFO] Processing patient ${patient_id}"
        do_load_sigs_data "${patient_id}" "${sample_dir}" "${database}" "${credentials}" "${patient_db_jar}" "$@"
        if [ $? -ne 0 ]; then
          echo "Error loading sigs data for patient ${patient_id}"
          exit 1
        fi
        processed_ids+=("${patient_id}")
  fi
done

echo "All patients processed successfully."