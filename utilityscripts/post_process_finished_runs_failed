#!/usr/bin/env bash

source lims_functions
source metadata_functions

health_check_processed_runs

echo "[INFO] Count HealthChecker succeeded"
health_check_processed_runs | grep -c "HealthChecker has succeeded for"

purple_version=2.53
echo "[INFO] Count analysis purple version v${purple_version}"
health_check_processed_runs | grep -c "This analysis has been done with purple version=${purple_version}"

source_path=/data/gcp/processed_runs

process_runs=$(ls ${source_path})
for run in ${process_runs}; do
    sample_id=$(load_tumor_sample_from_metadata ${source_path}/${run})
    tumor_barcode=$(find_barcode_for_sample_name ${sample_id})
    cohort=$(find_cohort_type ${tumor_barcode})

    # do some jobs for all runs
    echo "[INFO] Running PROTECT"
    run_protect_prod ${source_path}/${run}

    if [[ -z ${cohort} ]]; then
        echo "[ERROR] Could not resolve cohort for ${run}."
    elif [[ ${cohort} == "CORELR02" || ${cohort} == "CORERI02" ]]; then
        echo "[INFO] Moving set ${run} to /data/core/runs"
        mv ${source_path}/${run} /data/core/runs/
    elif [[ ${cohort} == "CPCT" || ${cohort} == "CPCTpancreas" || ${cohort} == "DRUP" || ${cohort} == "DRUPstage3" || ${cohort} == "COREDB" ]]; then
        echo "[INFO] Moving set ${run} to /data/cpct/runs"
        mv ${source_path}/${run} /data/cpct/runs/
        echo "[INFO] Copy run to hmf-crunch at GCP"
        copy_run_to_gcp_crunch_project ${run}
        echo "[INFO] Loading run into database"
        load_run_into_prod_db /data/cpct/runs/${run}
    else
        # This is for patients that require a summary (WIDE, some CORE)
        echo "[INFO] Copying set ${run} to /data/cpct/reportable_runs"
        cp -R ${source_path}/${run} /data/cpct/reportable_runs/
        if [[ ${cohort} == "CORE" || ${cohort} == "CORELR11" || ${cohort} == "CORESC11" ]]; then
            echo "[INFO] Moving set ${run} to /data/core/runs"
            mv ${source_path}/${run} /data/core/runs/
        elif [[ ${cohort} == "WIDE" ]]; then
            echo "[INFO] Moving set ${run} to /data/cpct/runs"
            mv ${source_path}/${run} /data/cpct/runs/
            echo "[INFO] Copy run to hmf-crunch at GCP"
            copy_run_to_gcp_crunch_project ${run}
            echo "[INFO] Loading run into database"
            load_run_into_prod_db /data/cpct/runs/${run}
        fi
    fi
done
