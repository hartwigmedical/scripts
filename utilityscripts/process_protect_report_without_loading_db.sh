#!/usr/bin/env bash

run_dir=$1 && shift

echo "[INFO] Re-run PROTECT with updated tumor location of run ${run_dir}"
run_protect_prod ${run_dir}

echo "[INFO] Create patient report of ${run_dir}"
create_patient_report_for_run ${run_dir}