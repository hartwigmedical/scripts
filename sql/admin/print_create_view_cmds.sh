#!/usr/bin/env bash

source message_functions || exit 1

# Purpose: helper script for creating views (eg when doing database instance upgrades)
# Description: Prints mysqlsh commands to execute (so send stdout to a "create_views.sh" for instance)

cnf=$1 # options file https://dev.mysql.com/doc/refman/8.0/en/option-files.html
repos_root_path=$2 # all sql files with CREATE VIEW code is assumed to be in locally cloned repos

REPOS_ROOT="${repos_root_path:-/Users/stef/IdeaProjects}"
CONFIG_FILE=${cnf:-default-sql-config.cnf}

[[ -f "${CONFIG_FILE}" ]] || die "Cannot find config file [${CONFIG_FILE}]"

HMFPATIENTS=(
  "${REPOS_ROOT}/scripts/sql/views/create_view_first_treatment_after_biopsy.sql"
  "${REPOS_ROOT}/scripts/sql/views/create_view_first_matched_treatment_response.sql"
  "${REPOS_ROOT}/scripts/sql/views/create_view_hpc.sql"
  "${REPOS_ROOT}/scripts/sql/views/create_view_clinical.sql"
  "${REPOS_ROOT}/scripts/sql/views/create_view_clinical_purity.sql"
  "${REPOS_ROOT}/scripts/sql/views/create_view_datarequest_all.sql"
  "${REPOS_ROOT}/scripts/sql/views/create_view_datarequest.sql"
  "${REPOS_ROOT}/scripts/sql/views/create_view_datarequest_all_pseudonymized.sql"
  "${REPOS_ROOT}/scripts/sql/views/create_view_datarequest_with_nkiChecked.sql"
  "${REPOS_ROOT}/datarequest/views/dr_clinical_views.sql"
  "${REPOS_ROOT}/datarequest/views/summary.sql"
  "${REPOS_ROOT}/datarequest/views/datasets.sql"
)

CATALOG=(
  "${REPOS_ROOT}/catalog-bridge/views/create_producer_schema.sql"
  "${REPOS_ROOT}/catalog-bridge/views/create_consumer_schema.sql"
)

ACTIN_EMC_PHASE1=(
  "${REPOS_ROOT}/actin/database/src/main/resources/generate_views.sql"
)

ACTIN_PERSONALIZATION=(
  "${REPOS_ROOT}/actin-personalization/database/src/main/resources/generate_views.sql"
)

VICC_DB=(
  "${REPOS_ROOT}/scripts/sql/vicc/jax_view.sql"
  "${REPOS_ROOT}/scripts/sql/vicc/oncokb_biological_view.sql"
  "${REPOS_ROOT}/scripts/sql/vicc/oncokb_clinical_view.sql"
  "${REPOS_ROOT}/scripts/sql/vicc/vicc_drug_view.sql"
  "${REPOS_ROOT}/scripts/sql/vicc/vicc_feature_view.sql"
  "${REPOS_ROOT}/scripts/sql/vicc/vicc_view.sql"
)

ICLUSION=(
  "${REPOS_ROOT}/scripts/sql/iclusion/view/studyBlacklistedTumorLocations.sql"
  "${REPOS_ROOT}/scripts/sql/iclusion/view/studyMutationConditions.sql"
  "${REPOS_ROOT}/scripts/sql/iclusion/view/studyTumorLocations.sql"
)

CKB=(
  "${REPOS_ROOT}/scripts/sql/ckb/views/actionabilty_view.sql"
  "${REPOS_ROOT}/scripts/sql/ckb/views/drugs_view.sql"
  "${REPOS_ROOT}/scripts/sql/ckb/views/genes_view.sql"
  "${REPOS_ROOT}/scripts/sql/ckb/views/global_approval_view.sql"
  "${REPOS_ROOT}/scripts/sql/ckb/views/treatment_approach_view.sql"
  "${REPOS_ROOT}/scripts/sql/ckb/views/trials_countries_view.sql"
  "${REPOS_ROOT}/scripts/sql/ckb/views/trials_extended_view.sql"
  "${REPOS_ROOT}/scripts/sql/ckb/views/trials_view.sql"
  "${REPOS_ROOT}/scripts/sql/ckb/views/variants_view.sql"
)

ALL_SQL_FILES=(
    "${HMFPATIENTS[@]}" "${CATALOG[@]}" "${ACTIN_EMC_PHASE1[@]}" "${ACTIN_PERSONALIZATION[@]}"
    "${VICC_DB[@]}" "${ICLUSION[@]}" "${CKB[@]}"
)

main () {

    echo "# Checking all SQL files for existence"
    for sql_file in "${ALL_SQL_FILES[@]}"; do
        [[ -f "${sql_file}" ]] || die "Does not exist [${sql_file}]. Exiting";
    done

    echo "# Printing all all mysqlsh cmds"
    print_create_view_cmds "hmfpatients" "${HMFPATIENTS[@]}"
    print_create_view_cmds "catalog" "${CATALOG[@]}"
    print_create_view_cmds "actin_emc_phase1" "${ACTIN_EMC_PHASE1[@]}"
    print_create_view_cmds "actin_personalization" "${ACTIN_PERSONALIZATION[@]}"
    print_create_view_cmds "vicc_db" "${VICC_DB[@]}"
    print_create_view_cmds "vicc_db_pilot" "${VICC_DB[@]}"
    print_create_view_cmds "iclusion_latest" "${ICLUSION[@]}"
    print_create_view_cmds "iclusion_production" "${ICLUSION[@]}"
    print_create_view_cmds "ckb_latest" "${CKB[@]}"
}

function print_create_view_cmds () {
    local database_name=$1 && shift
    local sql_files=("$@")
    for sql_file in "${sql_files[@]}"; do
        echo "echo 'Processing SQL file [${sql_file}]'"
        echo "mysqlsh --defaults-extra-file=${CONFIG_FILE} --database ${database_name} -f ${sql_file}"
    done
}

main