#!/usr/bin/env bash

source message_functions || exit 1

# Purpose: helper script for deleting views (eg when doing database instance upgrades)
# Description: Prints DROP VIEW statements to be executed (either manually in db editor or write stdout to "drop_views.sql")

sql_config_file=$1 # options file https://dev.mysql.com/doc/refman/8.0/en/option-files.html

[[ -n "${sql_config_file}" ]] || die "Provide path to config file with host/user/pw"
[[ -f "${sql_config_file}" ]] || die "Does not exist [${sql_config_file}]"

echo "# Retrieving database names [SHOW databases]"
databases_json=$(mysqlsh --defaults-extra-file="${sql_config_file}" --log-level=error --result-format=json -e 'SHOW databases;')
readarray -t databases_names < <(jq -r '.[]' <<< "${databases_json}" | grep -vE '^(Database|information_schema|mysql|performance_schema|sys)')

echo "# Found the following databases"
for db in "${databases_names[@]}"; do
    echo "# Database/schema: ${db}"
done

echo "# Printing DROP statement for each view in each database"
for db in "${databases_names[@]}"; do
    views_json=$(
        mysqlsh --defaults-extra-file="${sql_config_file}" --log-level=error --result-format=json \
        -e "SELECT table_name FROM information_schema.views WHERE table_schema = '$db';"
    )
    readarray -t views_names < <(jq -r '.[]' <<< "${views_json}" | grep -vE '^table_name')
    for view in "${views_names[@]}"; do
        echo "DROP VIEW ${db}.${view};"
    done
done
