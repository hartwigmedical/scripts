#!/usr/bin/env bash

source locate_files || exit 1
source message_functions || exit 1

sql_file_or_select_statement=$1 && shift

database_name="hmfpatients"
credentials=$(locate_prod_backup_database_credentials_cnf)

do_execute_sql_on_database "${sql_file_or_select_statement}" "${database_name}" "${credentials}"
