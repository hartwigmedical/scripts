#!/usr/bin/env bash

print_usage(){
    echo "-----"
    echo " Descr: Get relevant children of doids"
    echo " Usage: $(basename $0) -d <doids> "
    echo " Exmpl: $(basename $0) -d '1234,14534'"
    echo "-----"
    exit 1
}

while getopts ':d:' flag; do
    case "${flag}" in
        d) doids=${OPTARG} ;;
        *) print_usage
        exit 1 ;;
    esac
done

if [[  -z "${doids}" ]]; then
    echo "[ERROR] script get_relevant_children_of_doids did not run, check usage below:"
    print_usage
fi


echo ""
echo "[START] get_relevant_children_of_doids: $(date +"%y%m%d (%T)")"
echo ""

mkdir temp

execute_sql_on_prod /data/common/repos/scripts/datarequest/sql/doids_database.sql > temp/doids_database.txt



