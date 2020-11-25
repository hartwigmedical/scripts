#!/usr/bin/env bash

print_usage(){
    echo "-----"
    echo " Descr: Get relevant children of doids"
    echo " Usage: $(basename $0) -d <doids> "
    echo " Exmpl: $(basename $0) -d '1324,14534'"
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

all_doids=("")
for input_doid in $(echo ${doids} | sed "s/,/ /g")
do
    all_doids+=" ${input_doid} "
    url_doid="https://www.ebi.ac.uk/ols/api/ontologies/doid/terms/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FDOID_${input_doid}/children"
    json_doid=$( curl -s ${url_doid} )
    for row in $( echo ${json_doid} | jq -r '.[].terms[]? | .iri' )
    do
        row=${row#"http://purl.obolibrary.org/obo/DOID_"}
        all_doids+=" ${row} "
    done
done
all_doids=($(echo "${all_doids[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))


execute_sql_on_prod /data/common/repos/scripts/datarequest/sql/doids_database.sql > temp/doids_database.txt
echo ""

relevant_doids=("")
while IFS= read -r database_doids
do
    database_doids=$( echo ${database_doids} | sed 's/,/ /g' )
    database_doids=($(echo "${database_doids[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
    for doid in ${all_doids[@]}
    do
        for database_doid in ${database_doids[@]}
        do
            if [[ ${doid} == ${database_doid} ]]; then
                relevant_doids+=" ${doid}, "
            fi
        done
    done
done < temp/doids_database.txt
relevant_doids=${relevant_doids#" "}
relevant_doids=${relevant_doids%", "}
relevant_doids="WHERE doids IN (${relevant_doids})"
echo "SQL WHERE statement to be included for datarequest:"
echo ${relevant_doids}

rm -r temp