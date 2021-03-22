#!/usr/bin/env bash

print_usage(){
    echo "-----"
    echo " Descr: Get risks iso27001 from BizzMine"
    echo " Usage: $(basename $0)"
    echo " Exmpl: $(basename $0)"
    echo "-----"
    exit 1
}

#while getopts ':i:' flag; do
#    case "${flag}" in
#        i) dr_id=${OPTARG} ;;
#        *) print_usage
#        exit 1 ;;
#    esac
#done
#
#if [[ -z "${dr_id}" ]]; then
#    echo "[ERROR] script get_risks_iso27001_BizzMine did not run, check usage below:"
#    print_usage
#fi


api_url=$"https://api.bizzmine.cloud/collection/RiskAnalyses/"
api_token=$( cat /data/common/dbs/api_credentials/bizzmine/api_token )

## start with script
echo ""
echo "[START] get_risks_iso27001_BizzMine: $(date +"%y%m%d (%T)")"
echo ""

curl get -s -H "x-token: ${api_token}" -H "x-tenant: hartwigmedicalfoundation" ${api_url}/instances > temp.json
n_risks=$( jq -r '.[] | select(.RiskAnalyses_Regulation.text=="ISO27001") | .RiskAnalyses_RiskAnalysesID ' temp.json | sort |  uniq | wc -l )
echo ""
echo "--- total number of ISO27001 risks ---"
echo ${n_risks}

jq -r '.[] | select(.RiskAnalyses_Regulation.text=="ISO27001") | .ISO27001Annexes[]?.ISO27001Annexes_Annex' temp.json > annexes.txt
declare -a annexes=("A.5 Information security policies" "A.6 Organization of information security" "A.7 Human resource security" "A.8 Asset management" "A.9 Access control" "A.10 Cryptography" "A.11 Physical and environmental security" "A.12 Operations security" "A.13 Communications security" "A.14 System acquisition, development and maintenance" "A.15 Supplier relationships" "A.16 Information security incident management" "A.17 Information security aspects of business continuity management" "A.18 Compliance" )
for annex in "${annexes[@]}"; do
  echo ""
  echo ${annex}
  annex=$( echo ${annex:0:4} | sed -e 's/^[[:space:]]*//' )
  n_head_annex=$( cat annexes.txt | grep "${annex}" | wc -l )
  echo -n "- number of risks under head annex: "; echo ${n_head_annex}
  if [ $( cat annexes.txt | grep "${annex}" | wc -l ) != 0 ]; then
      echo "- divided over (sub)annexes"
      cat annexes.txt | grep "${annex}" | sort | uniq > sub_annexes.txt
      while read sub_annex; do
         echo -n "     "${sub_annex}; echo ": "$( cat annexes.txt | grep "${sub_annex}" | wc -l)
      done <sub_annexes.txt
  fi
done
echo ""

rm temp.json; rm annexes.txt; rm sub_annexes.txt
