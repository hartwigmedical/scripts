#!/usr/bin/env bash

source message_functions || exit 1

print_usage() {
    echo "-----"
    echo " Descr: Get risks iso27001 from BizzMine"
    echo " Usage: $(basename $0) -i <indication_to_show_seperate_risk_0_1>"
    echo " Exmpl: $(basename $0) -i 1"
    echo "-----"
    exit 1
}

while getopts ':i:' flag; do
    case "${flag}" in
        i) ind=${OPTARG} ;;
        *) print_usage
        exit 1 ;;
    esac
done

if [[ -z "${ind}" ]]; then
    warn "script get_risks_iso27001_BizzMine did not run, check usage below:"
    print_usage
fi

api_url=$"https://api.bizzmine.cloud/collection/RiskAnalyses/"
api_token=$( gcloud secrets versions access "latest" --secret=bizzmine-api-token-dr --project=hmf-secrets )

## start with script
echo ""
echo "[START] get_risks_iso27001_BizzMine: $(date +"%y%m%d (%T)")"
echo ""

curl get -s -H "x-token: ${api_token}" -H "x-tenant: hartwigmedicalfoundation" ${api_url}/instances > temp.json
n_risks=$( jq -r '.[] | select(.RiskAnalyses_Regulation.text=="ISO27001") | .RiskAnalyses_RiskAnalysesID ' temp.json | sort |  uniq | wc -l )
echo ""
echo "--- total number of ISO27001 risks ---"
n_risks=$(expr ${n_risks} )
echo ${n_risks}

jq -r '.[] | select(.RiskAnalyses_Regulation.text=="ISO27001") | .ISO27001Annexes[]?.ISO27001Annexes_Annex' temp.json > annexes.txt
jq -r '.[] | select(.RiskAnalyses_Regulation.text=="ISO27001") | .RiskAnalyses_Name + ","+ .ISO27001Annexes[]?.ISO27001Annexes_Annex' temp.json > annexes_risks.txt

declare -a annexes=("A.5 Information security policies" "A.6 Organization of information security" "A.7 Human resource security" "A.8 Asset management" "A.9 Access control" "A.10 Cryptography" "A.11 Physical and environmental security" "A.12 Operations security" "A.13 Communications security" "A.14 System acquisition, development and maintenance" "A.15 Supplier relationships" "A.16 Information security incident management" "A.17 Information security aspects of business continuity management" "A.18 Compliance" )
for annex in "${annexes[@]}"; do
    echo ""
    echo ${annex}
    annex=$( echo ${annex:0:4} | sed -e 's/^[[:space:]]*//' )
    n_head_annex=$( cat annexes.txt | grep "${annex}" | wc -l )
    n_head_annex=$(expr ${n_head_annex} )
    perc=$( awk 'BEGIN{printf "%.0f\n",'${n_head_annex}'/'${n_risks}'*100}' )
    echo -n "- number of risks under head annex: "; echo -n ${n_head_annex}; echo " ("${perc}"%)"
    if [[ ${ind} == 1 ]]; then
        cat annexes_risks.txt | grep "${annex}" | csvcut -c 1 | sort | uniq | sed 's/^/    /'
    fi
    if [[ $(cat annexes.txt | grep -c "${annex}") != 0 ]]; then
        echo "- divided over (sub)annexes"
        cat annexes.txt | grep "${annex}" | sort | uniq > sub_annexes.txt
        while read sub_annex; do
            n_sub_annex=$( cat annexes.txt | grep "${sub_annex}" | wc -l )
            n_sub_annex=$(expr ${n_sub_annex} )
            perc=$( awk 'BEGIN{printf "%.0f\n",'${n_sub_annex}'/'${n_risks}'*100}' )
            echo -n "     "${sub_annex}": "; echo -n ${n_sub_annex}; echo " ("${perc}"%)"
            if [[ ${ind} == 1 ]]; then
                cat annexes_risks.txt | grep "${sub_annex}" | csvcut -c 1 | sort | uniq | sed 's/^/       /'
            fi
        done <sub_annexes.txt
    fi
done
echo ""

rm temp.json; rm annexes.txt; rm annexes_risks.txt; rm sub_annexes.txt
