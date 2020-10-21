#!/usr/bin/env bash

print_usage(){
    echo "-----"
    echo " Descr: Add researcher to data on GCP for a specific DR"
    echo " Usage: $(basename $0) -i <dr-id> -s <suffix> -e <email>"
    echo " Exmpl: $(basename $0) -i 'DR-001' -s 'update1' -e 'john@doe.com,jaap@doe.com'"
    echo "-----"
    exit 1
}

while getopts ':i:s:e:' flag; do
    case "${flag}" in
        i) dr_id=${OPTARG} ;;
        s) dr_suffix=${OPTARG} ;;
        e) gcp_mail=${OPTARG} ;;
        *) print_usage
        exit 1 ;;
    esac
done

if [[ -z "${dr_id}" || -z "${gcp_mail}" ]]; then
    echo "[ERROR] script add_researcher_to_data_gcp did not run, check usage below:"
    print_usage
fi


api_url=$"https://api.hartwigmedicalfoundation.nl"
api_url_spec=$"${api_url}/hmf/v1"
api_key=$"/data/common/dbs/api_credentials/api.key"
api_cert=$"/data/common/dbs/api_credentials/api.crt"


## quick input checks
[[ ! -z "${dr_id}" && "${dr_id}" =~ ^DR ]] || die "dr-id incorrect (${dr_id})?"
[[ ! -z "${gcp_mail}" && "${gcp_mail}" =~ \@.+\. ]] || die "gcp_mail incorrect (${gcp_mail})"

## we need the index of DR
release_id=$"${dr_id}"
dr_index=$( echo "${dr_id}" | sed 's/^DR\-//')
request_id=$"${dr_index}"

## reset release/request ids to include suffix if given
if [[ "${dr_suffix}" != "" ]]; then
    release_id="${dr_id}-${dr_suffix}"
    request_id="${dr_index}-${dr_suffix}"
fi

bucket_name=$"hmf-dr-${request_id}"

## start with script
echo ""
echo "[START] add_researcher_to_data_gcp: $(date +"%y%m%d (%T)")"
echo ""

echo "[INFO] Accounts inputted in the script for which access to GCP data will tried to be given:"
echo ${gcp_mail}
echo ""

check_whether_GCPaccount_exists -e ${gcp_mail}

echo ""
echo "[INFO] ADD ${gcp_mail} TO IAM OF BUCKET: ${bucket_name} RELATED TO: ${release_id}."
echo ""

echo "[INFO] Current permissions of the bucket ${bucket_name}:"
gsutil -u hmf-share iam get gs://${bucket_name}/
echo ""

for email in $(echo ${gcp_mail} | sed "s/,/ /g")
do
    email_in_permissions=$( gsutil -u hmf-share iam get gs://${bucket_name}/ | grep $email )
    if [[ ${email_in_permissions} == "" ]]; then
        echo "....Adding account $email to the iam of the bucket...."
        gsutil -u hmf-share iam ch user:$email:objectViewer gs://${bucket_name}/
        email_in_permissions=$( gsutil -u hmf-share iam get gs://${bucket_name}/ | grep $email )
        if [[ ${email_in_permissions} == "" ]]; then
            echo "[ERROR] Account $email not added as user in the bucket."
            echo ""
        else
            echo "[INFO] Account $email correctly added as user:viewer to the gcp bucket."
            echo "[INFO] Updated permissions of the bucket ${bucket_name}:"
            gsutil -u hmf-share iam get gs://${bucket_name}/
            echo ""
        fi
        echo ""
    else
        echo "[ERROR] Account $email already present in current permissions of the bucket, so adding to iam not needed."
        echo ""
     fi
done
echo ""


group_id=$( curl -s --cert ${api_cert} --key ${api_key} ${api_url_spec}/groups | jq --arg release_id "$release_id" '.[] | select(.name==$release_id) | .id')
if [[ "${group_id}" == "" ]]; then
    echo "[INFO] NO API GROUP (for ACL permissions) RELATED TO: ${release_id}."
else
    echo ""
    echo "[INFO] ADD $gcp_mail to API GROUP (for ACL permissions) RELATED TO: ${release_id}."
    echo ""
    echo "[INFO] Current accounts in API group related to ${release_id}:"
    account_nrs=$( curl -s --cert ${api_cert} --key ${api_key} ${api_url_spec}/groups/${group_id}/members | jq '.[] | .account_id' )
    email_in_acl_group=$""
    for account_nr in $account_nrs
    do
        curl -s --cert ${api_cert} --key ${api_key} ${api_url_spec}/accounts/${account_nr} | jq '.email' | sed -e 's/^"//' -e 's/"$//'
        email_acl=$( curl -s --cert ${api_cert} --key ${api_key} ${api_url_spec}/accounts/${account_nr} | jq '.email' | sed -e 's/^"//' -e 's/"$//' )
        email_in_acl_group+=${email_acl}
    done
    echo ""

    for email in $(echo ${gcp_mail} | sed "s/,/ /g")
    do
        if [[ $( echo ${email_in_acl_group} | grep $email ) == "" ]]; then
            echo "....Adding $email to API group  ...."
            if [[ $( curl -s --cert ${api_cert} --key ${api_key} ${api_url_spec}/accounts | grep $email ) == "" ]]; then
                account_nr_new=$( curl -s --cert ${api_cert} --key ${api_key} -d '{"email": "$email"}' -H "Content-Type: application/json" -X POST ${api_url_spec}/accounts | jq .[] )
            else
                account_nr_new=$( curl -s --cert ${api_cert} --key ${api_key} ${api_url_spec}/accounts | jq --arg email_select "$email" '.[] | select(.email==$email_select) | .id' )
            fi
            curl -s --cert ${api_cert} --key ${api_key} -d '{"group_id": "$group_id", "account_id": $account_nr_new}' -H "Content-Type: application/json" -X POST ${api_url_spec}/groups/${group_id}/members
        else
            echo "[ERROR] Account $email already present in API group related to ${release_id}, so no adding needed."
        fi
    done
    echo ""

    echo "[INFO] Updated accounts in API group related to ${release_id}:"
    account_nrs=$( curl -s --cert ${api_cert} --key ${api_key} ${api_url_spec}/groups/${group_id}/members | jq '.[] | .account_id' )
    for account_nr in $account_nrs
    do
        curl -s --cert ${api_cert} --key ${api_key} ${api_url_spec}/accounts/${account_nr} | jq '.email' | sed -e 's/^"//' -e 's/"$//'
    done
    echo ""

    echo "[INFO] Current number of raw files with API group related to ${release_id} in the ACL:"
    curl -s --cert ${api_cert} --key ${api_key} ${api_url_spec}/groups/${group_id}/files | jq . | grep file_id | wc -l
    echo ""
fi

echo ""
echo ""

