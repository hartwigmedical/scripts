#!/usr/bin/env bash

print_usage(){
    echo "-----"
    echo " Descr: Checks created GCP bucket for DR"
    echo " Usage: $(basename $0) -b <bucket_name> -e <gcp_mail>"
    echo " Exmpl: $(basename $0) -b 'hmf-dr-XXX' -e 'john@doe.com jaap@doe.com'"
    echo "-----"
    exit 1
}

while getopts ':b:e:' flag; do
    case "${flag}" in
        b) bucket_name=${OPTARG} ;;
        e) gcp_mail=${OPTARG} ;;
        *) print_usage
        exit 1 ;;
    esac
done

if [[ -z "${bucket_name}" || -z "${gcp_mail}" ]]; then
    print_usage
fi

echo ""
echo "[INFO] Files in the bucket" ${bucket_name}":"
gsutil du -h gs://${bucket_name}/
echo ""

echo "[INFO] Persmissions of the bucket"  ${bucket_name}":"
gsutil -u hmf-share iam get gs://${bucket_name}/
echo ""

for email in ${gcp_mail}
do
email_in_persmissions=$( gsutil -u hmf-share iam get gs://${bucket_name}/ | grep $email )
if [[ ${email_in_persmissions} == "" ]]; then
    echo "[ERROR] Account $email not added to IAM GCP bucket. Please add manually in the GUI (role = storage object viewer)."
else echo "[INFO] Permissions account $email set correctly."
fi
done
echo ""

