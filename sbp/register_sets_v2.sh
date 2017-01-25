#!/bin/bash

#######
# Hmf Register script for the HMF API
# Written for HMF by A. Repton
# Requires key and cert that allows access through OpenAM
#
# Exit Codes:
#  2: Missing Files
#  3: Unknown argument or help
#  4: Missing Arguments
######


CheckFileExists(){
  if [ ! -f $1 ]; then
    echo "$1 is missing. Please put in the same directory as this script"
    exit 2
  fi
}

PrintHelp(){
  echo "Please use with -e EntityName -i IniName -s SetName -r ReferenceSampleName -t TumourSampleName"
  exit 3
}

ENTITY_NAME=''
INI_NAME='' 
SET_NAME='' 
REF_SAMPLE='' 
TUMOR_SAMPLE=''
CRT_DIR=/home/sbp/hmfupload/

while getopts 'e:i:s:r:t:h' flag; do
  case "${flag}" in
    e) ENTITY_NAME=${OPTARG} ;;
    i) INI_NAME="${OPTARG}" ;;
    s) SET_NAME="${OPTARG}" ;;
    t) TUMOR_SAMPLE="${OPTARG}" ;;
    r) REF_SAMPLE="${OPTARG}" ;;
    h) PrintHelp ;;
    *) error "Unknown option ${OPTARG}. Please use -h for usage" ;;
  esac
done

if [[ -z $ENTITY_NAME ]] || [[ -z $INI_NAME ]] || [[ -z $SET_NAME ]] || [[ -z $REF_SAMPLE ]] || [[ -z $TUMOR_SAMPLE ]]; then
  echo "You are missing one or more options"
  PrintHelp
fi

for FILE in $CRT_DIR/hmf.crt $CRT_DIR/hmf.key; do
  CheckFileExists $FILE
done

/usr/bin/curl -s -v \
  --cert-type pem \
  --cert $CRT_DIR/hmf.crt \
  --key $CRT_DIR/hmf.key \
  https://api.hartwigmedicalfoundation.nl/hmf/v1/sets \
  -XPOST \
  -H "Content-Type: application/json" \
  -d '{"entity": "'$ENTITY_NAME'", "ini": "'$INI_NAME'", "name": "'$SET_NAME'", "ref_sample": "'$REF_SAMPLE'", "tumor_sample": "'$TUMOR_SAMPLE'"}'

