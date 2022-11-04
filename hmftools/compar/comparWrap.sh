#!/usr/bin/env bash
# Author: Teoman Deger
# Version 0.3, 03-11-2022
# -----------------------

source ~/scripts/functions/message_functions || exit 1

#quickuse
#~/comparWRAP.sh -J ~/data/COMPAR/compar_v1.1beta.jar -R ~/data/TISPA/PT01/TV-PT01-WVT-DNA-WGS/ -N ~/data/TISPA/PT01/TV-PT01-WVTT-DNA-WGS/ -C ALL -L DETAILED -S -O ./COMPAR_OUT

#error message
print_usage(){
    echo ""
    echo "Descr: Function to pass options more directly unto COMPAR"
    echo "Usage: $(basename "$0") -J \${COMPAR-jar location} -R \${Reference Dir} -N \${New Dir} -C \${Category} -L \${Levels}"
    echo "------------------------------------------------------------"
    echo -e "Available Categories:\n -ALL\n -GERMLINE_DELETION\n -COPY_NUMBER\n -SOMATIC_VARIANT\n -DRIVER\n -CUPPA\n -PURITY\n -LILAC\n -FUSION\n -GERMLINE_SV\n -DISRUPTION\n -GERMLINE_VARIANT\n -GENE_COPY_NUMBER\n -CHORD"
    echo -e "\nAvailable Levels:\n -REPORTABLE\n -DETAILED"
    echo -e "\nOptions:"
    echo -e " -S Report every category seperately     only if (-C ALL)"
    echo -e " -X Set output name                      default: parent of RefDir)"
    echo -e " -O Set output directory                 default: current Dir"
    echo "------------------------------------------------------------"
    echo "Exmpl: $(basename "$0")  -J /pathTo/compar.jar  -R /pathTo/refData  -N /pathTo/newData  -C ALL              -L DETAILED    -S"
    echo "       $(basename "$0")  -J /pathTo/compar.jar  -R /pathTo/refData  -N /pathTo/refData  -C SOMATIC_VARIANT  -L REPORTABLE  -O /pathTo/OutputDir"
    echo ""
    exit 1
}

#option settings
sepOutput='false'
outputDir=./
while getopts ':J:R:N:C:L:O:X:S' flag; do
    case "${flag}" in
        J) comparLoc=${OPTARG} ;;
        R) refdir=${OPTARG} ;;
        N) newdir=${OPTARG} ;;
        C) runCats=${OPTARG} ;;
        L) level=${OPTARG} ;;
        O) outputDir=${OPTARG} ;;
        X) expName=${OPTARG} ;;
        S) sepOutput='TRUE' ;;
        *) print_usage >&2
        exit 1 ;;
    esac
done

#sanity checks
if [ ! -d ${outputDir} ]; then
  echo Output directory \"${outputDir}\" doesnt exist, creating
  mkdir -p ${outputDir};
fi

if [[ -z "${comparLoc}" || -z "${refdir}" || -z "${newdir}"  || -z "${runCats}" || -z "${level}" ]]; then
 print_usage; fi

if [[ ${level} != DETAILED ]] && [[ ${level} != REPORTABLE ]]; then
 print_usage; fi

if [[ ${runCats} == ALL ]] && [[ ${sepOutput} == "TRUE" ]]; then
 list=(GERMLINE_DELETION COPY_NUMBER SOMATIC_VARIANT DRIVER CUPPA PURITY LILAC FUSION GERMLINE_SV DISRUPTION GERMLINE_VARIANT GENE_COPY_NUMBER CHORD)
else
 list=$runCats; fi


filesource=\"REF\;sample_dir=$refdir\,NEW\;sample_dir=${newdir}\"
IFS='/' read -r -a arrayref <<< "$refdir"
IFS='/' read -r -a arraynew <<< "$newdir"
rlen=${#arrayref[@]}
nlen=${#arraynew[@]}

if [[ -z "${expName}" ]]; then
expName=${arrayref[$((rlen-2))]}; fi
echo "SampleId,RefSampleId,NewSampleId" > temp_sampID.csv; \
echo "${expName},${arrayref[$((rlen-1))]},${arraynew[$((nlen-1))]}" >> temp_sampID.csv

outid=${expName}_${arrayref[$((rlen-1))]}"_vs_"${arraynew[$((nlen-1))]}


for VALUE in "${list[@]}"
do
  echo -e "\n"
  echo ${VALUE}
  java -jar ${comparLoc}\
       -sample_id_file ./temp_sampID.csv\
       -categories ${VALUE}\
       -match_level ${level}\
       -file_sources ${filesource}\
       -output_dir ${outputDir}\
       -output_id ${VALUE} -log_debug
sed -i -e "s/REF_ONLY/${arrayref[$((rlen-1))]}_only/g" ${outputDir}${expName}.cmp.${VALUE}.combined.csv
sed -i -e "s/NEW_ONLY/${arraynew[$((nlen-1))]}_only/g" ${outputDir}${expName}.cmp.${VALUE}.combined.csv
done

if [[ ${runCats} == ALL ]] && [[ ${sepOutput} == "TRUE" ]]
then
echo -e "\n"
echo ALL COMBINED
java -jar ${comparLoc}\
     -sample_id_file ./temp_sampID.csv\
     -categories $runCats\
     -match_level ${level}\
     -file_sources ${filesource}\
     -output_dir ${outputDir}\
     -output_id ${VALUE} -log_debug
sed -i -e "s/REF_ONLY/${arrayref[$((rlen-1))]}_only/g" ${outputDir}${expName}.cmp.${VALUE}.combined.csv
sed -i -e "s/NEW_ONLY/${arraynew[$((nlen-1))]}_only/g" ${outputDir}${expName}.cmp.${VALUE}.combined.csv
fi

if [[ ${level} == REPORTABLE ]]
then
  echo -e "\n"
  cat ${outputDir}${expName}.cmp.${outid}.combined.csv
fi

rm temp_sampID.csv
