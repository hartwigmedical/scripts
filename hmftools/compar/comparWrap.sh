#!/usr/bin/env bash
# Author: Teoman Deger
# Version 0.3, 03-11-2022
# -----------------------

source ~/scripts/functions/message_functions || exit 1

#TO DO:
#-explain order of CSV and folder structure

#quickuse
#./comparWRAP.sh -J ./data/COMPAR/compar_v1.1beta.jar -I ./data/COMPAR/sample_id_mappings2.csv -D ./data/COLO-data/ -C ALL -L DETAILED -S -O ./COMPAR_OUT
#error message
print_usage(){
    echo ""
    echo "Descr: Function to pass options more directly unto COMPAR"
    echo "Usage: $(basename "$0") -J \${COMPAR-jar location} -I \${Sample-IDs} -C \${Category} -L \${Levels}"
    echo "------------------------------------------------------------"
    echo -e "Available Categories:\n -ALL\n -GERMLINE_DELETION\n -COPY_NUMBER\n -SOMATIC_VARIANT\n -DRIVER\n -CUPPA\n -PURITY\n -LILAC\n -FUSION\n -GERMLINE_SV\n -DISRUPTION\n -GERMLINE_VARIANT\n -GENE_COPY_NUMBER\n -CHORD"
    echo -e "\nAvailable Levels:\n -REPORTABLE\n -DETAILED"
    echo -e "\nOptions:\n -S Report every category seperately, only if (-C ALL)"
    echo "------------------------------------------------------------"
    echo "Exmpl: $(basename "$0")  -J /pathTo/compar.jar  -I /pathTo/SAMPLE_ID.csv  -D /pathTo/data  -C ALL              -L DETAILED    -S"
    echo "       $(basename "$0")  -J /pathTo/compar.jar  -I /pathTo/SAMPLE_ID.csv  -D /pathTo/data  -C SOMATIC_VARIANT  -L REPORTABLE  -O /pathTo/OutputDir"
    echo ""
    exit 1
}

#option settings
sepOutput='false'
outputDir=./
while getopts ':J:I:D:C:L:O:S' flag; do
    case "${flag}" in
        J) comparLoc=${OPTARG} ;;
        I) IDfile=${OPTARG} ;;
        D) pathtodir=${OPTARG} ;;
        C) runCats=${OPTARG} ;;
        L) level=${OPTARG} ;;
        O) outputDir=${OPTARG} ;;
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

if [[ -z "${comparLoc}" || -z "${IDfile}" || -z "${pathtodir}"  || -z "${runCats}" || -z "${level}" ]]; then
 print_usage; fi

if [[ ${level} != DETAILED ]] && [[ ${level} != REPORTABLE ]]; then
 print_usage; fi

if [[ ${runCats} == ALL ]] && [[ ${sepOutput} == "TRUE" ]]; then
 list=(GERMLINE_DELETION COPY_NUMBER SOMATIC_VARIANT DRIVER CUPPA PURITY LILAC FUSION GERMLINE_SV DISRUPTION GERMLINE_VARIANT GENE_COPY_NUMBER CHORD)
else
 list=$runCats; fi

expName=$(tail -n +2 $IDfile | cut -f 1 -d",")
refName=$(tail -n +2 $IDfile | cut -f 2 -d",")
newName=$(tail -n +2 $IDfile | cut -f 3 -d",")
refdir=${pathtodir}${expName}"/"${refName}"/"
newdir=${pathtodir}${expName}"/"${newName}"/"
filesource=\"REF\;sample_dir=$refdir\,NEW\;sample_dir=${newdir}\"
outid=${refName}"_vs_"${newName}

for VALUE in "${list[@]}"
do
  echo -e "\n"
  echo ${VALUE}
  java -jar ${comparLoc}\
       -sample_id_file ${IDfile}\
       -categories ${VALUE}\
       -match_level ${level}\
       -file_sources ${filesource}\
       -output_dir ${outputDir}\
       -output_id ${outid}"_"${VALUE} -log_debug
sed -i -e "s/REF_ONLY/${refName}_only/g" ${outputDir}${expName}.cmp.${outid}_${VALUE}.combined.csv
sed -i -e "s/NEW_ONLY/${newName}_only/g" ${outputDir}${expName}.cmp.${outid}_${VALUE}.combined.csv
done

if [[ ${runCats} == ALL ]] && [[ ${sepOutput} == "TRUE" ]]
then
echo -e "\n"
echo ALL COMBINED
java -jar ${comparLoc}\
     -sample_id_file ${IDfile}\
     -categories $runCats\
     -match_level ${level}\
     -file_sources ${filesource}\
     -output_dir ${outputDir}\
     -output_id ${outid} -log_debug
sed -i -e "s/REF_ONLY/${refName}_only/g" ${outputDir}${expName}.cmp.${outid}.combined.csv
sed -i -e "s/NEW_ONLY/${newName}_only/g" ${outputDir}${expName}.cmp.${outid}.combined.csv
fi

if [[ ${level} == REPORTABLE ]]
then
  echo -e "\n"
  cat ${outputDir}${expName}.cmp.${outid}.combined.csv
fi