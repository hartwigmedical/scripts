#!/usr/bin/env bash
# Author: Teoman Deger
# -----------------------

source message_functions || exit 1


workingDir=""
folderNames=()
hgVersion="hg38"
outputDir=""
bcftools="/data/tools/bcftools/1.9/bcftools"

script="$(basename "$0")"

print_usage() {
cat <<EOM
-----
 Descr: Compares pipeline runs and creates summarizing plots of differences
 Usage: $script [required arguments] [optional arguments]
 Required arguments:
   -i [s]    Input directories, separated by ';'. Provide either 2 or 4 such directories.
   -w [s]    Working directory
 Optional arguments:
   -r [s]    Reference genome version. Either h19 or h38. Defaults to hg38.
   -o [s]    Output directory. Defaults to 'expPlots' subdirectory of working directory.
   -b [s]    Path to bcftools. Defaults to '/data/tools/bcftools/1.9/bcftools'.
-----
EOM
}

if [[ $# -eq 0 || $1 == "-h" || $1 == "--help" ]]; then
    print_usage
    exit 1
fi

while getopts ":i:w:r:o:b:" opt; do
    case ${opt} in
        i)
            IFS=';' read -ra folderNames <<< "${OPTARG}";;
        w)
            workingDir="${OPTARG}";;
        r)
            hgVersion="${OPTARG}";;
        o)
            outputDir="${OPTARG}";;
        b)
            bcftools="${OPTARG}";;
        \?)
            print_usage
            die "Invalid option: -${OPTARG}";;
        :)
            print_usage
            die "Option -${OPTARG} requires an argument.";;
    esac
done

runMode=${#folderNames[@]}

## sanity checks
[[ -n "${workingDir}" ]] || (print_usage && die "No working directory provided.")
if [[ "${runMode}" != "2" && "${runMode}" != "4" ]]; then
    print_usage
    die "Invalid number of input directories provided: ${hgVersion}"
fi
if [[ "${hgVersion}" != "hg19" && "${hgVersion}" != "hg38" ]]; then
    print_usage
    die "Invalid reference genome version: ${hgVersion}"
fi

if [[ -z "${outputDir}" ]]; then
    outputDir="${workingDir}/expPlots"
fi

## realpath everything relevant
for i in "${!folderNames[@]}"; do
    folderNames[i]=$(realpath "${folderNames[$i]}")
done
workingDir=$(realpath "${workingDir}")
outputDir=$(realpath "${outputDir}")
scriptDir=$(dirname "$(realpath "$0")")

echo "################################################################################################"
echo "#                                                                                              #"
echo "#   Cobalt plotter does not work for runMode=4                                                 #"
echo "#   All R-scripts are silenced, error messages wont show, it makes the output better to read   #"
echo "#   The R-scripts have some folders hard-coded                                                 #"
echo "#   esvee and gridss are currently not compatible for germline SVs                             #"
echo "#                                                                                              #"
echo "################################################################################################"

file1=${folderNames[0]}
file2=${folderNames[1]}
name1=$(echo "${folderNames[0]}" | awk -F '/' '{ print $(NF-1) }')
name2=$(echo "${folderNames[1]}" | awk -F '/' '{ print $(NF-1) }')

isecDir="${workingDir}/isec"

if [ ! -d "${outputDir}" ]; then
  echo "Output directory '${outputDir}' doesnt exist, creating"
  mkdir -p "${outputDir}";
fi

echo "runMode is ${runMode}"

if [ "${runMode}" -eq 4 ]; then
  file3=${folderNames[2]}
  file4=${folderNames[3]}
  name3=$(echo "${folderNames[2]}" | awk -F '/' '{ print $(NF-1) }')
  name4=$(echo "${folderNames[3]}" | awk -F '/' '{ print $(NF-1) }')
  echo "cobalt-plotter wont run in runmode=4"
fi


if [ "${runMode}" -eq 2 ]; then
  cobDat1="${file1}/cobalt/*.cobalt.ratio.tsv.gz"
  cobDat1L=$(ls "$cobDat1")
  cdName=$(basename "${cobDat1L}" .cobalt.ratio.tsv.gz)
  cobSeg1="${file1}/cobalt/${cdName}.cobalt.ratio.pcf"

  cobDat2="${file2}/cobalt/*.cobalt.ratio.tsv.gz"
  cobDat2L=$(ls "$cobDat2")
  cdName=$(basename "${cobDat2L}" .cobalt.ratio.tsv.gz)
  cobSeg2="${file2}/cobalt/${cdName}.cobalt.ratio.pcf"

  echo "Creating a Cobalt distribution plot"
  "${scriptDir}/cnpPlot.R" "$cobDat1" "$cobSeg1" "$cobDat2" "$cobSeg2" "${hgVersion}" "${workingDir}" \
        "${scriptDir}/genome_length.txt" > /dev/null 2>&1
  mv "${workingDir}/CobaltDistPlot.png" "${workingDir}/CobaltPlot_${name1}vs${name2}.png"
  mv "${workingDir}/CobaltPlot_${name1}vs${name2}.png" "${outputDir}"
fi

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

echo -e "\nSVs"
echo "${name1} vs ${name2}"
$bcftools isec -f "PASS" "${file1}/purple/"*.purple.sv.vcf.gz "${file2}/purple/"*.purple.sv.vcf.gz -p "${isecDir}"

count_event_types() {
  local file=$1
  local prefix=$2

  # Determine whether to use EVENTTYPE or SVTYPE based on VCF header
  if grep -q "##esveeVersion=" "${file}"; then
    #echo "Using SVTYPE for event counting"
    event_type="SVTYPE"
  elif grep -q "##gridssVersion=" "${file}"; then
    #echo "Using EVENTTYPE for event counting"
    event_type="EVENTTYPE"
  else
    echo "Neither esveeVersion nor gridssVersion found, defaulting to SVTYPE"
    event_type="SVTYPE"
  fi

  pull0=$(grep -c -v "^#" "${file}")
  pull1=$(grep -c -e "${event_type}=DEL" "${file}")
  pull2=$(grep -c -e "${event_type}=SGL" "${file}")
  pull3=$(grep -c -e "${event_type}=BND" "${file}")
  pull4=$(grep -c -e "${event_type}=INV" "${file}")
  pull5=$(grep -c -e "${event_type}=DUP" "${file}")
  pull6=$(grep -c -e "${event_type}=INS" "${file}")
  echo -e "$prefix\t${pull1}\t${pull2}\t${pull3}\t${pull4}\t${pull5}\t${pull6}\t${pull0}"
}

echo -e "Name\tDel\tSgl\tBnd\tInv\tDup\tIns\tTotal"
count_event_types "${isecDir}/0002.vcf" "Shared"
count_event_types "${isecDir}/0000.vcf" "$name1"
count_event_types "${isecDir}/0001.vcf" "$name2"


# Generate and index the shared VCF
$bcftools view "${isecDir}/0002.vcf" -Oz -o "${workingDir}/shared_sv_12.vcf.gz"
$bcftools index "${workingDir}/shared_sv_12.vcf.gz"

if [ "${runMode}" -eq 4 ]; then
  # Process additional comparisons
  $bcftools isec -f "PASS" "${file3}/purple/"*.purple.sv.vcf.gz "${file4}/purple/"*.purple.sv.vcf.gz -p "${isecDir}"
  echo -e "Name\tDel\tSgl\tBnd\tInv\tDup\tIns\tTotal"
  count_event_types "${isecDir}/0002.vcf" "Shared"
  count_event_types "${isecDir}/0000.vcf" "$name3"
  count_event_types "${isecDir}/0001.vcf" "$name4"

  $bcftools view "${isecDir}/0002.vcf" -Oz -o "${workingDir}/shared_sv_34.vcf.gz"
  $bcftools index "${workingDir}/shared_sv_34.vcf.gz"

  echo "shared_12 vs shared_34"
  $bcftools isec -f "PASS" "${workingDir}/shared_sv_12.vcf.gz" "${workingDir}/shared_sv_34.vcf.gz" -p "${isecDir}"
  echo -e "Name\tDel\tSgl\tBnd\tInv\tDup\tIns\tTotal"
  count_event_types "${isecDir}/0002.vcf" "Shared"
  count_event_types "${isecDir}/0000.vcf" "$name3"
  count_event_types "${isecDir}/0001.vcf" "$name4"
fi

rm "${workingDir}/shared_"*

#######################################################################################
#######################################################################################
#######################################################################################

echo -e "\nGermline SVs"
echo "${name1} vs ${name2}"
$bcftools isec -f "PASS" "${file1}/purple/"*.purple.sv.germline.vcf.gz "${file2}/purple/"*.purple.sv.germline.vcf.gz -p "${isecDir}"

echo -e "Name\tDel\tSgl\tBnd\tInv\tDup\tIns\tTotal"
count_event_types "${isecDir}/0002.vcf" "Shared"
count_event_types "${isecDir}/0000.vcf" "$name1"
count_event_types "${isecDir}/0001.vcf" "$name2"

# Generate and index the shared VCF
$bcftools view "${isecDir}/0002.vcf" -Oz -o "${workingDir}/shared_sv_12.vcf.gz"
$bcftools index "${workingDir}/shared_sv_12.vcf.gz"

if [ "${runMode}" -eq 4 ]; then
  # Process additional comparisons
  $bcftools isec -f "PASS" "${file3}/purple/"*.purple.sv.germline.vcf.gz "${file4}/purple/"*.purple.sv.germline.vcf.gz -p "${isecDir}"
  echo -e "Name\tDel\tSgl\tBnd\tInv\tDup\tIns\tTotal"
  count_event_types "${isecDir}/0002.vcf" "Shared"
  count_event_types "${isecDir}/0000.vcf" "$name3"
  count_event_types "${isecDir}/0001.vcf" "$name4"

  $bcftools view "${isecDir}/0002.vcf" -Oz -o "${workingDir}/shared_sv_34.vcf.gz"
  $bcftools index "${workingDir}/shared_sv_34.vcf.gz"

  echo "shared_12 vs shared_34"
  $bcftools isec -f "PASS" "${workingDir}/shared_sv_12.vcf.gz" "${workingDir}/shared_sv_34.vcf.gz" -p "${isecDir}"
  echo -e "Name\tDel\tSgl\tBnd\tInv\tDup\tIns\tTotal"
  count_event_types "${isecDir}/0002.vcf" "Shared"
  count_event_types "${isecDir}/0000.vcf" "$name3"
  count_event_types "${isecDir}/0001.vcf" "$name4"
fi

rm "${workingDir}/shared_"*

#######################################################################################
#######################################################################################
#######################################################################################

#identify uniques and shared 1&2 + plots
echo -e "\nGermline Variants"
$bcftools isec -f "PASS" "${file1}/purple/"*.purple.germline.vcf.gz "${file2}/purple/"*.purple.germline.vcf.gz -p "${isecDir}"
refOnly=$(grep -c -v "^#" "${isecDir}/0000.vcf")  #refOnly
shared=$(grep -c -v "^#" "${isecDir}/0002.vcf")  #shared_refdata
newOnly=$(grep -c -v "^#" "${isecDir}/0001.vcf")  #newOnly
echo -e "${name1}_only\tShared\t${name2}_only"
echo -e "${refOnly}\t${shared}\t${newOnly}"
#cat ~/coloOldvNew/isec/0003.vcf | grep -c -v "^#"  #shared_newdata
$bcftools view "${isecDir}/0002.vcf" -Oz -o "${workingDir}/shared_gl_12.vcf.gz"
$bcftools index "${workingDir}/shared_gl_12.vcf.gz"
"${scriptDir}/plotter.R" "$name1" "$name2" GLVars "${isecDir}" > /dev/null 2>&1

if [ "${runMode}" -eq 4 ]; then
  #identify uniques and shared 3&4 + plots
  $bcftools isec -f "PASS" "${file3}/purple/"*.purple.germline.vcf.gz "${file4}/purple/"*.purple.germline.vcf.gz -p "${isecDir}"
  refOnly=$(grep -c -v "^#" "${isecDir}/0000.vcf")  #refOnly
  shared=$(grep -c -v "^#" "${isecDir}/0002.vcf")  #shared_refdata
  newOnly=$(grep -c -v "^#" "${isecDir}/0001.vcf")  #newOnly
  echo -e "${name3}_only\tShared\t${name4}_only"
  echo -e "${refOnly}\t${shared}\t${newOnly}"
  #cat ~/coloOldvNew/isec/0003.vcf | grep -c -v "^#"  #shared_newdata
  $bcftools view "${isecDir}/0002.vcf" -Oz -o "${workingDir}/shared_gl_34.vcf.gz"
  $bcftools index "${workingDir}/shared_gl_34.vcf.gz"
  "${scriptDir}/plotter.R" "$name3" "$name4" GLVars "${isecDir}" > /dev/null 2>&1

  #identify batch_specific_uniques and total_shareds
  $bcftools isec -f "PASS" "${workingDir}/shared_gl_12.vcf.gz" "${workingDir}/shared_gl_34.vcf.gz" -p "${isecDir}"
  refOnly=$(grep -c -v "^#" "${isecDir}/0000.vcf")  #refOnly
  shared=$(grep -c -v "^#" "${isecDir}/0002.vcf")  #shared_refdata
  newOnly=$(grep -c -v "^#" "${isecDir}/0001.vcf")  #newOnly
  echo -e "${name1}/${name2}_only\tShared all\t${name3}/${name4}_only"
  echo -e "${refOnly}\t${shared}\t${newOnly}"

  #extract all final Vars
  $bcftools view "${isecDir}/0000.vcf"-Oz -o "${workingDir}/only_12.vcf.gz"
  $bcftools index "${workingDir}/only_12.vcf.gz"
  $bcftools view "${isecDir}/0001.vcf"-Oz -o "${workingDir}/only_34.vcf.gz"
  $bcftools index "${workingDir}/only_34.vcf.gz"
  $bcftools view "${isecDir}/0002.vcf" -Oz -o "${workingDir}/shared_all.vcf.gz"
  $bcftools index "${workingDir}/shared_all.vcf.gz"

  #extract old-only regions from old-files
  $bcftools isec -f "PASS" "${workingDir}/only_12.vcf.gz" "${file1}/purple/"*.purple.germline.vcf.gz -p "${isecDir}"
  $bcftools view "${isecDir}/0003.vcf"-Oz -o "${workingDir}/only_12_1.vcf.gz"
  $bcftools index "${workingDir}/only_12_1.vcf.gz"
  $bcftools isec -f "PASS" "${workingDir}/only_12.vcf.gz" "${file2}/purple/"*.purple.germline.vcf.gz -p "${isecDir}"
  $bcftools view "${isecDir}/0003.vcf"-Oz -o "${workingDir}/only_12_2.vcf.gz"
  $bcftools index "${workingDir}/only_12_2.vcf.gz"

  #extract new-only regions from new-files
  $bcftools isec -f "PASS" "${workingDir}/only_34.vcf.gz" "${file3}/purple/"*.purple.germline.vcf.gz -p "${isecDir}"
  $bcftools view "${isecDir}/0003.vcf"-Oz -o "${workingDir}/only_34_3.vcf.gz"
  $bcftools index "${workingDir}/only_34_3.vcf.gz"
  $bcftools isec -f "PASS" "${workingDir}/only_34.vcf.gz" "${file4}/purple/"*.purple.germline.vcf.gz -p "${isecDir}"
  $bcftools view "${isecDir}/0003.vcf"-Oz -o "${workingDir}/only_34_4.vcf.gz"
  $bcftools index "${workingDir}/only_34_4.vcf.gz"

  #extract shared lists from all
  $bcftools isec -f "PASS" "${workingDir}/shared_all.vcf.gz" "${file1}/purple/"*.purple.germline.vcf.gz -p "${isecDir}"
  $bcftools view "${isecDir}/0003.vcf"-Oz -o "${workingDir}/shared_all_1.vcf.gz"
  $bcftools index "${workingDir}/shared_all_1.vcf.gz"
  $bcftools isec -f "PASS" "${workingDir}/shared_all.vcf.gz" "${file2}/purple/"*.purple.germline.vcf.gz -p "${isecDir}"
  $bcftools view "${isecDir}/0003.vcf"-Oz -o "${workingDir}/shared_all_2.vcf.gz"
  $bcftools index "${workingDir}/shared_all_2.vcf.gz"
  $bcftools isec -f "PASS" "${workingDir}/shared_all.vcf.gz" "${file3}/purple/"*.purple.germline.vcf.gz -p "${isecDir}"
  $bcftools view "${isecDir}/0003.vcf"-Oz -o "${workingDir}/shared_all_3.vcf.gz"
  $bcftools index "${workingDir}/shared_all_3.vcf.gz"
  $bcftools isec -f "PASS" "${workingDir}/shared_all.vcf.gz" "${file4}/purple/"*.purple.germline.vcf.gz -p "${isecDir}"
  $bcftools view "${isecDir}/0003.vcf"-Oz -o "${workingDir}/shared_all_4.vcf.gz"
  $bcftools index "${workingDir}/shared_all_4.vcf.gz"

  "${scriptDir}/shrdPlotter.R" GLVars "${workingDir}" > /dev/null 2>&1
fi

mv "${workingDir}/"*.png "${outputDir}" 2>/dev/null
mv "${isecDir}/"*.png "${outputDir}" 2>/dev/null
rm "${workingDir}/shared_"* 2>/dev/null
rm "${workingDir}/only_"* 2>/dev/null


#######################################################################################
#######################################################################################

#identify uniques and shared 1&2 + plots
echo -e "\nSomatic Variants"
$bcftools isec -f "PASS" "${file1}/purple/"*.purple.somatic.vcf.gz "${file2}/purple/"*.purple.somatic.vcf.gz -p "${isecDir}"
refOnly=$(grep -c -v "^#" "${isecDir}/0000.vcf")  #refOnly
shared=$(grep -c -v "^#" "${isecDir}/0002.vcf")  #shared_refdata
newOnly=$(grep -c -v "^#" "${isecDir}/0001.vcf")  #newOnly
echo -e "${name1}_only\tShared\t${name2}_only"
echo -e "${refOnly}\t${shared}\t${newOnly}"
#cat ~/coloOldvNew/isec/0003.vcf | grep -c -v "^#"  #shared_newdata
$bcftools view "${isecDir}/0002.vcf" -Oz -o "${workingDir}/shared_som_12.vcf.gz"
$bcftools index "${workingDir}/shared_som_12.vcf.gz"
"${scriptDir}/plotter.R" "$name1" "$name2" SomVars "${isecDir}" > /dev/null 2>&1
"${scriptDir}/triNucPlot.R" "$name1" "$name2" SomVars "$hgVersion" "${isecDir}" > /dev/null 2>&1

if [ "${runMode}" -eq 4 ]; then
  #identify uniques and shared 3&4 + plots
  $bcftools isec -f "PASS" "${file3}/purple/"*.purple.somatic.vcf.gz "${file4}/purple/"*.purple.somatic.vcf.gz -p "${isecDir}"
  refOnly=$(grep -c -v "^#" "${isecDir}/0000.vcf")  #refOnly
  shared=$(grep -c -v "^#" "${isecDir}/0002.vcf")  #shared_refdata
  newOnly=$(grep -c -v "^#" "${isecDir}/0001.vcf")  #newOnly
  echo -e "${name3}_only\tShared\t${name4}_only"
  echo -e "${refOnly}\t${shared}\t${newOnly}"
  #cat ~/coloOldvNew/isec/0003.vcf | grep -c -v "^#"  #shared_newdata
  $bcftools view "${isecDir}/0002.vcf" -Oz -o "${workingDir}/shared_som_34.vcf.gz"
  $bcftools index "${workingDir}/shared_som_34.vcf.gz"
  "${scriptDir}/plotter.R" "$name3" "$name4" SomVars "${isecDir}" > /dev/null 2>&1
  "${scriptDir}/triNucPlot.R" "$name3" "$name4" SomVars "$hgVersion" "${isecDir}" > /dev/null 2>&1

  #identify batch_specific_uniques and total_shareds
  $bcftools isec -f "PASS" "${workingDir}/shared_som_12.vcf.gz" "${workingDir}/shared_som_34.vcf.gz" -p "${isecDir}"
  refOnly=$(grep -c -v "^#" "${isecDir}/0000.vcf")  #refOnly
  shared=$(grep -c -v "^#" "${isecDir}/0002.vcf")  #shared_refdata
  newOnly=$(grep -c -v "^#" "${isecDir}/0001.vcf")  #newOnly
  echo -e "${name1}/${name2}_only\tShared all\t${name3}/${name4}_only"
  echo -e "${refOnly}\t${shared}\t${newOnly}"
  "${scriptDir}/triNucPlot.R" Old_shrd New_shrd SomVars "$hgVersion" "${isecDir}" > /dev/null 2>&1

  #extract all final Vars
  $bcftools view "${isecDir}/0000.vcf"-Oz -o "${workingDir}/only_12.vcf.gz"
  $bcftools index "${workingDir}/only_12.vcf.gz"
  $bcftools view "${isecDir}/0001.vcf"-Oz -o "${workingDir}/only_34.vcf.gz"
  $bcftools index "${workingDir}/only_34.vcf.gz"
  $bcftools view "${isecDir}/0002.vcf"-Oz -o "${workingDir}/shared_all.vcf.gz"
  $bcftools index "${workingDir}/shared_all.vcf.gz"

  #extract old-only regions from old-files
  $bcftools isec -f "PASS" "${workingDir}/only_12.vcf.gz" "${file1}/purple/"*.purple.somatic.vcf.gz -p "${isecDir}"
  $bcftools view "${isecDir}/0003.vcf"-Oz -o "${workingDir}/only_12_1.vcf.gz"
  $bcftools index "${workingDir}/only_12_1.vcf.gz"
  $bcftools isec -f "PASS" "${workingDir}/only_12.vcf.gz" "${file2}/purple/"*.purple.somatic.vcf.gz -p "${isecDir}"
  $bcftools view "${isecDir}/0003.vcf"-Oz -o "${workingDir}/only_12_2.vcf.gz"
  $bcftools index "${workingDir}/only_12_2.vcf.gz"

  #extract new-only regions from new-files
  $bcftools isec -f "PASS" "${workingDir}/only_34.vcf.gz" "${file3}/purple/"*.purple.somatic.vcf.gz -p "${isecDir}"
  $bcftools view "${isecDir}/0003.vcf"-Oz -o "${workingDir}/only_34_3.vcf.gz"
  $bcftools index "${workingDir}/only_34_3.vcf.gz"
  $bcftools isec -f "PASS" "${workingDir}/only_34.vcf.gz" "${file4}/purple/"*.purple.somatic.vcf.gz -p "${isecDir}"
  $bcftools view "${isecDir}/0003.vcf"-Oz -o "${workingDir}/only_34_4.vcf.gz"
  $bcftools index "${workingDir}/only_34_4.vcf.gz"

  #extract shared lists from all
  $bcftools isec -f "PASS" "${workingDir}/shared_all.vcf.gz" "${file1}/purple/"*.purple.somatic.vcf.gz -p "${isecDir}"
  $bcftools view "${isecDir}/0003.vcf"-Oz -o "${workingDir}/shared_all_1.vcf.gz"
  $bcftools index "${workingDir}/shared_all_1.vcf.gz"
  $bcftools isec -f "PASS" "${workingDir}/shared_all.vcf.gz" "${file2}/purple/"*.purple.somatic.vcf.gz -p "${isecDir}"
  $bcftools view "${isecDir}/0003.vcf"-Oz -o "${workingDir}/shared_all_2.vcf.gz"
  $bcftools index "${workingDir}/shared_all_2.vcf.gz"
  $bcftools isec -f "PASS" "${workingDir}/shared_all.vcf.gz" "${file3}/purple/"*.purple.somatic.vcf.gz -p "${isecDir}"
  $bcftools view "${isecDir}/0003.vcf"-Oz -o "${workingDir}/shared_all_3.vcf.gz"
  $bcftools index "${workingDir}/shared_all_3.vcf.gz"
  $bcftools isec -f "PASS" "${workingDir}/shared_all.vcf.gz" "${file4}/purple/"*.purple.somatic.vcf.gz -p "${isecDir}"
  $bcftools view "${isecDir}/0003.vcf"-Oz -o "${workingDir}/shared_all_4.vcf.gz"
  $bcftools index "${workingDir}/shared_all_4.vcf.gz"

  "${scriptDir}/shrdPlotter.R" SomVars "${workingDir}" > /dev/null 2>&1
  "${scriptDir}/shrdNucPlotter.R" "$name1" "$name2" "$name3" "$name4" SomVars "$hgVersion" "${workingDir}" > /dev/null 2>&1
fi

mv "${workingDir}/"*.png "${outputDir}" 2>/dev/null
mv "${isecDir}/"*.png "${outputDir}" 2>/dev/null
rm "${workingDir}/shared_"* 2>/dev/null
rm "${workingDir}/only_"* 2>/dev/null


#######################################################################################
#######################################################################################
#######################################################################################

echo -e "\nSomatic Drivers"
$bcftools filter -o "${workingDir}/shared_file1_somdr.vcf.gz" -O z -i 'INFO/REPORTED="0"' "${file1}/purple/"*.purple.somatic.vcf.gz
$bcftools index "${workingDir}/shared_file1_somdr.vcf.gz"
$bcftools filter -o "${workingDir}/shared_file2_somdr.vcf.gz" -O z -i 'INFO/REPORTED="0"' "${file2}/purple/"*.purple.somatic.vcf.gz
$bcftools index "${workingDir}/shared_file2_somdr.vcf.gz"
$bcftools isec "${workingDir}/shared_file1_somdr.vcf.gz" "${workingDir}/shared_file2_somdr.vcf.gz" -p "${isecDir}"
SomDrR=$(grep -c -v "^#" "${isecDir}/0000.vcf")  #refOnly
SomDrN=$(grep -c -v "^#" "${isecDir}/0001.vcf")  #newOnly
SomDrS=$(grep -c -v "^#" "${isecDir}/0002.vcf")  #shared_refdata
echo -e "${name1}_only\tShared\t${name2}_only"
echo -e "${SomDrR}\t${SomDrS}\t${SomDrN}"
#cat ~/coloOldvNew/isec/0003.vcf | grep -c -v "^#"  #shared_newdata
$bcftools view "${isecDir}/0002.vcf"-Oz -o "${workingDir}/shared_somdri_12.vcf.gz"
$bcftools index "${workingDir}/shared_somdri_12.vcf.gz"

if [ "${runMode}" -eq 4 ]; then
  $bcftools filter -o "${workingDir}/shared_file3_somdr.vcf.gz" -O z -i 'INFO/REPORTED="0"' "${file3}/purple/"*.purple.somatic.vcf.gz
  $bcftools index "${workingDir}/shared_file3_somdr.vcf.gz"
  $bcftools filter -o "${workingDir}/shared_file4_somdr.vcf.gz" -O z -i 'INFO/REPORTED="0"' "${file4}/purple/"*.purple.somatic.vcf.gz
  $bcftools index "${workingDir}/shared_file4_somdr.vcf.gz"
  $bcftools isec "${workingDir}/shared_file3_somdr.vcf.gz" "${workingDir}/shared_file4_somdr.vcf.gz" -p "${isecDir}"
  SomDrR=$(grep -c -v "^#" "${isecDir}/0000.vcf")  #refOnly
  SomDrN=$(grep -c -v "^#" "${isecDir}/0001.vcf")  #newOnly
  SomDrS=$(grep -c -v "^#" "${isecDir}/0002.vcf")  #shared_refdata
  echo -e "${name1}_only\tShared\t${name2}_only"
  echo -e "${SomDrR}\t${SomDrS}\t${SomDrN}"
  #cat ~/coloOldvNew/isec/0003.vcf | grep -c -v "^#"  #shared_newdata
  $bcftools view "${isecDir}/0002.vcf"-Oz -o "${workingDir}/shared_somdri_34.vcf.gz"
  $bcftools index "${workingDir}/shared_somdri_34.vcf.gz"

  $bcftools isec -f "PASS" "${workingDir}/shared_somdri_12.vcf.gz" "${workingDir}/shared_somdri_34.vcf.gz" -p "${isecDir}"
  SomDrR=$(grep -c -v "^#" "${isecDir}/0000.vcf")  #refOnly
  SomDrN=$(grep -c -v "^#" "${isecDir}/0001.vcf")  #newOnly
  SomDrS=$(grep -c -v "^#" "${isecDir}/0002.vcf")  #shared_refdata
  echo -e "${name1}/${name2}_only\tShared all\t${name3}/${name4}_only"
  echo -e "${SomDrR}\t${SomDrS}\t${SomDrN}"

fi

mv "${workingDir}/"*.png "${outputDir}" 2>/dev/null
mv "${isecDir}/"*.png "${outputDir}" 2>/dev/null
rm "${workingDir}/shared_"* 2>/dev/null
rm "${workingDir}/only_"* 2>/dev/null

#######################################################################################
#######################################################################################
#######################################################################################

echo -e "\nCN Drivers"
cnd1="${workingDir}/cnd1.txt"
cnd2="${workingDir}/cnd2.txt"
cat "${file1}/purple/"*.driver.catalog.somatic.tsv | grep -v MUTATION | tail -n +2 | awk '{ print $1, $2, $3, $4, $6 }' | sort > "${cnd1}"
cat "${file2}/purple/"*.driver.catalog.somatic.tsv | grep -v MUTATION | tail -n +2 | awk '{ print $1, $2, $3, $4, $6 }' | sort > "${cnd2}"

# Find unique and shared entries using comm
CNDr=$(comm -23 "${cnd1}" "${cnd2}" | wc -l)  # Unique in CND1
CNDn=$(comm -13 "${cnd1}" "${cnd2}" | wc -l)  # Unique in CND2
CNDs=$(comm -12 "${cnd1}" "${cnd2}" | wc -l)  # Shared entries
shared_12="${workingDir}/shared_12.txt"
comm -12 "${cnd1}" "${cnd2}" > "${shared_12}"

echo -e "${name1}\tShared\t ${name2}"
echo -e "${CNDr}\t${CNDs}\t${CNDn}"
rm "${cnd1}" "${cnd2}"

if [ "${runMode}" -eq 4 ]; then
  cnd3="${workingDir}/cnd3.txt"
  cnd4="${workingDir}/cnd4.txt"
  cat "${file3}/purple/"*.driver.catalog.somatic.tsv | grep -v MUTATION | tail -n +2 | awk '{ print $1, $2, $3, $4, $6 }' | sort > "${cnd3}"
  cat "${file4}/purple/"*.driver.catalog.somatic.tsv | grep -v MUTATION | tail -n +2 | awk '{ print $1, $2, $3, $4, $6 }' | sort > "${cnd4}"

  # Find unique and shared entries using comm
  CNDr=$(comm -23 "${cnd3}" "${cnd4}" | wc -l)  # Unique in CND1
  CNDn=$(comm -13 "${cnd3}" "${cnd4}" | wc -l)  # Unique in CND2
  CNDs=$(comm -12 "${cnd3}" "${cnd4}" | wc -l)  # Shared entries
  shared_34="${workingDir}/shared_34.txt"
  comm -12 "${cnd3}" "${cnd4}" > "${shared_34}"

  echo -e "${name3}\tShared\t ${name4}"
  echo -e "${CNDr}\t${CNDs}\t${CNDn}"
  rm "${cnd3}" "${cnd4}"

  # Find unique and shared entries using comm
  CNDr=$(comm -23 "${shared_12}" "${shared_34}" | wc -l)  # Unique in CND1
  CNDn=$(comm -13 "${shared_12}" "${shared_34}" | wc -l)  # Unique in CND2
  CNDs=$(comm -12 "${shared_12}" "${shared_34}" | wc -l)  # Shared entries
  echo -e "${name1}/${name2}_only\tShared all\t${name3}/${name4}_only"
  echo -e "${CNDr}\t${CNDs}\t${CNDn}"

fi

mv "${workingDir}/"*.png "${outputDir}" 2>/dev/null
mv "${isecDir}/"*.png "${outputDir}" 2>/dev/null
rm "${workingDir}/shared_"* 2>/dev/null
rm "${workingDir}/only_"* 2>/dev/null


#######################################################################################
#######################################################################################
#######################################################################################
