#!/usr/bin/env bash
# Author: Teoman Deger
# -----------------------

bcftools="/data/tools/bcftools/1.9/bcftools"
myDir="/home/tdeger/vCompar"
isecDir="${myDir}/isec"
hgVersion="hg19"
runMode=2

#folderNames=(
# "/home/tdeger/coloOldvNew/COLOv2-5.31/"
# "/home/tdeger/coloOldvNew/COLOv3-5.31/"
# "/home/tdeger/coloOldvNew/COLOv6-5.31/"
# "/home/tdeger/coloOldvNew/COLOv8-5.31/"
#)

folderNames=(
"/home/dkoetsier/genomeScan/ref/COLO829v012T-genomescan/"
"/home/dkoetsier/genomeScan/new/COLO829-saas-genomescan/"
)


echo "################################################################################################"
echo "#                                                                                              #"
echo "#   Cobalt plotter does not work for runMode=4                                                 #"
echo "#   All R-scripts are silenced, error messages wont show, it makes the output better to read   #"
echo "#   The R-scripts have some folders hard-coded                                                 #"
echo "#   esvee and gridds are currently not compatible for germline SVs                             #"
echo "#                                                                                              #"
echo "################################################################################################"

file1=$(echo ${folderNames[0]})
file2=$(echo ${folderNames[1]})
name1=$(echo ${folderNames[0]} | awk -F '/' '{ print $(NF-1) }')
name2=$(echo ${folderNames[1]} | awk -F '/' '{ print $(NF-1) }')


outputDir="expPlots"
if [ ! -d ${outputDir} ]; then
  echo Output directory \"${outputDir}\" doesnt exist, creating
  mkdir -p ${outputDir};
fi

echo "runMode is ${runMode}"
if [ ! ${runMode} -eq 2 ] && [ ! ${runMode} -eq 4 ]; then
  echo "runMode must be 2 or 4"
fi

if [ ${runMode} -eq 4 ]; then
  file3=$(echo ${folderNames[2]})
  file4=$(echo ${folderNames[3]})
  name3=$(echo ${folderNames[2]} | awk -F '/' '{ print $(NF-1) }')
  name4=$(echo ${folderNames[3]} | awk -F '/' '{ print $(NF-1) }')
  echo "cobalt-plotter wont run in runmode=4"
fi


if [ ${runMode} -eq 2 ]; then
  cobDat1="${file1}cobalt/*.cobalt.ratio.tsv.gz"
  cobDat1L=$(ls $cobDat1)
  cdName=$(basename "${cobDat1L}" .cobalt.ratio.tsv.gz)
  cobSeg1="${file1}cobalt/${cdName}.cobalt.ratio.pcf"

  cobDat2="${file2}cobalt/*.cobalt.ratio.tsv.gz"
  cobDat2L=$(ls $cobDat2)
  cdName=$(basename "${cobDat2L}" .cobalt.ratio.tsv.gz)
  cobSeg2="${file2}cobalt/${cdName}.cobalt.ratio.pcf"

  echo "Creating a Cobalt distribution plot"
#  ./cnpPlot.R $cobDat1 $cobSeg1 $cobDat2 $cobSeg2 ${hgVersion} ${myDir} > /dev/null 2>&1
#  mv CobaltDistPlot.png CobaltPlot_${name1}vs${name2}.png
#  mv CobaltPlot_${name1}vs${name2}.png ${outputDir}
fi

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

echo -e "\nSVs"
echo "${name1} vs ${name2}"
$bcftools isec -f "PASS" ${file1}/purple/*.purple.sv.vcf.gz ${file2}/purple/*.purple.sv.vcf.gz -p isec

count_event_types() {
  local file=$1
  local prefix=$2

  # Determine whether to use EVENTTYPE or SVTYPE based on VCF header
  if grep -q "##esveeVersion=" ${file}; then
    #echo "Using SVTYPE for event counting"
    event_type="SVTYPE"
  elif grep -q "##gridssVersion=" ${file}; then
    #echo "Using EVENTTYPE for event counting"
    event_type="EVENTTYPE"
  else
    echo "Neither esveeVersion nor gridssVersion found, defaulting to SVTYPE"
    event_type="SVTYPE"
  fi

  pull0=$(cat ${file} | grep -v "^#" | wc -l)
  pull1=$(cat ${file} | grep -c -e "${event_type}=DEL")
  pull2=$(cat ${file} | grep -c -e "${event_type}=SGL")
  pull3=$(cat ${file} | grep -c -e "${event_type}=BND")
  pull4=$(cat ${file} | grep -c -e "${event_type}=INV")
  pull5=$(cat ${file} | grep -c -e "${event_type}=DUP")
  pull6=$(cat ${file} | grep -c -e "${event_type}=INS")
  echo -e "$prefix\t${pull1}\t${pull2}\t${pull3}\t${pull4}\t${pull5}\t${pull6}\t${pull0}"
}

echo -e "Name\tDel\tSgl\tBnd\tInv\tDup\tIns\tTotal"
count_event_types "${myDir}/isec/0002.vcf" "Shared"
count_event_types "${myDir}/isec/0000.vcf" "$name1"
count_event_types "${myDir}/isec/0001.vcf" "$name2"


# Generate and index the shared VCF
$bcftools view ${myDir}/isec/0002.vcf -Oz -o ${myDir}/shared_sv_12.vcf.gz
$bcftools index ${myDir}/shared_sv_12.vcf.gz

if [ ${runMode} -eq 4 ]; then
  # Process additional comparisons
  $bcftools isec -f "PASS" ${file3}/purple/*.purple.sv.vcf.gz ${file4}/purple/*.purple.sv.vcf.gz -p isec
  echo -e "Name\tDel\tSgl\tBnd\tInv\tDup\tIns\tTotal"
  count_event_types "${myDir}/isec/0002.vcf" "Shared"
  count_event_types "${myDir}/isec/0000.vcf" "$name3"
  count_event_types "${myDir}/isec/0001.vcf" "$name4"

  $bcftools view ${myDir}/isec/0002.vcf -Oz -o ${myDir}/shared_sv_34.vcf.gz
  $bcftools index ${myDir}/shared_sv_34.vcf.gz

  echo "shared_12 vs shared_34"
  $bcftools isec -f "PASS" ${myDir}/shared_sv_12.vcf.gz ${myDir}/shared_sv_34.vcf.gz -p isec
  echo -e "Name\tDel\tSgl\tBnd\tInv\tDup\tIns\tTotal"
  count_event_types "${myDir}/isec/0002.vcf" "Shared"
  count_event_types "${myDir}/isec/0000.vcf" "$name3"
  count_event_types "${myDir}/isec/0001.vcf" "$name4"
fi

rm shared_*

#######################################################################################
#######################################################################################
#######################################################################################

echo -e "\nGermline SVs"
echo "${name1} vs ${name2}"
$bcftools isec -f "PASS" ${file1}/purple/*.purple.sv.germline.vcf.gz ${file2}/purple/*.purple.sv.germline.vcf.gz -p isec

echo -e "Name\tDel\tSgl\tBnd\tInv\tDup\tIns\tTotal"
count_event_types "${myDir}/isec/0002.vcf" "Shared"
count_event_types "${myDir}/isec/0000.vcf" "$name1"
count_event_types "${myDir}/isec/0001.vcf" "$name2"

# Generate and index the shared VCF
$bcftools view ${myDir}/isec/0002.vcf -Oz -o ${myDir}/shared_sv_12.vcf.gz
$bcftools index ${myDir}/shared_sv_12.vcf.gz

if [ ${runMode} -eq 4 ]; then
  # Process additional comparisons
  $bcftools isec -f "PASS" ${file3}/purple/*.purple.sv.germline.vcf.gz ${file4}/purple/*.purple.sv.germline.vcf.gz -p isec
  echo -e "Name\tDel\tSgl\tBnd\tInv\tDup\tIns\tTotal"
  count_event_types "${myDir}/isec/0002.vcf" "Shared"
  count_event_types "${myDir}/isec/0000.vcf" "$name3"
  count_event_types "${myDir}/isec/0001.vcf" "$name4"

  $bcftools view ${myDir}/isec/0002.vcf -Oz -o ${myDir}/shared_sv_34.vcf.gz
  $bcftools index ${myDir}/shared_sv_34.vcf.gz

  echo "shared_12 vs shared_34"
  $bcftools isec -f "PASS" ${myDir}/shared_sv_12.vcf.gz ${myDir}/shared_sv_34.vcf.gz -p isec
  echo -e "Name\tDel\tSgl\tBnd\tInv\tDup\tIns\tTotal"
  count_event_types "${myDir}/isec/0002.vcf" "Shared"
  count_event_types "${myDir}/isec/0000.vcf" "$name3"
  count_event_types "${myDir}/isec/0001.vcf" "$name4"
fi

rm shared_*

#######################################################################################
#######################################################################################
#######################################################################################

#identify uniques and shared 1&2 + plots
echo -e "\nGermline Variants"
$bcftools isec -f "PASS" ${file1}/purple/*.purple.germline.vcf.gz ${file2}/purple/*.purple.germline.vcf.gz -p isec
refonly=$(cat ${myDir}/isec/0000.vcf | grep -v "^#" | wc -l)  #refonly
shared=$(cat ${myDir}/isec/0002.vcf | grep -v "^#" | wc -l)  #shared_refdata
newonly=$(cat ${myDir}/isec/0001.vcf | grep -v "^#" | wc -l)  #newonly
echo -e "${name1}_only\tShared\t${name2}_only"
echo -e "${refonly}\t${shared}\t${newonly}"
#cat ~/coloOldvNew/isec/0003.vcf | grep -v "^#" | wc -l  #shared_newdata
$bcftools view ${myDir}/isec/0002.vcf -Oz -o ${myDir}/shared_gl_12.vcf.gz
$bcftools index ${myDir}/shared_gl_12.vcf.gz
./plotter.R $name1 $name2 GLVars ${isecDir} #> /dev/null 2>&1

if [ ${runMode} -eq 4 ]; then
  #identify uniques and shared 3&4 + plots
  $bcftools isec -f "PASS" ${file3}/purple/*.purple.germline.vcf.gz ${file4}/purple/*.purple.germline.vcf.gz -p isec
  refonly=$(cat ${myDir}/isec/0000.vcf | grep -v "^#" | wc -l)  #refonly
  shared=$(cat ${myDir}/isec/0002.vcf | grep -v "^#" | wc -l)  #shared_refdata
  newonly=$(cat ${myDir}/isec/0001.vcf | grep -v "^#" | wc -l)  #newonly
  echo -e "${name3}_only\tShared\t${name4}_only"
  echo -e "${refonly}\t${shared}\t${newonly}"
  #cat ~/coloOldvNew/isec/0003.vcf | grep -v "^#" | wc -l  #shared_newdata
  $bcftools view ${myDir}/isec/0002.vcf -Oz -o ${myDir}/shared_gl_34.vcf.gz
  $bcftools index ${myDir}/shared_gl_34.vcf.gz
  ./plotter.R $name3 $name4 GLVars ${isecDir} > /dev/null 2>&1

  #identify batch_specific_uniques and total_shareds
  $bcftools isec -f "PASS" ${myDir}/shared_gl_12.vcf.gz ${myDir}/shared_gl_34.vcf.gz -p isec
  refonly=$(cat ${myDir}/isec/0000.vcf | grep -v "^#" | wc -l)  #refonly
  shared=$(cat ${myDir}/isec/0002.vcf | grep -v "^#" | wc -l)  #shared_refdata
  newonly=$(cat ${myDir}/isec/0001.vcf | grep -v "^#" | wc -l)  #newonly
  echo -e "${name1}/${name2}_only\tShared all\t${name3}/${name4}_only"
  echo -e "${refonly}\t${shared}\t${newonly}"

  #extract all final Vars
  $bcftools view ${myDir}/isec/0000.vcf -Oz -o ${myDir}/only_12.vcf.gz
  $bcftools index ${myDir}/only_12.vcf.gz
  $bcftools view ${myDir}/isec/0001.vcf -Oz -o ${myDir}/only_34.vcf.gz
  $bcftools index ${myDir}/only_34.vcf.gz
  $bcftools view ${myDir}/isec/0002.vcf -Oz -o ${myDir}/shared_all.vcf.gz
  $bcftools index ${myDir}/shared_all.vcf.gz

  #extract old-only regions from old-files
  $bcftools isec -f "PASS" ${myDir}/only_12.vcf.gz ${file1}/purple/*.purple.germline.vcf.gz -p isec
  $bcftools view ${myDir}/isec/0003.vcf -Oz -o ${myDir}/only_12_1.vcf.gz
  $bcftools index ${myDir}/only_12_1.vcf.gz
  $bcftools isec -f "PASS" ${myDir}/only_12.vcf.gz ${file2}/purple/*.purple.germline.vcf.gz -p isec
  $bcftools view ${myDir}/isec/0003.vcf -Oz -o ${myDir}/only_12_2.vcf.gz
  $bcftools index ${myDir}/only_12_2.vcf.gz

  #extract new-only regions from new-files
  $bcftools isec -f "PASS" ${myDir}/only_34.vcf.gz ${file3}/purple/*.purple.germline.vcf.gz -p isec
  $bcftools view ${myDir}/isec/0003.vcf -Oz -o ${myDir}/only_34_3.vcf.gz
  $bcftools index ${myDir}/only_34_3.vcf.gz
  $bcftools isec -f "PASS" ${myDir}/only_34.vcf.gz ${file4}/purple/*.purple.germline.vcf.gz -p isec
  $bcftools view ${myDir}/isec/0003.vcf -Oz -o ${myDir}/only_34_4.vcf.gz
  $bcftools index ${myDir}/only_34_4.vcf.gz

  #extract shared lists from all
  $bcftools isec -f "PASS" ${myDir}/shared_all.vcf.gz ${file1}/purple/*.purple.germline.vcf.gz -p isec
  $bcftools view ${myDir}/isec/0003.vcf -Oz -o ${myDir}/shared_all_1.vcf.gz
  $bcftools index ${myDir}/shared_all_1.vcf.gz
  $bcftools isec -f "PASS" ${myDir}/shared_all.vcf.gz ${file2}/purple/*.purple.germline.vcf.gz -p isec
  $bcftools view ${myDir}/isec/0003.vcf -Oz -o ${myDir}/shared_all_2.vcf.gz
  $bcftools index ${myDir}/shared_all_2.vcf.gz
  $bcftools isec -f "PASS" ${myDir}/shared_all.vcf.gz ${file3}/purple/*.purple.germline.vcf.gz -p isec
  $bcftools view ${myDir}/isec/0003.vcf -Oz -o ${myDir}/shared_all_3.vcf.gz
  $bcftools index ${myDir}/shared_all_3.vcf.gz
  $bcftools isec -f "PASS" ${myDir}/shared_all.vcf.gz ${file4}/purple/*.purple.germline.vcf.gz -p isec
  $bcftools view ${myDir}/isec/0003.vcf -Oz -o ${myDir}/shared_all_4.vcf.gz
  $bcftools index ${myDir}/shared_all_4.vcf.gz

  ./shrdPlotter.R GLVars ${myDir} #> /dev/null 2>&1
fi

mv ./*.png expPlots 2>/dev/null
mv ./isec/*.png expPlots 2>/dev/null
rm shared_* 2>/dev/null
rm only_* 2>/dev/null


#######################################################################################
#######################################################################################

#identify uniques and shared 1&2 + plots
echo -e "\nSomatic Variants"
$bcftools isec -f "PASS" ${file1}/purple/*.purple.somatic.vcf.gz ${file2}/purple/*.purple.somatic.vcf.gz -p isec
refonly=$(cat ${myDir}/isec/0000.vcf | grep -v "^#" | wc -l)  #refonly
shared=$(cat ${myDir}/isec/0002.vcf | grep -v "^#" | wc -l)  #shared_refdata
newonly=$(cat ${myDir}/isec/0001.vcf | grep -v "^#" | wc -l)  #newonly
echo -e "${name1}_only\tShared\t${name2}_only"
echo -e "${refonly}\t${shared}\t${newonly}"
#cat ~/coloOldvNew/isec/0003.vcf | grep -v "^#" | wc -l  #shared_newdata
$bcftools view ${myDir}/isec/0002.vcf -Oz -o ${myDir}/shared_som_12.vcf.gz
$bcftools index ${myDir}/shared_som_12.vcf.gz
./plotter.R $name1 $name2 SomVars ${isecDir} > /dev/null 2>&1
./triNucPlot.R $name1 $name2 SomVars $hgVersion ${isecDir} > /dev/null 2>&1

if [ ${runMode} -eq 4 ]; then
  #identify uniques and shared 3&4 + plots
  $bcftools isec -f "PASS" ${file3}/purple/*.purple.somatic.vcf.gz ${file4}/purple/*.purple.somatic.vcf.gz -p isec
  refonly=$(cat ${myDir}/isec/0000.vcf | grep -v "^#" | wc -l)  #refonly
  shared=$(cat ${myDir}/isec/0002.vcf | grep -v "^#" | wc -l)  #shared_refdata
  newonly=$(cat ${myDir}/isec/0001.vcf | grep -v "^#" | wc -l)  #newonly
  echo -e "${name3}_only\tShared\t${name4}_only"
  echo -e "${refonly}\t${shared}\t${newonly}"
  #cat ~/coloOldvNew/isec/0003.vcf | grep -v "^#" | wc -l  #shared_newdata
  $bcftools view ${myDir}/isec/0002.vcf -Oz -o ${myDir}/shared_som_34.vcf.gz
  $bcftools index ${myDir}/shared_som_34.vcf.gz
  ./plotter.R $name3 $name4 SomVars ${isecDir} #> /dev/null 2>&1
  ./triNucPlot.R $name3 $name4 SomVars $hgVersion ${isecDir} > /dev/null 2>&1

  #identify batch_specific_uniques and total_shareds
  $bcftools isec -f "PASS" ${myDir}/shared_som_12.vcf.gz ${myDir}/shared_som_34.vcf.gz -p isec
  refonly=$(cat ${myDir}/isec/0000.vcf | grep -v "^#" | wc -l)  #refonly
  shared=$(cat ${myDir}/isec/0002.vcf | grep -v "^#" | wc -l)  #shared_refdata
  newonly=$(cat ${myDir}/isec/0001.vcf | grep -v "^#" | wc -l)  #newonly
  echo -e "${name1}/${name2}_only\tShared all\t${name3}/${name4}_only"
  echo -e "${refonly}\t${shared}\t${newonly}"
  ./triNucPlot.R Old_shrd New_shrd SomVars $hgVersion ${isecDir} #> /dev/null 2>&1

  #extract all final Vars
  $bcftools view ${myDir}/isec/0000.vcf -Oz -o ${myDir}/only_12.vcf.gz
  $bcftools index ${myDir}/only_12.vcf.gz
  $bcftools view ${myDir}/isec/0001.vcf -Oz -o ${myDir}/only_34.vcf.gz
  $bcftools index ${myDir}/only_34.vcf.gz
  $bcftools view ${myDir}/isec/0002.vcf -Oz -o ${myDir}/shared_all.vcf.gz
  $bcftools index ${myDir}/shared_all.vcf.gz

  #extract old-only regions from old-files
  $bcftools isec -f "PASS" ${myDir}/only_12.vcf.gz ${file1}/purple/*.purple.somatic.vcf.gz -p isec
  $bcftools view ${myDir}/isec/0003.vcf -Oz -o ${myDir}/only_12_1.vcf.gz
  $bcftools index ${myDir}/only_12_1.vcf.gz
  $bcftools isec -f "PASS" ${myDir}/only_12.vcf.gz ${file2}/purple/*.purple.somatic.vcf.gz -p isec
  $bcftools view ${myDir}/isec/0003.vcf -Oz -o ${myDir}/only_12_2.vcf.gz
  $bcftools index ${myDir}/only_12_2.vcf.gz

  #extract new-only regions from new-files
  $bcftools isec -f "PASS" ${myDir}/only_34.vcf.gz ${file3}/purple/*.purple.somatic.vcf.gz -p isec
  $bcftools view ${myDir}/isec/0003.vcf -Oz -o ${myDir}/only_34_3.vcf.gz
  $bcftools index ${myDir}/only_34_3.vcf.gz
  $bcftools isec -f "PASS" ${myDir}/only_34.vcf.gz ${file4}/purple/*.purple.somatic.vcf.gz -p isec
  $bcftools view ${myDir}/isec/0003.vcf -Oz -o ${myDir}/only_34_4.vcf.gz
  $bcftools index ${myDir}/only_34_4.vcf.gz

  #extract shared lists from all
  $bcftools isec -f "PASS" ${myDir}/shared_all.vcf.gz ${file1}/purple/*.purple.somatic.vcf.gz -p isec
  $bcftools view ${myDir}/isec/0003.vcf -Oz -o ${myDir}/shared_all_1.vcf.gz
  $bcftools index ${myDir}/shared_all_1.vcf.gz
  $bcftools isec -f "PASS" ${myDir}/shared_all.vcf.gz ${file2}/purple/*.purple.somatic.vcf.gz -p isec
  $bcftools view ${myDir}/isec/0003.vcf -Oz -o ${myDir}/shared_all_2.vcf.gz
  $bcftools index ${myDir}/shared_all_2.vcf.gz
  $bcftools isec -f "PASS" ${myDir}/shared_all.vcf.gz ${file3}/purple/*.purple.somatic.vcf.gz -p isec
  $bcftools view ${myDir}/isec/0003.vcf -Oz -o ${myDir}/shared_all_3.vcf.gz
  $bcftools index ${myDir}/shared_all_3.vcf.gz
  $bcftools isec -f "PASS" ${myDir}/shared_all.vcf.gz ${file4}/purple/*.purple.somatic.vcf.gz -p isec
  $bcftools view ${myDir}/isec/0003.vcf -Oz -o ${myDir}/shared_all_4.vcf.gz
  $bcftools index ${myDir}/shared_all_4.vcf.gz

  ./shrdPlotter.R SomVars ${myDir} > /dev/null 2>&1
  ./shrdNucPlotter.R $name1 $name2 $name3 $name4 SomVars $hgVersion ${myDir} #> /dev/null 2>&1
fi

mv ./*.png expPlots 2>/dev/null
mv ./isec/*.png expPlots 2>/dev/null
rm shared_* 2>/dev/null
rm only_* 2>/dev/null


#######################################################################################
#######################################################################################
#######################################################################################

echo -e "\nSomatic Drivers"
$bcftools filter -o shared_file1_somdr.vcf.gz -O z -i 'INFO/REPORTED="0"' ${file1}/purple/*.purple.somatic.vcf.gz
$bcftools index shared_file1_somdr.vcf.gz
$bcftools filter -o shared_file2_somdr.vcf.gz -O z -i 'INFO/REPORTED="0"' ${file2}/purple/*.purple.somatic.vcf.gz
$bcftools index shared_file2_somdr.vcf.gz
$bcftools isec shared_file1_somdr.vcf.gz shared_file2_somdr.vcf.gz -p isec
SomDrR=$(cat ${myDir}/isec/0000.vcf | grep -v "^#" | wc -l)  #refonly
SomDrN=$(cat ${myDir}/isec/0001.vcf | grep -v "^#" | wc -l)  #newonly
SomDrS=$(cat ${myDir}/isec/0002.vcf | grep -v "^#" | wc -l)  #shared_refdata
echo -e "${name1}_only\tShared\t${name2}_only"
echo -e "${SomDrR}\t${SomDrS}\t${SomDrN}"
#cat ~/coloOldvNew/isec/0003.vcf | grep -v "^#" | wc -l  #shared_newdata
$bcftools view ${myDir}/isec/0002.vcf -Oz -o ${myDir}/shared_somdri_12.vcf.gz
$bcftools index ${myDir}/shared_somdri_12.vcf.gz

if [ ${runMode} -eq 4 ]; then
  $bcftools filter -o shared_file3_somdr.vcf.gz -O z -i 'INFO/REPORTED="0"' ${file3}/purple/*.purple.somatic.vcf.gz
  $bcftools index shared_file3_somdr.vcf.gz
  $bcftools filter -o shared_file4_somdr.vcf.gz -O z -i 'INFO/REPORTED="0"' ${file4}/purple/*.purple.somatic.vcf.gz
  $bcftools index shared_file4_somdr.vcf.gz
  $bcftools isec shared_file3_somdr.vcf.gz shared_file4_somdr.vcf.gz -p isec
  SomDrR=$(cat ${myDir}/isec/0000.vcf | grep -v "^#" | wc -l)  #refonly
  SomDrN=$(cat ${myDir}/isec/0001.vcf | grep -v "^#" | wc -l)  #newonly
  SomDrS=$(cat ${myDir}/isec/0002.vcf | grep -v "^#" | wc -l)  #shared_refdata
  echo -e "${name1}_only\tShared\t${name2}_only"
  echo -e "${SomDrR}\t${SomDrS}\t${SomDrN}"
  #cat ~/coloOldvNew/isec/0003.vcf | grep -v "^#" | wc -l  #shared_newdata
  $bcftools view ${myDir}/isec/0002.vcf -Oz -o ${myDir}/shared_somdri_34.vcf.gz
  $bcftools index ${myDir}/shared_somdri_34.vcf.gz

  $bcftools isec -f "PASS" ${myDir}/shared_somdri_12.vcf.gz ${myDir}/shared_somdri_34.vcf.gz -p isec
  SomDrR=$(cat ${myDir}/isec/0000.vcf | grep -v "^#" | wc -l)  #refonly
  SomDrN=$(cat ${myDir}/isec/0001.vcf | grep -v "^#" | wc -l)  #newonly
  SomDrS=$(cat ${myDir}/isec/0002.vcf | grep -v "^#" | wc -l)  #shared_refdata
  echo -e "${name1}/${name2}_only\tShared all\t${name3}/${name4}_only"
  echo -e "${SomDrR}\t${SomDrS}\t${SomDrN}"

fi

mv ./*.png expPlots 2>/dev/null
mv ./isec/*.png expPlots 2>/dev/null
rm shared_* 2>/dev/null
rm only_* 2>/dev/null

#######################################################################################
#######################################################################################
#######################################################################################

echo -e "\nCN Drivers"
cat ${file1}/purple/*.driver.catalog.somatic.tsv | grep -v MUTATION | tail -n +2 | awk '{ print $1, $2, $3, $4, $6 }' | sort > cnd1.txt
cat ${file2}/purple/*.driver.catalog.somatic.tsv | grep -v MUTATION | tail -n +2 | awk '{ print $1, $2, $3, $4, $6 }' | sort > cnd2.txt

# Find unique and shared entries using comm
CNDr=$(comm -23 cnd1.txt cnd2.txt | wc -l)  # Unique in CND1
CNDn=$(comm -13 cnd1.txt cnd2.txt | wc -l)  # Unique in CND2
CNDs=$(comm -12 cnd1.txt cnd2.txt | wc -l)  # Shared entries
comm -12 cnd1.txt cnd2.txt > shared_12.txt

echo -e "${name1}\tShared\t ${name2}"
echo -e "${CNDr}\t${CNDs}\t${CNDn}"
rm cnd1.txt cnd2.txt

if [ ${runMode} -eq 4 ]; then
  cat ${file3}/purple/*.driver.catalog.somatic.tsv | grep -v MUTATION | tail -n +2 | awk '{ print $1, $2, $3, $4, $6 }' | sort > cnd3.txt
  cat ${file4}/purple/*.driver.catalog.somatic.tsv | grep -v MUTATION | tail -n +2 | awk '{ print $1, $2, $3, $4, $6 }' | sort > cnd4.txt

  # Find unique and shared entries using comm
  CNDr=$(comm -23 cnd3.txt cnd4.txt | wc -l)  # Unique in CND1
  CNDn=$(comm -13 cnd3.txt cnd4.txt | wc -l)  # Unique in CND2
  CNDs=$(comm -12 cnd3.txt cnd4.txt | wc -l)  # Shared entries
  comm -12 cnd3.txt cnd4.txt > shared_34.txt

  echo -e "${name3}\tShared\t ${name4}"
  echo -e "${CNDr}\t${CNDs}\t${CNDn}"
  rm cnd3.txt cnd4.txt

  # Find unique and shared entries using comm
  CNDr=$(comm -23 shared_12.txt shared_34.txt | wc -l)  # Unique in CND1
  CNDn=$(comm -13 shared_12.txt shared_34.txt | wc -l)  # Unique in CND2
  CNDs=$(comm -12 shared_12.txt shared_34.txt | wc -l)  # Shared entries
  echo -e "${name1}/${name2}_only\tShared all\t${name3}/${name4}_only"
  echo -e "${CNDr}\t${CNDs}\t${CNDn}"

fi

mv ./*.png expPlots 2>/dev/null
mv ./isec/*.png expPlots 2>/dev/null
rm shared_* 2>/dev/null
rm only_* 2>/dev/null


#######################################################################################
#######################################################################################
#######################################################################################
