#!/usr/bin/env bash
# -----------------------

#list of cram-output, aligned and sorted
myList=("gs://ultima-ppmseq-analysis/analysis/fr36584653/sorter/sorter_output/sample_output-/sample_output-sample/sample.cram"
"gs://ultima-ppmseq-analysis/analysis/fr34221950/sorter/sorter_output/sample_output-/sample_output-sample/sample.cram"
"gs://ultima-ppmseq-analysis/analysis/fr34240372/sorter/sorter_output/sample_output-/sample_output-sample/sample.cram")

#directories
workDir="/home/tdeger/GATK-try"
outDir="/home/tdeger/GATK-try/outputSNPs"

#sanitycheck directories
if [ ! -d ${outDir} ]; then
  echo Output directory \"${outDir}\" doesnt exist, creating
  mkdir -p ${outDir}; fi

for num in ${!myList[@]}
do
  #extract filename, extention and name without extention
  file=$(echo ${myList[$num]} | awk -F '/' '{ print $NF }' )
  type=$(echo ${myList[$num]} | awk -F '.' '{ print $NF }' )
  #name=$(echo ${file} | awk -F '.' '{ print $1 }' ) #use this one for files if they have a unique name e.g. COREDB0123456T.cram
  name=$(echo ${myList[$num]} | awk -F '/' '{ print $5 }' )
  echo "Filename: ${name}"
  echo "Filetype: ${type}"
  #downloads file to local
  gsutil  -o 'GSUtil:parallel_thread_count=1' -o GSUtil:sliced_object_download_max_components=$(nproc) -m cp ${myList[$num]} ${workDir}
  if [[ ${type} == cram ]]; then
    gsutil  -o 'GSUtil:parallel_thread_count=1' -o GSUtil:sliced_object_download_max_components=$(nproc) -m cp ${myList[$num]}.crai ${workDir}; fi
  if [[ ${type} == bam ]]; then
    gsutil  -o 'GSUtil:parallel_thread_count=1' -o GSUtil:sliced_object_download_max_components=$(nproc) -m cp ${myList[$num]}.bai ${workDir}; fi
  #runs docker image to do actually SNP genotyping
  docker run \
    -v ${workDir}:${workDir} \
    --entrypoint /gatk/gatk broadinstitute/gatk:4.6.0.0 HaplotypeCaller \
    -L ${workDir}/31SNPsmMIP_design_hg38.vcf \
    -R ${workDir}/hg38ult/Homo_sapiens_assembly38.fasta \
    -I ${workDir}/${file} \
    -O ${workDir}/${name}.vcf \
    -ERC BP_RESOLUTION
  #moves outputs and removes downloaded files
  mv ${name}.vcf ${outDir}
  mv ${name}.vcf.idx ${outDir}
  rm ${workDir}/${file}
  if [[ ${type} == cram ]]; then
    rm ${workDir}/${file}.crai; fi
  if [[ ${type} == bam ]]; then
    rm ${workDir}/${file}.bai; fi
done

#sanitycheck SNPcheck output
for file in *.vcf
do
  cat ${file} | tail -n 31 | awk -F ' ' '{print $10}' | awk -F ':' '{print $2}' > tempVCF.txt
  #bcftools query -f'[%AD\n]' ${file} > tempVCF.txt
  passed=0
  while read p; do
    num1=$(echo "$p" | awk -F ',' '{print $1}')
    num2=$(echo "$p" | awk -F ',' '{print $2}')
    num3=$(echo "$p" | awk -F ',' '{print $3}')
    if [ -z "$num3" ]
    then
      count=$(echo "${num1}+${num2}" | bc) #if output has 2 counts e.g. 23,15
    else
      count=$(echo "${num1}+${num2}+${num3}" | bc) #if output has 3 counts e.g. 16,13,1
    fi
    if (( ${count} > 9 )) #if more than 10x coverage in a location, that location is passed (just as in AMBER)
    then
      passed=$(($passed+1))
    fi
  done <tempVCF.txt
  if (( ${passed} > 17 )) #if at least 18/31 SNP locations are passed, that sample is succesfull
  then
    touch ${file}_SUCCESS
  else
    touch ${file}_WARNING
  fi
  rm tempVCF.txt
done

exit