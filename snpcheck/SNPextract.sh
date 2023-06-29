#!/usr/bin/env bash
# Author: Teoman Deger
# Purpose: make bam from MIP-data and extract wanted data
# Version 0.2, 03-02-2023
# -----------------------

bwadir='/data/tools/bwa/0.7.17/'
samtdir='/data/tools/samtools/1.14/'
source locate_files || exit 1
source message_functions || exit 1

print_usage(){
    echo ""
    echo "Descr: Extract SNP-types from .fastq or .bam obtained from smMIP-experiments"
    echo "Usage: $(basename "$0") -I /pathTo/inputDir  -L /pathTo/smMIP-list -R SNP/SVs (optional: -O /pathTo/outputDir  -X ExportFile-name"
    echo "------------------------------------------------------------"
    exit 1
}

#option settings
outputDir=./
while getopts ':I:L:R:O:X' flag; do
    case "${flag}" in
        I) inDir=${OPTARG} ;;
        L) smMIPfile=${OPTARG} ;;
        R) runMode=${OPTARG} ;;
        O) outputDir=${OPTARG} ;;
        X) expFile=${OPTARG} ;;
        *) print_usage >&2
        exit 1 ;;
    esac
done

#sanity checks
if [[ -z "${inDir}" || -z "${smMIPfile}" || -z "${runMode}" ]]; then
 print_usage; fi

if [[ -z "${expFile}" ]]; then
 expFile="temp"; fi

if [[ ${runMode} != SNP ]] && [[ ${runMode} != SVs ]]; then
  echo "Run-mode required: either '-R SNP' or '-R SVs'"
  print_usage; fi

if [[ $# -eq 0 ]] ; then
    print_usage; fi

if [ ! -d ${outputDir} ]; then
  echo Output directory \"${outputDir}\" doesnt exist, creating
  mkdir -p ${outputDir}; fi

if [[ ! -f ${smMIPfile} ]]; then
 echo "List of smMIP sequences is not a file, terminating.";
 exit 1; fi

if [[ ! -d ${inDir} ]]; then
 echo "Given input-dir is not a directory";
 exit 1; fi

if [[  -z ${inDir} ]]; then
  echo "Given input-dir contains no files"
  exit 1; fi

if [[ -n ${inDir} ]] && [[ $(ls -lR ${inDir}*.fastq.gz | wc -l) != 0 ]]; then
  for i in *1_R1_001.fastq.gz
  do
     header=$(zcat $i | head -n 1)
     id=$(echo $header | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g')
     lb=${i:0:11}
     sm=$(echo $header | head -n 1 | grep -Eo "[ATGCN]+$")
     pu=$(echo $header | head -n 1 | cut -f 3 -d":" | sed 's/@//' | sed 's/:/_/g')
     rg="@RG\tID:$id\tLB:$lb\tPL:ILLUMINA\tPU:$pu\tSM:$sm"
     j=${i::-17}
     if [[ ! -f ${lb}_sorted_merged.bam ]]; then
     echo combining ${lb} files
     ${bwadir}/bwa mem -R $rg -t 8 ./ref-snMIP/reference-smMIP.fa ${j}1_R1_001.fastq.gz ${j}1_R2_001.fastq.gz | ${samtdir}/samtools sort -@ 8 -o ${lb}_L1_sorted.bam -
     ${bwadir}/bwa mem -R $rg -t 8 ./ref-snMIP/reference-smMIP.fa ${j}2_R1_001.fastq.gz ${j}2_R2_001.fastq.gz | ${samtdir}/samtools sort -@ 8 -o ${lb}_L2_sorted.bam -
     #${bwadir}/bwa mem -R $rg -t 8 ./ref-snMIP/reference-smMIP.fa ${j}3_R1_001.fastq.gz ${j}3_R2_001.fastq.gz | ${samtdir}/samtools sort -@ 8 -o ${lb}_L3_sorted.bam -
     #${bwadir}/bwa mem -R $rg -t 8 ./ref-snMIP/reference-smMIP.fa ${j}4_R1_001.fastq.gz ${j}4_R2_001.fastq.gz | ${samtdir}/samtools sort -@ 8 -o ${lb}_L4_sorted.bam -
     echo converting sam-file to sorted bam-file and merge
     ${samtdir}/samtools merge ${lb}_sorted_merged.bam ${lb}_L1_sorted.bam ${lb}_L2_sorted.bam #${lb}_L3_sorted.bam ${lb}_L4_sorted.bam
     ${samtdir}/samtools index -b -@ 4 ${lb}_sorted_merged.bam
     echo "---------------------------------"
     rm ${lb}_L*
     else
       echo "${naming}_sorted_merged.bam already exists, continuing with next"
     fi
  done
fi

echo "continuing"

#Sample ID,Plate Barcode,Gene Symbol,NCBI SNP Reference,Assay Name or ID,Allele 1 Call,Allele 2 Call
if [[ ${runMode} == SNP ]]; then
echo "SampleID, SNP-ID, A, C, T, G, noSNP" > ${expFile}.csv
for k in *.bam
do
  l=${k::-18}
  while read -r ID extSeq ligSeq; do
    if [[ ${extSeq} == "Extension_Sequence" ]]; then
      continue
    fi
    aCounts=$(${samtdir}/samtools view -@ 6 ${k} $ID | grep -c ${extSeq}A${ligSeq})
    cCounts=$(${samtdir}/samtools view -@ 6 ${k} $ID | grep -c ${extSeq}C${ligSeq})
    tCounts=$(${samtdir}/samtools view -@ 6 ${k} $ID | grep -c ${extSeq}T${ligSeq})
    gCounts=$(${samtdir}/samtools view -@ 6 ${k} $ID | grep -c ${extSeq}G${ligSeq})
    noSNP=$(${samtdir}/samtools view -@ 6 ${k} $ID | grep -c ${extSeq}${ligSeq})
    echo "${l}, ${ID}, ${aCounts}, ${cCounts}, ${tCounts}, ${gCounts}, ${noSNP}" >> ${expFile}.csv
  done < <(sed 's/\r$//' ${smMIPfile})
done
fi


if [[ ${runMode} == SVs ]]; then
echo "SampleID, SV-ID" > ${expFile}.csv
for k in *.bam
do
  l=${k::-18}
  while read -r ID extSeq ligSeq; do
    if [[ ${extSeq} == "Extension_Sequence" ]]; then
      continue
    fi
    SVcounts=$(${samtdir}/samtools view ${k} | grep -c ${extSeq}${ligSeq})
    echo "${l}, ${ID}, ${SVcounts}" >> ${expFile}.csv
  done < <(sed 's/\r$//' ${smMIPfile})
done
fi