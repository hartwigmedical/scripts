#!/usr/bin/env bash
# Purpose: Confirm SV nt-sequences deduced from .purple.sv.vcf file
# Author: Teoman Deger
# Version 0.1, 20-09-2022
# -----------------------./

INPUT=~/data/SV-Test.csv
SAMTOOLSDIR="/opt/tools/samtools/1.14"
BAMLOC=~/data/COLO829v003T.bam
IFS=';'

[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read -r SVnum sample SVTYPE seqnameso starto endo widtho strando REFo ALTo seqnamesh starth endh widthh strandh REFh ALTh svSequence Extension Ligation
do
  Ext_noQuot=${Extension:1:-1}
  Lig_noQuot=${Ligation:1:-1}
  greptxt=${Ext_noQuot}[[:alpha:]]\\{0,0\\}${Lig_noQuot}
  echo SV-check ${SVnum}: Checking presence of ${SVTYPE} of ${seqnameso}:${starto} and ${seqnamesh}:${starth}
  FOUND=$(${SAMTOOLSDIR}/samtools view -@ 25 ${BAMLOC} ${seqnameso}:${starto}-${endo} ${seqnamesh}:${starth}-${endh} | grep -c ${greptxt})
  if [ ${FOUND} -ne 0 ]
  then
    echo ${FOUND} perfectly matching sequence reads found
  else
    echo No perfectly matching sequence reads found
    i=0
    while [ ${FOUND} -eq 0 ]
    do
     ((i++))
     if [[ $i -eq 50 ]]
     then
       break
     fi
     greptxt=${Ext_noQuot}[[:alpha:]]\\{0,${i}\\}${Lig_noQuot}
     FOUND=$(${SAMTOOLSDIR}/samtools view -@ 25 ${BAMLOC} ${seqnameso}:${starto}-${endo} ${seqnamesh}:${starth}-${endh} | grep -c ${greptxt})
    done
    echo ${FOUND} matching reads found within ${i}nt distance
  fi
  echo -e "\n"
done < <(tail -n +2 $INPUT)