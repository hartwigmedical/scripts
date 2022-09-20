#!/usr/bin/env bash
# Purpose: Confirm SV nt-sequences deduced from .purple.sv.vcf file
# Author: Teoman Deger
# Version 0.1, 20-09-2022
# -----------------------

INPUT=~/data/SV-Test.csv
IFS=','
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
#while read line
while read -r SVnum sample SVTYPE seqnameso starto endo widtho strando REFo ALTo seqnamesh starth endh widthh strandh REFh ALTh svSequence Extension Ligation
do
  Ext_noQuot=${Extension:1:-1}
  Lig_noQuot=${Ligation:1:-1}
  #/opt/tools/samtools/1.14/samtools view -@ 25 COLO829v003T.bam ${seqnameso}:${starto}-${endo} ${seqnamesh}:${starth}-${endh} | grep -c ${Ext_noQuot}[[:alpha:]]\{0,0\}${Lig_noQuot}
  /opt/tools/samtools/1.14/samtools view -@ 25 COLO829v003T.bam ${seqnameso}:${starto}-${endo} ${seqnamesh}:${starth}-${endh} | grep -c ${Ext_noQuot}
done < <(tail -n +2 $INPUT)