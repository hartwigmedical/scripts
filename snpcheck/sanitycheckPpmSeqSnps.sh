#!/usr/bin/env bash
# -----------------------

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
      count=$(echo "${num1}+${num2}" | bc)
    else
      count=$(echo "${num1}+${num2}+${num3}" | bc)
    fi
    if (( ${count} > 9 ))
    then
      passed=$(($passed+1))
    fi
  done <tempVCF.txt
  if (( ${passed} > 17 ))
  then
    touch ${file}_SUCCESS
  else
    touch ${file}_WARNING
  fi
  rm tempVCF.txt
done

exit

