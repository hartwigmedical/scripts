#!/bin/bash
#input_vcf=$1
#output=$2
if [[ ! -f input.list ]] ; then
        echo "Generating input list"
        find /data/cpct/runs/ -name '*.gridss.vcf.gz' > input.list
fi
for input_vcf in $(cat input.list) ; do
        output=/data/experiments/germline_tis/candidates/$(basename $input_vcf).ti.vcf
        log=/data/experiments/germline_tis/logs/$(basename $input_vcf).log
        filtered_vcf=${output}.tmp.filtered.vcf
        gunzip -c $input_vcf \
        | awk ' { if (length($0) >= 4000) { gsub(":0.00:", ":0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000:")} ; print $0  } ' \
        | awk '$1 ~ /^#/ || $6 > 350 && $5 ~ /[\[\]]/' > $filtered_vcf
        Rscript find_germline_TIs.R --input $filtered_vcf --output $output 2>&1 > $log
        rm $filtered_vcf
done
