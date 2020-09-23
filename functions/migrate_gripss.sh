#!/usr/bin/env bash

run_dir=$1

input_dir=${run_dir}/gridss
output_dir=${run_dir}/gripss

if [[ ! -d ${input_dir} ]]; then
  echo "Input dir does not exist ${run_dir}"
#  exit 0
fi

if [[ -d ${output_dir} ]]; then
  echo "Already migrated ${run_dir}"
#  exit 0
fi

mkdir ${output_dir}

gridss_filtered=$(find ${input_dir} -name *.gridss.somatic.filtered.vcf.gz)
gridss_filtered_index=${gridss_filtered}.tbi

gripss_filtered=${gridss_filtered//gridss/gripss}
gripss_filtered_index=${gridss_filtered_index//gridss/gripss}

mv ${gridss_filtered} ${gripss_filtered}
mv ${gridss_filtered_index} ${gripss_filtered_index}

gridss_somatic=$(find ${input_dir} -name *.gridss.somatic.vcf.gz)
gridss_somatic_index=${gridss_somatic}.tbi

gripss_somatic=${gridss_somatic//gridss/gripss}
gripss_somatic_index=${gridss_somatic_index//gridss/gripss}

mv ${gridss_somatic} ${gripss_somatic}
mv ${gridss_somatic_index} ${gripss_somatic_index}

if [ ! "$(ls -A ${input_dir})" ]; then
    rm -rf ${input_dir}
fi

