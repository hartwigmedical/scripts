#!/bin/bash
#
#

rna_root=/data2/rna
tools_root=$rna_root/tools
ref_root=$rna_root/ref

star_version=2.7.3a
star_root=$tools_root/STAR-${star_version}
arriba_version=1.1.0
arriba_root=$tools_root/arriba_v${arriba_version}/
starfusion_version=1.8.1
starfusion_root=$tools_root/STAR-Fusion-v${starfusion_version}

# Add tools to path
export PATH=$arriba_root:$star_root/bin/Linux_x86_64:$starfusion_root:$PATH

##################################
# Once-off installation/setup
# Install tools
mkdir -p $rna_root $tools_root $ref_root
cd $tools_root
if [[ ! -d $arriba_root ]] ; then
  wget https://github.com/suhrig/arriba/releases/download/v${arriba_version}/arriba_v${arriba_version}.tar.gz
  tar -xzf arriba_v${arriba_version}.tar.gz
  rm arriba_v${arriba_version}.tar.gz
fi
if [[ ! -d $starfusion_root ]] ; then
  wget https://github.com/STAR-Fusion/STAR-Fusion/releases/download/STAR-Fusion-v${starfusion_version}/STAR-Fusion-${starfusion_version}.FULL.tar.gz
  tar -xzf STAR-Fusion-v${starfusion_version}.FULL.tar.gz
  rm STAR-Fusion-v${starfusion_version}.FULL.tar.gz
  # https://github.com/STAR-Fusion/STAR-Fusion/wiki/installing-star-fusion#tools-required
  #perl -MCPAN -e shell
  #install DB_File
  #install URI::Escape
  #install Set::IntervalTree
  #install Carp::Assert
  #install JSON::XS
  #install PerlIO::gzip
fi
if [[ ! -d $star_root ]] ; then
  wget https://github.com/alexdobin/STAR/archive/${star_version}.tar.gz
  tar -xzf ${star_version}.tar.gz
  rm ${star_version}.tar.gz
fi
# Set up reference and indexes
# arriba
if [[ ! -d $ref_root/hs37d5_GENCODE19 ]] ; then
  mkdir $ref_root/hs37d5_GENCODE19
  cd $ref_root/hs37d5_GENCODE19
  $arriba_root/download_references.sh hs37d5+GENCODE19
fi
# star fusion
if [[ ! -d $ref_root/GRCh37_gencode_v19_CTAT_lib_Oct012019.plug-n-play ]] ; then
  wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh37_gencode_v19_CTAT_lib_Oct012019.plug-n-play.tar.gz
  tar -xzf GRCh37_gencode_v19_CTAT_lib_Oct012019.plug-n-play.tar.gz
  rm GRCh37_gencode_v19_CTAT_lib_Oct012019.plug-n-play.tar.gz
fi



##################################
# Per sample processing
input_dir=$1
# Per sample inputs
if [[ ! -d $input_dir ]] ; then
  echo "Input directory not found. Please specify the directory with the RNA-Seq sequencing data on the command line" 1>&2
  #exit 1
fi

cd $input_dir
r1=$(find $input_dir/*_R1_001.fastq.gz | tr '\n' ',')
r2=${r1//_R1_001.fastq.gz/_R2_001.fastq.gz}
# arriba (and star)
( $arriba_root/run_arriba.sh \
  $ref_root/hs37d5_GENCODE19/STAR_index_hs37d5_GENCODE19/ \
  $ref_root/hs37d5_GENCODE19/GENCODE19.gtf \
  $ref_root/hs37d5_GENCODE19/hs37d5.fa \
  $arriba_root/database/blacklist_hg19_hs37d5_GRCh37_2018-11-04.tsv.gz \
  $r1 \
  $r2 \
  $(nproc) )

# STAR-Fusion
# TODO: do we wantt to run FusionInspector as well?
( STAR-Fusion \
  --genome_lib_dir $ref_root/GRCh37_gencode_v19_CTAT_lib_Oct012019 \
  --left_fq $r1 \
  --right_fq $r2 \
  --output_dir $input_dir/star_fusion_outdir )

  


