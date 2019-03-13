#!/bin/bash
#
# Stand-alone chromaticHMF/GRIDSS pipeline
#
set -o errexit -o pipefail -o noclobber -o nounset
! getopt --test > /dev/null 
if [[ ${PIPESTATUS[0]} -ne 4 ]]; then
    echo '`getopt --test` failed in this environment.'
    exit 1
fi

OPTIONS=vo:t:n:s:
LONGOPTS=verbose,output_dir:tumour_bam:,normal_bam,sample,threads:
! PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    # e.g. return value is 1
    #  then getopt has complained about wrong arguments to stdout
    exit 2
fi
eval set -- "$PARSED"
run_dir=""
tumour_bam=""
normal_bam=""
verbose="n"
threads=$(nproc)
while true; do
    case "$1" in
        -v|--verbose)
            verbose=y
            shift
            ;;
		-n|--normal_bam)
            normal_bam=y
            shift 2
            ;;
        -o|--output_dir)
            run_dir="$2"
            shift 2
            ;;
		-t|--tumour_bam)
            tumour_bam="$2"
            shift 2
            ;;
		-s|--sample)
            sample="$2"
            shift 2
            ;;
        --)
            shift
            break
            ;;
		--threads)
			threads=$(printf -v int '%d\n' "$2" 2>/dev/null)
			shift 2
			;;
        *)
            echo "Programming error"
            exit 3
            ;;
    esac
done
if [[ ! -f "$tumour_bam" ]] ; then
	echo "Missing tumour BAM"
	exit 1
fi
if [[ ! -f "$normal_bam" ]] ; then
	echo "Missing normal BAM"
	exit 1
fi
mkdir -p "$run_dir"
if [[ ! -d "$run_dir" ]] ; then
	echo "Unable to create $run_dir"
	exit 1
fi
if [[ -z "$sample" ]] ; then
	sample_name=$(basename $tumour_bam .bam)
fi
if [[ $threads -lt 1 ]] ; then
	echo "Illegal thread count: $threads"
	exit 1
fi

joint_sample_name=$sample
ref_sample=${sample}_N
tumor_sample=${sample}_T



base_path=$(dirname $(readlink $0 || echo $0))
###
# Reference data
ref_genome=$base_path/refgenomes/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta
viral_ref_genome=$base_path/refgenomes/human_virus/human_virus.fa
encode_blacklist=$base_path/dbs/encode/ENCFF001TDO.bed
repeatmasker=$base_path/dbs/repeatmasker/hg19.fa.out
gc_profile=${base_path}/dbs/gc/GC_profile.1000bp.cnp
jar_dir=$base_path/jars
hmf_scripts=$base_path/hmfscripts/
gridss_jar=$(ls -1 $base_path/jars/gridss-2.2.3-gridss-jar-with-dependencies.jar)
purple_jar=$(ls -1 $base_path/jars/*purity-ploidy-estimator*.jar)
# TODO Download R libraries!
###
# External tools
export PATH=${base_path}/tools/circos_v0.69.6/bin/circos:$PATH

mkdir -p $run_dir/logs
log_prefix=$run_dir/logs/$(+%Y%m%d_%H%M%S).$HOSTNAME.$$

echo ############################################
echo # Running GRIDSS
echo ############################################
gridss_dir=$run_dir/gridss/
assembly_bam=$gridss_dir/$joint_sample_name.assembly.bam
gridss_raw_vcf=$gridss_dir/${joint_sample_name}.gridss.vcf.gz
gridss_somatic_full_vcf=$gridss_dir/${tumor_sample}.gridss.full.somatic.vcf
gridss_somatic_vcf=$gridss_dir/${tumor_sample}.gridss.somatic.vcf
if [[ ! -f $gridss_somatic_vcf ]] ; then
	mkdir -p $gridss_dir
	if ! which Rscript >/dev/null 2>&1 ; then
		echo "Missing R installation. Please add Rscript to PATH"
		exit 1
	fi
	if ! which java >/dev/null 2>&1 ; then
		echo "Missing java. Please add java 1.8 or later to PATH"
		exit 1
	fi
	if ! which bwa >/dev/null 2>&1 ; then
		echo "Missing bwa. Please add to PATH"
		exit 1
	fi
	gridss_jvm_args="
		-ea
		-Dreference_fasta=$ref_genome
		-Dsamjdk.create_index=true
		-Dsamjdk.use_async_io_read_samtools=true
		-Dsamjdk.use_async_io_write_samtools=true
		-Dsamjdk.use_async_io_write_tribble=true
		-Dsamjdk.buffer_size=2097152
		-Dgridss.gridss.output_to_temp_file=true
		-cp $gridss_jar "

	java -Xmx31g $gridss_jvm_args gridss.CallVariants \
		TMP_DIR=$gridss_dir \
		WORKING_DIR=$gridss_dir \
		REFERENCE_SEQUENCE="$ref_genome" \
		INPUT="$normal_bam" \
		INPUT="$tumor_sample" \
		OUTPUT="tmp.raw.$gridss_raw_vcf" \
		ASSEMBLY="$assembly_bam" \
		BLACKLIST="$encode_blacklist" \
		2>&1 | tee -a $log_prefix.gridss.CallVariants.log
		
	java -Xmx1G $gridss_jvm_args \
		gridss.AnnotateUntemplatedSequence \
		REFERENCE_SEQUENCE=$ref_genome \
		INPUT=tmp.raw.$gridss_raw_vcf \
		OUTPUT=tmp.human.$gridss_raw_vcf \
		WORKER_THREADS=$threads \
		2>&1 | tee -a $log_prefix.gridss.AnnotateUntemplatedSequence.human.log
		
	java -Xmx1G $gridss_jvm_args \
		gridss.AnnotateUntemplatedSequence \
		REFERENCE_SEQUENCE=$viral_ref_genome \
		INPUT=tmp.human.$gridss_raw_vcf \
		OUTPUT=$gridss_raw_vcf \
		WORKER_THREADS=$threads \
		2>&1 | tee -a $log_prefix.gridss.AnnotateUntemplatedSequence.viral.log
		
	rm tmp*$gridss_raw_vcf* 2>/dev/null
	unzipped_vcf=$(dirname ${original_vcf})/$(basename -s .gz ${original_vcf})
	# workaround for https://github.com/Bioconductor/VariantAnnotation/issues/19
	gunzip -c ${gridss_raw_vcf} | awk ' { if (length($0) >= 4000) { gsub(":0.00:", ":0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000:")} ; print $0  } ' > tmp.$gridss_raw_vcf.vcf
	Rscript $hmf_scripts/gridss/gridss_somatic_filter.R \
		-p ${gridss_pon} \
		-i tmp.$gridss_raw_vcf.vcf \
		-o ${gridss_somatic_vcf} \
		-f ${gridss_somatic_full_vcf} \
		-s $hmf_scripts/gridss/ \
		--gc
	rm tmp.$gridss_raw_vcf.vcf  2>/dev/null
	mv ${gridss_somatic_vcf}.bgz ${gridss_somatic_vcf}.gz
	mv ${gridss_somatic_vcf}.bgz.tbi ${gridss_somatic_vcf}.gz.tbi
	mv ${gridss_somatic_full_vcf}.bgz ${gridss_somatic_full_vcf}.gz
	mv ${gridss_somatic_full_vcf}.bgz.tbi ${gridss_somatic_full_vcf}.gz.tbi
else
	echo "Found $gridss_somatic_vcf, skipping GRIDSS" 
fi
if [[ ! -s $gridss_somatic_vcf.gz ]] ; then
	echo "Failed to generate GRIDSS VCF"
	exit 1
fi

echo ############################################
echo # Running PURPLE
echo ############################################
purple_output=${run_dir}/purple

somatic_vcf="????TODO?????"

java -Dorg.jooq.no-logo=true -Xmx16G -Xms4G \
    -jar ${purple_jar} \
    -somatic_vcf $somatic_vcf \
    -structural_vcf $gridss_somatic_vcf.gz \
    -circos circos \
    -run_dir ${run_dir} \
    -ref_genome hg19 \
    -output_dir ${purple_output} \
    -gc_profile ${gc_profile} \
    -sv_recovery_vcf $gridss_somatic_full_vcf.gz
purple_raw_vcf=$purple_output/????TODO????*.purple.sv.vcf.gz

echo ############################################
echo # Running repeatmasker annotation
echo ############################################
purple_annotated_vcf=${purple_raw_vcf/.purple.sv.vcf/.purple.ann.sv.vcf}

Rscript $hmf_scripts/gridss/gridss_annotate_insertions_repeatmaster.R \
	--input $purple_raw_vcf.gz \
	--output $purple_annotated_vcf \
	--repeatmasker $repeatmasker \
	--scriptdir $hmf_scripts/gridss/
mv ${purple_annotated_vcf}.bgz ${purple_annotated_vcf}.gz
mv ${purple_annotated_vcf}.bgz.tbi ${purple_annotated_vcf}.gz.tbi

echo ############################################
echo # Running FUCHSIA
echo ############################################







































