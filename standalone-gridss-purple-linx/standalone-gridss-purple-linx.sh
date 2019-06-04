#!/bin/bash
#
# Stand-alone GRIDSS-PURPLE-Linx pipeline
#
set -o errexit -o pipefail -o noclobber -o nounset
! getopt --test > /dev/null 
if [[ ${PIPESTATUS[0]} -ne 4 ]]; then
    echo '`getopt --test` failed in this environment.'
    exit 1
fi

OPTIONS=v:o:t:n:s:
LONGOPTS=vcf,output_dir:tumour_bam:,normal_bam,sample,threads:
! PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    # e.g. return value is 1
    #  then getopt has complained about wrong arguments to stdout
    exit 2
fi
eval set -- "$PARSED"
run_dir=""
ref_dir=""
tumour_bam=""
normal_bam=""
snvindel_vcf=""
threads=$(nproc)
cleanup="y"
while true; do
    case "$1" in
        -v|--vcf)
            snvindel_vcf="$2"
            shift 2
            ;;
		-n|--normal_bam)
            normal_bam="$2"
            shift 2
            ;;
        -o|--output_dir)
            run_dir="$2"
            shift 2
            ;;
		-r|--ref_dir)
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
        --snvindel_vcf)
		snvindel_vcf="$2"
            shift
            ;;
		--threads)
			threads=$(printf -v int '%d\n' "$2" 2>/dev/null)
			shift 2
			;;
		--)
            shift
            break
            ;;
        *)
            echo "Programming error"
            exit 3
            ;;
    esac
done
if [[ ! -f "$snvindel_vcf" ]] ; then
	echo "Missing SNV VCF. A SNV VCF with the AD genotype field populated is required."
	echo "Use the script for generating this VCF with strelka if you have not already generated a compatible VCF."
	exit 1
fi
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
if [[ ! -d "$ref_dir" ]] ; then
	echo "Could not find Hartwig reference directory $ref_dir"
	exit 1
fi
if [[ -z "$sample" ]] ; then
	sample_name=$(basename $tumour_bam .bam)
fi
if [[ "$threads" -lt 1 ]] ; then
	echo "Illegal thread count: $threads"
	exit 1
fi

joint_sample_name=$sample
ref_sample=${sample}_N
tumor_sample=${sample}_T



base_path=$(dirname $(readlink $0 || echo $0))
###
# Reference data
ref_genome=$ref_dir/refgenomes/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta
viral_ref_genome=$ref_dir/refgenomes/human_virus/human_virus.fa
encode_blacklist=$ref_dir/dbs/encode/ENCFF001TDO.bed
repeatmasker=$ref_dir/dbs/repeatmasker/hg19.fa.out
baf_snps=$ref_dir/dbs/amber/GermlineHetPon.hg19.bed
gc_profile=$ref_dir/dbs/gc/GC_profile.1000bp.cnp
gridss_pon=$ref_dir/dbs/gridss/pon/
strelka_config=$ref_dir/settings/strelka/strelka_config_bwa_genome.ini
hmf_scripts=$ref_dir/hmfscripts/
gridss_jar=$(ls -1 $ref_dir/jars/*gridss*.jar)
purple_jar=$(ls -1 $ref_dir/jars/*purity-ploidy-estimator*.jar)
amber_jar=$(ls -1 $ref_dir/jars/*amber*.jar)
cobalt_jar=$(ls -1 $ref_dir/jars/*cobalt*.jar)

for program in bwa sambamba samtools circos Rscript java ; do
	if ! which $program > /dev/null ; then
		echo "Missing required dependency $program. $program must be on PATH"
		exit 1
	fi
done

for rpackage in tidyverse devtools assertthat testthat NMF stringdist stringr argparser R.cache "copynumber" StructuralVariantAnnotation "VariantAnnotation" "rtracklayer" "BSgenome" "org.Hs.eg.db" "TxDb.Hsapiens.UCSC.hg19.knownGene" "BSgenome.Hsapiens.UCSC.hg19" ; do
	if ! Rscript -e "installed.packages()" | grep $rpackage > /dev/null ; then
		echo "Missing R package $rpackage"
		echo "All required R packages can be installed by running Rscript install_rpackages.R with appropriate permissions."
		exit 1
	fi
done

mkdir -p $run_dir/logs
log_prefix=$run_dir/logs/$(+%Y%m%d_%H%M%S).$HOSTNAME.$$

JVM_MEMORY_USAGE_ARGS="
	-XX:+UnlockExperimentalVMOptions
	-XX:+UseCGroupMemoryLimitForHeap
	-XX:MaxRAMFraction=1
	-XshowSettings:vm"

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
	gridss_jvm_args="
		-ea
		$JVM_MEMORY_USAGE_ARGS
		-Dreference_fasta=$ref_genome
		-Dsamjdk.create_index=true
		-Dsamjdk.use_async_io_read_samtools=true
		-Dsamjdk.use_async_io_write_samtools=true
		-Dsamjdk.use_async_io_write_tribble=true
		-Dsamjdk.buffer_size=$((4 * 1024 * 1024))
		-Dgridss.gridss.output_to_temp_file=true
		-cp $gridss_jar "
		
	java $gridss_jvm_args gridss.CallVariants \
		TMP_DIR=$gridss_dir \
		WORKING_DIR=$gridss_dir \
		REFERENCE_SEQUENCE="$ref_genome" \
		INPUT="$normal_bam" \
		INPUT="$tumor_sample" \
		OUTPUT="tmp.raw.$gridss_raw_vcf" \
		ASSEMBLY="$assembly_bam" \
		BLACKLIST="$encode_blacklist" \
		2>&1 | tee -a $log_prefix.gridss.CallVariants.log
		
	java $gridss_jvm_args \
		gridss.AnnotateUntemplatedSequence \
		REFERENCE_SEQUENCE=$ref_genome \
		INPUT=tmp.raw.$gridss_raw_vcf \
		OUTPUT=tmp.human.$gridss_raw_vcf \
		WORKER_THREADS=$threads \
		2>&1 | tee -a $log_prefix.gridss.AnnotateUntemplatedSequence.human.log
		
	java $gridss_jvm_args \
		gridss.AnnotateUntemplatedSequence \
		REFERENCE_SEQUENCE=$viral_ref_genome \
		INPUT=tmp.human.$gridss_raw_vcf \
		OUTPUT=$gridss_raw_vcf \
		WORKER_THREADS=$threads \
		2>&1 | tee -a $log_prefix.gridss.AnnotateUntemplatedSequence.viral.log
	# workaround for https://github.com/Bioconductor/VariantAnnotation/issues/19
	gunzip -c ${gridss_raw_vcf} | awk ' { if (length($0) >= 4000) { gsub(":0.00:", ":0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000:")} ; print $0  } ' > tmp.decompressed.$gridss_raw_vcf.vcf
	Rscript $hmf_scripts/gridss/gridss_somatic_filter.R \
		-p ${gridss_pon} \
		-i tmp.decompressed.$gridss_raw_vcf.vcf \
		-o ${gridss_somatic_vcf} \
		-f ${gridss_somatic_full_vcf} \
		-s $hmf_scripts/gridss/ \
		--gc \
		2>&1 | tee -a $log_prefix.gridss.somatic_filter.log
	if [[ $cleanup == "y" ]] ; then
		rm tmp.$gridss_raw_vcf.vcf* tmp.human.$gridss_raw_vcf* tmp.decompressed.$gridss_raw_vcf.vcf* 2>/dev/null
	fi
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
echo # Running Amber
echo ############################################
mkdir -p $run_dir/amber
amber_pileup_N=$run_dir/amber/$ref_sample.amber.pileup
amber_pileup_T=$run_dir/amber/$ref_sample.amber.pileup

sambamba mpileup \
	-t $threads \
	--tmpdir=$run_dir/amber \
	-L $baf_snps \
	$normal_bam \
	--samtools "-q 1 -f $ref_genome" > $amber_pileup_N
	2>&1 | tee $log_prefix.amber.N.log

sambamba mpileup \
	-t $threads \
	--tmpdir=$run_dir/amber \
	-L $baf_snps \
	$tumor_bam \
	--samtools "-q 1 -f $ref_genome" > $amber_pileup_T
	2>&1 | tee $log_prefix.amber.T.log

java $JVM_MEMORY_USAGE_ARGS \
	-jar $amber_jar \
	-sample $tumor_sample \
	-reference $amber_pileup_N \
	-tumor $amber_pileup_T \
	-output_dir $run_dir/amber \
	2>&1 | tee $log_prefix.amber.log

if [[ $cleanup == "y" ]] ; then
	rm $amber_pileup_N $amber_pileup_T $run_dir/amber/*.pcf1
fi

echo ############################################
echo # Running Cobalt
echo ############################################
mkdir -p $run_dir/cobalt
java $JVM_MEMORY_USAGE_ARGS -jar $cobalt_jar \
    -threads $threads \
    -reference $ref_sample \
    -reference_bam $normal_bam \
    -tumor $tumor_sample \
    -tumor_bam $tumor_bam \
    -output_dir $run_dir/cobalt \
    -gc_profile $gc_profile \
	2>&1 | tee $log_prefix.cobalt.log

if [[ $cleanup == "y" ]] ; then
	rm -f "[% dirs.cobalt %]"/*.pcf1
	rm -f "[% dirs.cobalt %]"/*.ratio
	rm -f "[% dirs.cobalt %]"/*.ratio.gc
	rm -f "[% dirs.cobalt %]"/*.count
fi

echo ############################################
echo # Running PURPLE
echo ############################################
purple_output=${run_dir}/purple

java -Dorg.jooq.no-logo=true $JVM_MEMORY_USAGE_ARGS \
    -jar ${purple_jar} \
    -somatic_vcf $somatic_vcf \
    -structural_vcf $gridss_somatic_vcf.gz \
    -circos circos \
    -run_dir ${run_dir} \
    -ref_genome hg19 \
	-cobalt $run_dir/cobalt \
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
echo # Running Linx
echo ############################################







































