#!/bin/bash

# Definition of tools used
sullivan=</path/to/sullivan.sh>
picard=</path/to/picard.sh>

# Params
hmfreg=<HMFregXXXX>
externalid=<CPCTxxxxxxxxR/T>
internalid=<FRxxxxxxxx>

sequencerrun1=<DATE_RUNNAME1_FLOWCELL>
sequencerrun2=<DATE_RUNNAME2_FLOWCELL>
mergerun=<true/false>
appendreg=<true/false>

if [ ${mergerun} == true ]; then
	pipelinerun=${sequencerrun1}_merge_${sequencerrun2}
else
	pipelinerun=${sequencerrun1}
fi

if [ ${appendreg} == true ]; then
	pipelinerun=${pipelinerun}_${hmfreg}
fi

# Paths
bamfile=/data/processed_archive/${pipelinerun}/${externalid}/mapping/${externalid}_dedup.realigned.bam
sequencerundir=/data/sequencer_archive/${sequencerrun1}/Data/Intensities/BaseCalls/${hmfreg}/${internalid}
if [ $mergerun == true ]; then
	sequencerundirmerge=/data/sequencer_archive/${sequencerrun2}/Data/Intensities/BaseCalls/${hmfreg}/${internalid}
fi
tmpfastqdir=/data/tmp/fastq

# Parameters
numRecords=1000000

# Run tools
rm ${tmpfastqdir}/*
$picard SamToFastq I=${bamfile} ODIR=${tmpfastqdir} OPRG=true RGT=ID NON_PF=true RC=true
if [ $mergerun  == true ]; then
	$sullivan -origpath ${sequencerundir} -recpath ${tmpfastqdir} -mergeorigpath ${sequencerundirmerge} -n ${numRecords} -d
else
	$sullivan -origpath ${sequencerundir} -recpath ${tmpfastqdir} -n ${numRecords} -d -a
fi

# For easy check whether everything looks sane:
ls -lah ${tmpfastqdir}
# Add the commands to clean up to the screen, should be executed in case QC succeeded based on manual inspection
echo "rm ${sequencerundir}/*.fastq.gz"
if [ $mergerun == true ]; then
	echo "rm ${sequencerundirmerge}/*.fastq.gz"
fi
