#!/bin/bash
for sheet in /data1/illumina_data/*/SampleSheet.csv; do 
	currHost=`hostname`
	runBase=`dirname $sheet`
	readmeFile=$runBase"/README"
	echo ""
	#echo "# ----- HMFreg";
	echo "# "$currHost;
	echo "# [SSHEET] "$runBase;
	
	## check images presence
	if [ -d "$runBase/Thumbnail_Images/L001" ]; then
		echo "# [NOTE] IMAGES still present";
	fi

	## print readme file path if present
	if [ -e $readmeFile ]; then
		echo "# [README] $readmeFile";
	fi

	## print internal HMF run-name
	cat $sheet | grep "ExperimentName" | cut -d',' --output-delimiter ": " -f 1,2
	## print sample-id, sample_name, submission-id, description
	cat $sheet | sed -e '1,/Lane,Sample_ID/d' | cut -d',' --output-delimiter " " -f2,3,8,9 | sort | uniq; 
done
