#!/bin/bash

BAM_SLICER=/data/common/tools/bam-slicer-v1.0/bam-slicer.jar

SAMPLE=$1
BED_FILE=$2

if [ -z "$SAMPLE" ] || [ -z "$BED_FILE" ];
  then
    echo "Usage: $ ./bam-slicer.sh sample bed_file"
    echo "   sample	sample to search for. e.g CPCT11111111T"
    echo "   bed_file	bed file containing regions to be sliced"
    exit 1
fi

# get validated somatic.ini/cpct.ini runs for this sample, sort in descending order of pipeline version then take first row
RUN_ROW=`query_sbp_api -type runs | grep -P "${SAMPLE}\t" | grep "Somatic.ini\|CPCT.ini" | grep Validated | sort -nrk7,7 | head -1`
RUN_NAME=`echo $RUN_ROW | cut -d ' ' -f1`
RUN_BUCKET=`echo $RUN_ROW | cut -d ' ' -f5`
BAM_KEY="$RUN_NAME/$SAMPLE/mapping/${SAMPLE}_dedup.realigned.bam"
BAM_INDEX_KEY="${BAM_KEY}.bai"
SLICED_BAM="${SAMPLE}_dedup.realigned.sliced.bam"

java -Dsamjdk.buffer_size=0 \
    -Xmx4G \
    -jar $BAM_SLICER \
    -s3 \
    -bucket $RUN_BUCKET \
    -input $BAM_KEY \
    -index $BAM_INDEX_KEY \
    -bed $BED_FILE \
    -output $SLICED_BAM \
    -max_chunks 2000
