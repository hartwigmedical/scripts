#!/bin/bash

## This script has to be executed from the
## bcl2fastq converter output directory
## (there were the Undetermined fastqs are)

echo "# Current dir"
echo $PWD

echo ""
echo "# size of Undetermined fqs"
du -hc Undetermined*fastq.gz | tail -1

echo ""
echo "# size of HMF project fqs 2 deep"
du -hc *HMF*/*/*fastq.gz | tail -1

echo ""
echo "# overview of all fqs 2 deep"
for sample in *HMF*/*/; do echo "----- $sample -----"; du -shc $sample/*fastq.gz; done

echo ""
echo "# size of HMF project fqs 1 deep"
du -hc *HMF*/*fastq.gz | tail -1

echo ""
echo "# overview of all fqs 1 deep"
# example: HMFproject/CPCT02010234R_S5_L001_R1_001.fastq.gz
#for sampleFQ in *HMF*/*fastq.gz; do 
du -shc *HMF*/*fastq.gz

echo ""
echo 'Done'
