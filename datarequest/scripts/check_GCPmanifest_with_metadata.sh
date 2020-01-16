#!/usr/bin/env bash

#put metadata file path at $1.
#put manifest file path at $2.
# note: comm function
#comm [-1] [-2] [-3 ] file1 file2
# -1 Suppress the output column of lines unique to file1.
# -2 Suppress the output column of lines unique to file2.
# -3 Suppress the output column of lines duplicated in file1 and file2.

# make dir for temporary files
mkdir temp_GCP_check

### Read in patients and samples selected for DR
awk -F"\t" 'FNR > 1 {print $1}' $1 | sort | uniq > temp_GCP_check/patientId_metadata.tsv
awk -F"\t" 'FNR > 1 {print $2}' $1 | sort | uniq > temp_GCP_check/samplesId_metadata.tsv


echo ''
###############
echo 'LOOK AT OVERLAP BETWEEN PATIENTS SELECTED FOR THE DATA REQUEST AND DNA FILES SHARED IN THE GCP MANIFEST'
echo ''

wc -l temp_GCP_check/patientId_metadata.tsv

jq '.patients | select(.[] | .samples | .[] | .data | .[] | .tags | .moleculetype=="DNA") | .[] | .patientId'  $2 > temp_GCP_check/patientId_DNA_GCP.tsv
sed 's/\"//g' temp_GCP_check/patientId_DNA_GCP.tsv  > testfile.tmp && mv testfile.tmp temp_GCP_check/patientId_DNA_GCP.tsv
cat temp_GCP_check/patientId_DNA_GCP.tsv |  sort | uniq  > testfile.tmp && mv testfile.tmp temp_GCP_check/patientId_DNA_GCP.tsv
wc -l temp_GCP_check/patientId_DNA_GCP.tsv

echo 'Overlap between two files:'
comm -1 -2 <(sort temp_GCP_check/patientId_metadata.tsv) <(sort temp_GCP_check/patientId_DNA_GCP.tsv) | wc -l
echo 'extra samples patientId_metadata.tsv:'
comm -2 -3 <(sort temp_GCP_check/patientId_metadata.tsv) <(sort temp_GCP_check/patientId_DNA_GCP.tsv)
echo 'extra samples patientId_DNA_GCP.tsv:'
comm -1 -3 <(sort temp_GCP_check/patientId_metadata.tsv) <(sort temp_GCP_check/patientId_DNA_GCP.tsv)

#######
echo ''
#######

wc -l temp_GCP_check/samplesId_metadata.tsv

jq '.patients | select(.[] | .samples | .[] | .data | .[] | .tags | .moleculetype=="DNA") | .[] | .samples | .[] | .sampleId' $2 > temp_GCP_check/sampleId_DNA_GCP.tsv
sed 's/\"//g' temp_GCP_check/sampleId_DNA_GCP.tsv  > testfile.tmp && mv testfile.tmp temp_GCP_check/sampleId_DNA_GCP.tsv
cat temp_GCP_check/sampleId_DNA_GCP.tsv |  sort | uniq  > testfile.tmp && mv testfile.tmp temp_GCP_check/sampleId_DNA_GCP.tsv
wc -l temp_GCP_check/sampleId_DNA_GCP.tsv

echo 'Overlap between two files:'
comm -1 -2 <(sort temp_GCP_check/samplesId_metadata.tsv) <(sort temp_GCP_check/sampleId_DNA_GCP.tsv) | wc -l
echo 'extra samples patientId_metadata.tsv:'
comm -2 -3 <(sort temp_GCP_check/samplesId_metadata.tsv) <(sort temp_GCP_check/sampleId_DNA_GCP.tsv)
echo 'extra samples patientId_DNA_GCP.tsv:'
comm -1 -3 <(sort temp_GCP_check/samplesId_metadata.tsv) <(sort temp_GCP_check/sampleId_DNA_GCP.tsv)

#######
echo ''
#######


jq '.patients | .[] | .samples | .[] | .data | .[] | .tags | .[]' $2 > temp_GCP_check/files_GCP.tsv

echo 'number of seperate DNA BAM files shared within GCP:'
grep 'DNA' temp_GCP_check/files_GCP.tsv | wc -l



echo ''
echo ''
###############
echo 'LOOK AT OVERLAP BETWEEN PATIENTS SELECTED FOR THE DATA REQUEST AND RNA FILES SHARED IN THE GCP MANIFEST'
echo ''

wc -l temp_GCP_check/patientId_metadata.tsv

jq '.patients | select(.[] | .samples | .[] | .data | .[] | .tags | .moleculetype=="RNA") | .[] | .patientId'  $2 > temp_GCP_check/patientId_RNA_GCP.tsv
sed 's/\"//g' temp_GCP_check/patientId_RNA_GCP.tsv  > testfile.tmp && mv testfile.tmp temp_GCP_check/patientId_RNA_GCP.tsv
cat temp_GCP_check/patientId_RNA_GCP.tsv |  sort | uniq  > testfile.tmp && mv testfile.tmp temp_GCP_check/patientId_RNA_GCP.tsv
wc -l temp_GCP_check/patientId_RNA_GCP.tsv

echo 'Overlap between two files:'
comm -1 -2 <(sort temp_GCP_check/patientId_metadata.tsv) <(sort temp_GCP_check/patientId_RNA_GCP.tsv) | wc -l
echo 'extra samples patientId_metadata.tsv:'
#comm -2 -3 <(sort temp_GCP_check/patientId_metadata.tsv) <(sort temp_GCP_check/patientId_RNA_GCP.tsv)
comm -2 -3 <(sort temp_GCP_check/patientId_metadata.tsv) <(sort temp_GCP_check/patientId_RNA_GCP.tsv) | wc -l
echo 'extra samples patientId_DNA_GCP.tsv:'
comm -1 -3 <(sort temp_GCP_check/patientId_metadata.tsv) <(sort temp_GCP_check/patientId_RNA_GCP.tsv)


#######
echo ''
#######

wc -l temp_GCP_check/samplesId_metadata.tsv

jq '.patients | select(.[] | .samples | .[] | .data | .[] | .tags | .moleculetype=="RNA") | .[] | .samples | .[] | .sampleId' $2 > temp_GCP_check/sampleId_RNA_GCP.tsv
sed 's/\"//g' temp_GCP_check/sampleId_RNA_GCP.tsv  > testfile.tmp && mv testfile.tmp temp_GCP_check/sampleId_RNA_GCP.tsv
cat temp_GCP_check/sampleId_RNA_GCP.tsv |  sort | uniq  > testfile.tmp && mv testfile.tmp temp_GCP_check/sampleId_RNA_GCP.tsv
wc -l temp_GCP_check/sampleId_RNA_GCP.tsv

echo 'Overlap between two files:'
comm -1 -2 <(sort temp_GCP_check/samplesId_metadata.tsv) <(sort temp_GCP_check/sampleId_RNA_GCP.tsv) | wc -l
echo 'extra samples patientId_metadata.tsv:'
#comm -2 -3 <(sort temp_GCP_check/samplesId_metadata.tsv) <(sort temp_GCP_check/sampleId_RNA_GCP.tsv)
comm -2 -3 <(sort temp_GCP_check/samplesId_metadata.tsv) <(sort temp_GCP_check/sampleId_RNA_GCP.tsv) | wc -l
echo 'extra samples patientId_DNA_GCP.tsv:'
comm -1 -3 <(sort temp_GCP_check/samplesId_metadata.tsv) <(sort temp_GCP_check/sampleId_RNA_GCP.tsv)

#######
echo ''
#######


jq '.patients | .[] | .samples | .[] | .data | .[] | .tags | .[]' $2 > temp_GCP_check/files_GCP.tsv

echo 'number of seperate RNA FASTQ files shared within GCP:'
grep 'RNA' temp_GCP_check/files_GCP.tsv | wc -l


# remove all temporary files
rm -r temp_GCP_check