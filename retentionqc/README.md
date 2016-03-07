This script checks whether all data from original fastq files have been retained in the final BAM, by taking the folowing steps:
 - Run picard SamToFastq to convert the final BAM into fastq files
 - Run sullivan to compare the recreated fastq files with the original fastq files.
