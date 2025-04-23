# Requirements
## wget awk grep
## bedtools 2.31.1 (mask fasta)
## bwa-mem2 2.2.1 (index)
## bwa-mem 0.7.17 (index)
## gatk 4.6.0.0 (index)
## picard tools 2.18.27 (sequence dictionary)

## runtime +- 4 hours on n2d-highmem-16 instance

# Download v0 of GRCh38
## https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0
wget https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta

# Download masking definitions
## https://www.nature.com/articles/s41587-021-01158-1
## https://doi.org/10.1016/j.jmoldx.2021.10.013
# -- wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_GRC_exclusions_T2Tv2.bed
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_GRC_exclusions.bed
cat GCA_000001405.15_GRCh38_GRC_exclusions.bed > exclusions.bed
### #sequence	sequenceStart	sequenceEnd	exclusion
### chr10_KI270825v1_alt	0	188315	contamination
grep -e "^>" Homo_sapiens_assembly38.fasta > contigs.txt

# PhiX
## Download sequence
wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/PhiX/Illumina/RTA/PhiX_Illumina_RTA.tar.gz
## Unzip
tar -xzvf PhiX_Illumina_RTA.tar.gz
## Add to fasta reference
echo -e ">phix\tAC:NC_001422\tgi:0000\tLN:5386\trl:unplaced" > phix_genome.fa
tail -n +2 PhiX/Illumina/RTA/Sequence/Chromosomes/phix.fa >> phix_genome.fa
cat phix_genome.fa >> Homo_sapiens_assembly38.fasta

# Add ALTs to the exclusions
## >chr1_KI270762v1_alt  AC:KI270762.1  gi:568335926  LN:354444  rg:chr1:2448811-2791270  rl:alt-scaffold  M5:b0397179e5a92bb7a3300b68e45a9f72  AS:GRCh38
grep "_alt" contigs.txt | awk -v OFS='\t' '{sub(/LN:/, "", $4) sub("^>", "", $1)} {print $1, 0, $4, "alt"}' >> exclusions.bed

# Add HLAs to the exclusions
## >HLA-A*01:01:01:01	HLA00001 3503 bp
grep "HLA" contigs.txt | awk -v OFS='\t' '{sub("^>", "", $1)} {print $1, 0, $3, "HLA"}' >> exclusions.bed

# Apply hard masks
## https://bioweb.pasteur.fr/docs/modules/bedtools/2.32.1/content/maskfastafromBed.html
bedtools2/bin/maskFastaFromBed \
-fi Homo_sapiens_assembly38.fasta \
-bed exclusions.bed \
-fo GRCh38_masked_exclusions_alts_hlas.fasta

# md5 GRCh38_masked_exclusions_alts_hlas.fasta > masking.md5sum
## MD5 (GRCh38_masked_exclusions_alts_hlas.fasta) = 66fe40b3a5666dedb670f34396decdfc
md5sum GRCh38_masked_exclusions_alts_hlas.fasta > masking.md5sum
echo "66fe40b3a5666dedb670f34396decdfc  GRCh38_masked_exclusions_alts_hlas.fasta" > tomatch.md5sum
diff masking.md5sum tomatch.md5sum

if [[ $(diff masking.md5sum tomatch.md5sum) ]]; then
    echo "There is an issue in the md5 checksums please check file integrity"
else
    echo "all good to go"
fi

bwa-mem2 index GRCh38_masked_exclusions_alts_hlas.fasta
gatk BwaMemIndexImageCreator -I ./GRCh38_masked_exclusions_alts_hlas.fasta -O GRCh38_masked_exclusions_alts_hlas.fasta.img
samtools faidx GRCh38_masked_exclusions_alts_hlas.fasta
samtools index GRCh38_masked_exclusions_alts_hlas.fasta
bwa index GRCh38_masked_exclusions_alts_hlas.fasta
java -jar picard.jar CreateSequenceDictionary R=./GRCh38_masked_exclusions_alts_hlas.fasta O=GRCh38_masked_exclusions_alts_hlas.fasta.dict
cp GRCh38_masked_exclusions_alts_hlas.fasta.dict GRCh38_masked_exclusions_alts_hlas.dict
