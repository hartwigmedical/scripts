#!/usr/bin/env bash

HOMOSAPIENS_TAX_ID=9606
HUMAN_REF=/data/common/refgenomes/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fa
REF=hg19_virus.fa

wget ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.genomic.fna.gz
wget ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.tsv
gunzip virushostdb.genomic.fna.gz


HUMAN_VIRUS_REFSEQ=$(cut -f 4,8 virushostdb.tsv  | grep $HOMOSAPIENS_TAX_ID | cut -f 1 | tr '\n' ',' | tr -d '\t ')
# filterbyname.sh is part of BBMap
filterbyname.sh  in=virushostdb.genomic.fna out=human_virus.fa names=$HUMAN_VIRUS_REFSEQ include=t
RepeatMasker -no_is -pa $(nproc) -s -noint -norna -species human -html human_virus.fa
cp $HUMAN_REF $REF
cat human_virus.fa.masked >> $REF

samtools faidx human_virus.fa
# Generate .alt file used by bwa to determine ALT contig mappings
samtools view -H ${HUMAN_REF/.fa/.dict} | grep "@SQ" > $REF.alt
for REFSEQ in $(echo $HUMAN_VIRUS_REFSEQ | tr , " ") ; do 
	echo "$REFSEQ	4	*	0	0	$(grep $REFSEQ human_virus.fa.fai | cut -f 2)M	*	0	0	*	*" >> $REF.alt
done
samtools faidx hg19_virus.fa
bwa index hg19_virus.fa

for REFSEQ in $(echo $HUMAN_VIRUS_REFSEQ | tr , " ") ; do 
	echo "$REFSEQ,$(grep $REFSEQ virushostdb.tsv | cut -f 2,3,11,12)" >> hg19_virus_refseq_to_human_readable.tsv
done
