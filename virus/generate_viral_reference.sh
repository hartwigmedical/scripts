#!/bin/bash
HOMOSAPIENS_TAX_ID=9606

wget ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.genomic.fna.gz
wget ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.tsv
gunzip virushostdb.genomic.fna.gz


HUMAN_VIRUS_REFSEQ=$(cut -f 4,8 virushostdb.tsv  | grep $HOMOSAPIENS_TAX_ID | cut -f 1 | tr '\n' ',' | tr -d '\t ')
# filterbyname.sh is part of BBMap
filterbyname.sh  in=virushostdb.genomic.fna out=human_virus.fa names=$HUMAN_VIRUS_REFSEQ include=t
RepeatMasker -no_is -pa $(nproc) -s -noint -norna -species human -html human_virus.fa

rm human_virus.fa
ln -s human_virus.fa.masked human_virus.fa

samtools faidx human_virus.fa
bwa index hg19_virus.fa

for REFSEQ in $(echo $HUMAN_VIRUS_REFSEQ | tr , " ") ; do 
	echo "$REFSEQ,$(grep $REFSEQ virushostdb.tsv | cut -f 2,3,11,12)" >> virus_refseq_to_human_readable.tsv
done
