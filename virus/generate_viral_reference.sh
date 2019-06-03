#!/bin/bash
HOMOSAPIENS_TAX_ID=9606

if [[ ! -f virushostdb.genomic.fna ]] ; then
	wget ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.genomic.fna.gz
	gunzip virushostdb.genomic.fna.gz
fi
if [[ ! -f virushostdb.tsv ]] ; then
	wget ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.tsv
fi

HUMAN_VIRUS_REFSEQ=$(cut -f 4,8 virushostdb.tsv  | grep $HOMOSAPIENS_TAX_ID | cut -f 1 | tr '\n' ',' | tr -d '\t ')
# # filterbyname.sh is part of BBMap
filterbyname.sh  in=virushostdb.genomic.fna out=human_virus.fa names=$HUMAN_VIRUS_REFSEQ include=t
RepeatMasker -no_is -pa $(nproc) -s -noint -norna -species human -html human_virus.fa

rm human_virus.fa
ln -s human_virus.fa.masked human_virus.fa

samtools faidx human_virus.fa
bwa index human_virus.fa

Rscript extract.R