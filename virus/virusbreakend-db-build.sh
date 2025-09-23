#!/bin/bash

# gcloud compute ssh --zone "australia-southeast1-a" "dev-aus-virusbreakend-db-build" --project "hmf-aus" --tunnel-through-iap
# gcloud compute ssh --zone "europe-west4-a" "dev-syd-virusbreakend-db-build" --project "hmf-crunch" --tunnel-through-iap

## Set up standard directories ================================
sudo mkdir -p /opt/tools
sudo chmod -R 777 /opt/tools

sudo mkdir -p /data/output/
sudo chmod -R 777 /data/output


## Install standard dependencies ================================
sudo apt update && sudo apt upgrade
sudo apt install -y git
sudo apt install -y build-essential
sudo apt install -y libz-dev
sudo apt install -y rsync
sudo apt install -y curl
sudo apt install -y unzip
sudo apt install -y tree


## Install other dependencies ================================
## kraken2 --------------------------------
cd /opt/tools/
git clone https://github.com/DerrickWood/kraken2/
bash /opt/tools/kraken2/install_kraken2.sh /opt/tools/kraken2/

## edirect --------------------------------
sudo sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)" ## installs to /root/edirect/
sudo mv /root/edirect/ /opt/tools

## curl with openssl 1.1.1 --------------------------------
## edirect doesn't work with newer openssl versions
sudo apt update
sudo apt install -y build-essential autoconf libtool pkg-config \
    libz-dev libnghttp2-dev libidn2-dev libpsl-dev libssh2-1-dev libbrotli-dev \
    automake git

cd /tmp
wget https://www.openssl.org/source/openssl-1.1.1u.tar.gz
tar xzf openssl-1.1.1u.tar.gz
cd openssl-1.1.1u
./config --prefix=/opt/openssl-1.1.1 --openssldir=/opt/openssl-1.1.1
make -j$(nproc)
sudo make install

cd /tmp
git clone https://github.com/curl/curl.git
cd curl
autoreconf -fi
./configure --with-ssl=/opt/openssl-1.1.1 --prefix=/opt/curl-openssl1.1.1
make -j$(nproc)
sudo make install

## ncbi-datasets-cli --------------------------------
mkdir -p /opt/tools/ncbi-datasets-cli
cd /opt/tools/ncbi-datasets-cli
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat'
chmod +x datasets dataformat

## FASTA indexing --------------------------------
sudo apt install samtools

sudo apt install default-jre
mkdir -p /opt/tools/picard
wget -O /opt/tools/picard/picard-3.4.0.jar https://github.com/broadinstitute/picard/releases/download/3.4.0/picard.jar

## Build database ================================
export PATH=/opt/tools/kraken2/:/opt/tools/edirect/:/opt/curl-openssl1.1.1/bin:/opt/tools/ncbi-datasets-cli:/usr/local/bin:/usr/bin:/bin
export LD_LIBRARY_PATH=/opt/openssl-1.1.1/lib

DEST_DB_PATH=/data/output/virusbreakend-db_2025-09-09/
mkdir -p $DEST_DB_PATH
cd $DEST_DB_PATH

## Download taxonomy and sequences from kraken --------------------------------
mkdir -p $DEST_DB_PATH/kraken2-db
kraken2-build --download-taxonomy      --db $DEST_DB_PATH/kraken2-db --threads 16
kraken2-build --download-library viral --db $DEST_DB_PATH/kraken2-db --threads 16
kraken2-build --download-library human --db $DEST_DB_PATH/kraken2-db --threads 16
## Human sequences need to be added. Otherwise the virusbreakend output VCF will be empty (no integrations found)

## Output of the above kraken2-build commands
#
# kraken2-db/taxonomy/
# ├── accmap.dlflag
# ├── citations.dmp
# ├── delnodes.dmp
# ├── division.dmp
# ├── gc.prt
# ├── gencode.dmp
# ├── images.dmp
# ├── merged.dmp
# ├── names.dmp
# ├── nodes.dmp
# ├── nucl_gb.accession2taxid
# ├── nucl_wgs.accession2taxid
# ├── readme.txt
# ├── taxdump.dlflag
# ├── taxdump.tar.gz
# └── taxdump.untarflag
# 
# kraken2-db/library/viral/
# ├── assembly_summary.txt
# ├── library.fna
# ├── library.fna.masked
# ├── manifest.txt
# └── prelim_map.txt

## Download more sequences NCBI --------------------------------
mkdir -p $DEST_DB_PATH/supplementary

## !WARNING!
## This query will download ~2GB of data and can take several hours
QUERY='Viruses[Organism]'
QUERY+=' NOT cellular organisms[ORGN]'
QUERY+=' NOT wgs[PROP]'
QUERY+=' NOT AC_000001:AC_999999[pacc]' 
QUERY+=' NOT gbdiv syn[prop]'
QUERY+='AND (srcdb_refseq[PROP] OR nuccore genome samespecies[Filter])'
QUERY+=' AND "complete genome"[Title] NOT "UNVERIFIED"[Title]'

esearch -db nuccore -query "$QUERY" | efetch -format fasta > $DEST_DB_PATH/supplementary/ncbi_human_virus_sequences.fa

## !WARNING!
## The RefSeq viral database is continuously updated and can be unstable.
## As of 2024-06-18 (date of building the previous virusbreakend DB), HPV33 was missing from RefSeq
## As of 2025-09-11 (date of building the latest virusbreakend DB), HPV33 was restored in RefSeq

## Download table of viral neighbour genomes and their associated RefSeq genome --------------------------------

## Old command
# # Issue with neighbours: serotypes (e.g. HPV-45) are listed as neighbours of the species taxid
# # Use Neighbours for host taxid filtering. TODO: get a better source
# wget --output-document=$DEST_DB_PATH/supplementary/taxid10239.nbr "https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=10239&cmd=download2"

## The above wget command no longer works. We therefore use the file from the original virusbreakend database: 
## https://github.com/PapenfussLab/gridss/blob/master/VIRUSBreakend_Readme.md#reference-data
# gsutil cp gs://hmf-crunch-resources/virus/virusbreakend_db/virusbreakenddb_20210401/taxid10239.nbr $DEST_DB_PATH/supplementary/
gsutil cp gs://hmf-crunch-resources/virus/virusbreakend_db/2024-06-18_virus-human-minimal/taxid10239.nbr $DEST_DB_PATH/supplementary/

## Build virusbreakend database --------------------------------

## Link kraken2 database files to preserve the original (for debugging purposes)
mkdir -p $DEST_DB_PATH/virusbreakend-db/
mkdir -p $DEST_DB_PATH/virusbreakend-db/taxonomy
mkdir -p $DEST_DB_PATH/virusbreakend-db/library/viral
mkdir -p $DEST_DB_PATH/virusbreakend-db/library/human

ln -sf $DEST_DB_PATH/kraken2-db/taxonomy/*      $DEST_DB_PATH/virusbreakend-db/taxonomy/
ln -sf $DEST_DB_PATH/kraken2-db/library/viral/* $DEST_DB_PATH/virusbreakend-db/library/viral/
ln -sf $DEST_DB_PATH/kraken2-db/library/human/* $DEST_DB_PATH/virusbreakend-db/library/human/

# Add supplementary files
ln -sf $DEST_DB_PATH/supplementary/taxid10239.nbr $DEST_DB_PATH/virusbreakend-db

## Add additional sequences
kraken2-build --add-to-library $DEST_DB_PATH/supplementary/ncbi_human_virus_sequences.fa --db $DEST_DB_PATH/virusbreakend-db

## Build
kraken2-build --build --db $DEST_DB_PATH/virusbreakend-db --threads 16

## Make indexes
for f in $(find $DEST_DB_PATH/ -name '*.fna') ; do
	rm -f ${f}.dict
	java -jar /opt/tools/picard/picard-3.4.0.jar CreateSequenceDictionary R=${f} O=${f}.dict

	samtools faidx $f
done

## Structure of the virusbreakend DB after the above commands
# virusbreakend-db
# ├── hash.k2d
# ├── library
# │   ├── added
# │   │   ├── hFZ5FM3Zgd.fna
# │   │   ├── hFZ5FM3Zgd.fna.dict
# │   │   ├── hFZ5FM3Zgd.fna.fai
# │   │   ├── hFZ5FM3Zgd.fna.masked
# │   │   ├── prelim_map.txt
# │   │   └── prelim_map_VvONseDUPx.txt
# │   ├── human
# │   │   ├── assembly_summary.txt -> /data/output/virusbreakend-db_2025-09-09//kraken2-db/library/human/assembly_summary.txt
# │   │   ├── library.fna -> /data/output/virusbreakend-db_2025-09-09//kraken2-db/library/human/library.fna
# │   │   ├── library.fna.dict
# │   │   ├── library.fna.fai
# │   │   ├── library.fna.masked -> /data/output/virusbreakend-db_2025-09-09//kraken2-db/library/human/library.fna.masked
# │   │   ├── manifest.txt -> /data/output/virusbreakend-db_2025-09-09//kraken2-db/library/human/manifest.txt
# │   │   ├── prelim_map.txt -> /data/output/virusbreakend-db_2025-09-09//kraken2-db/library/human/prelim_map.txt
# │   │   └── x -> /data/output/virusbreakend-db_2025-09-09//kraken2-db/library/human/x
# │   └── viral
# │       ├── assembly_summary.txt -> /data/output/virusbreakend-db_2025-09-09//kraken2-db/library/viral/assembly_summary.txt
# │       ├── library.fna -> /data/output/virusbreakend-db_2025-09-09//kraken2-db/library/viral/library.fna
# │       ├── library.fna.dict
# │       ├── library.fna.fai -> /data/output/virusbreakend-db_2025-09-09//kraken2-db/library/viral/library.fna.fai
# │       ├── library.fna.masked -> /data/output/virusbreakend-db_2025-09-09//kraken2-db/library/viral/library.fna.masked
# │       ├── manifest.txt -> /data/output/virusbreakend-db_2025-09-09//kraken2-db/library/viral/manifest.txt
# │       └── prelim_map.txt -> /data/output/virusbreakend-db_2025-09-09//kraken2-db/library/viral/prelim_map.txt
# ├── opts.k2d
# ├── seqid2taxid.map
# ├── taxid10239.nbr -> /data/output/virusbreakend-db_2025-09-09//supplementary/taxid10239.nbr
# ├── taxo.k2d
# └── taxonomy
#     ├── accmap.dlflag -> /data/output/virusbreakend-db_2025-09-09//kraken2-db/taxonomy/accmap.dlflag
#     ├── citations.dmp -> /data/output/virusbreakend-db_2025-09-09//kraken2-db/taxonomy/citations.dmp
#     ├── delnodes.dmp -> /data/output/virusbreakend-db_2025-09-09//kraken2-db/taxonomy/delnodes.dmp
#     ├── division.dmp -> /data/output/virusbreakend-db_2025-09-09//kraken2-db/taxonomy/division.dmp
#     ├── gc.prt -> /data/output/virusbreakend-db_2025-09-09//kraken2-db/taxonomy/gc.prt
#     ├── gencode.dmp -> /data/output/virusbreakend-db_2025-09-09//kraken2-db/taxonomy/gencode.dmp
#     ├── images.dmp -> /data/output/virusbreakend-db_2025-09-09//kraken2-db/taxonomy/images.dmp
#     ├── merged.dmp -> /data/output/virusbreakend-db_2025-09-09//kraken2-db/taxonomy/merged.dmp
#     ├── names.dmp -> /data/output/virusbreakend-db_2025-09-09//kraken2-db/taxonomy/names.dmp
#     ├── nodes.dmp -> /data/output/virusbreakend-db_2025-09-09//kraken2-db/taxonomy/nodes.dmp
#     ├── nucl_gb.accession2taxid -> /data/output/virusbreakend-db_2025-09-09//kraken2-db/taxonomy/nucl_gb.accession2taxid
#     ├── nucl_wgs.accession2taxid -> /data/output/virusbreakend-db_2025-09-09//kraken2-db/taxonomy/nucl_wgs.accession2taxid
#     ├── prelim_map.txt
#     ├── readme.txt -> /data/output/virusbreakend-db_2025-09-09//kraken2-db/taxonomy/readme.txt
#     ├── taxdump.dlflag -> /data/output/virusbreakend-db_2025-09-09//kraken2-db/taxonomy/taxdump.dlflag
#     ├── taxdump.tar.gz -> /data/output/virusbreakend-db_2025-09-09//kraken2-db/taxonomy/taxdump.tar.gz
#     └── taxdump.untarflag -> /data/output/virusbreakend-db_2025-09-09//kraken2-db/taxonomy/taxdump.untarflag

## Persist on bucket --------------------------------
# gcloud auth login

## Selectively copy files as not all are required
BUCKET_PATH=gs://hmf-crunch-resources/virus/virusbreakend_db/2025-09-11_fixed-missing-hpv33
gsutil -m cp "$DEST_DB_PATH/virusbreakend-db/*" $BUCKET_PATH
gsutil -m cp "$DEST_DB_PATH/virusbreakend-db/library/viral/*.fna*" "$BUCKET_PATH/library/viral/"
gsutil -m cp "$DEST_DB_PATH/virusbreakend-db/library/added/*.fna*" "$BUCKET_PATH/library/added/"
gsutil -m cp "$DEST_DB_PATH/virusbreakend-db/taxonomy/nodes.dmp" "$BUCKET_PATH/taxonomy/"
gsutil -m cp "$DEST_DB_PATH/virusbreakend-db/taxonomy/taxdump.tar.gz" "$BUCKET_PATH/taxonomy/"
