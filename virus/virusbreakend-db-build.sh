# ## Download kraken2
# git clone https://github.com/DerrickWood/kraken2/

## Copy tools to lnguyen@dev-syd-testing ================================
## Run locally
gcloud compute ssh --zone "europe-west4-a" "dev-syd-testing" --tunnel-through-iap --project "hmf-crunch"

##
sudo gcloud compute scp --recurse lnguyen@data-vm-prod-2:/data/tools/samtools/ /data/tools/
sudo gcloud compute scp --recurse lnguyen@data-vm-prod-2:/data/tools/gridss/ /data/tools/
sudo gcloud compute scp --recurse lnguyen@data-vm-prod-2:/data/tools/ncbi-blast/ /data/tools/
sudo gcloud compute scp --recurse lnguyen@data-vm-prod-2:/data/tools/bcftools/ /data/tools/
sudo gcloud compute scp --recurse lnguyen@data-vm-prod-2:/data/tools/repeatmasker/ /data/tools/
sudo gcloud compute scp --recurse lnguyen@data-vm-prod-2:/data/tools/bwa/ /data/tools/
sudo gcloud compute scp --recurse lnguyen@data-vm-prod-2:/data/experiments/virus/20240611_run_virus_breakend/tools/kraken2/ /data/tools/
sudo gcloud compute scp --recurse lnguyen@data-vm-prod-2:/data/experiments/virus/20240611_run_virus_breakend/tools/edirect/ /data/tools/

sudo gcloud compute scp lnguyen@data-vm-prod-2:/data/resources/bucket/virus_reference_genome/human_virus.fa /data/resources/bucket/virus_reference_genome/
sudo gcloud compute scp --recurse lnguyen@data-vm-prod-2:/data/benchmark_samples/colo_mini_v37 /data/samples/

sudo gcloud compute scp --recurse lnguyen@data-vm-prod-2:/data/experiments/virus/20240611_run_virus_breakend/ /data/experiments/virus/


# ## Install NCBI edirect 
# sudo sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)" ## installs to /root/edirect/
# sudo mv /root/edirect/ /data/experiments/virus/20240611_run_virus_breakend/tools/

## !!! Need to edit rsync_from_ncbi.pl: https://www.bpetersen.dk/post/kraken2-rsync_from_ncbi-pl-unexpected-ftp-path-new-server-for
## Change: if (! ($full_path =~ s#^ftp://${qm_server}${qm_server_path}/##))
## To: if (! ($full_path =~ s#^https://${qm_server}${qm_server_path}/##))

## Modified virusbreakend-build bash script to build reduced database ================================
## Run
export PATH=$PATH:/data/tools/edirect/:/data/tools/ncbi-blast/2.14.1/ncbi-blast/bin/:/data/tools/samtools/1.20/:/data/tools/kraken2/2.1.0/:/data/tools/gridss/2.13.2-modified/:/data/tools/bcftools/1.9/:/data/tools/bwa/0.7.17/:/data/tools/repeatmasker/4.1.2/RepeatMasker/

DEST_DB_PATH=/data/experiments/virus/20240611_run_virus_breakend/virusbreakend-db/virus-human/

## Download kraken databases
kraken2-build --download-taxonomy --db $DEST_DB_PATH
kraken2-build --download-library viral --db $DEST_DB_PATH
kraken2-build --download-library human --db $DEST_DB_PATH

## Download neighbours sequences. !THIS TAKES AGES (>2h)!
esearch -db nuccore -query "Viruses[Organism] NOT cellular organisms[ORGN] NOT wgs[PROP] NOT AC_000001:AC_999999[pacc] NOT gbdiv syn[prop] AND (srcdb_refseq[PROP] OR nuccore genome samespecies[Filter])" \
| efetch -format fasta > $DEST_DB_PATH/neighbour.fa

## Add neighbours to 
kraken2-build --add-to-library $DEST_DB_PATH/neighbour.fa --db $DEST_DB_PATH

# Issue with neighbours: serotypes (e.g. HPV-45) are listed as neighbours of the species taxid
# Use Neighbours for host taxid filtering. TODO: get a better source
wget --output-document=$DEST_DB_PATH/taxid10239.nbr "https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=10239&cmd=download2"

kraken2-build --build --db $DEST_DB_PATH

for f in $(find $DEST_DB_PATH/ -name '*.fna') ; do
	samtools faidx $f
	rm -f $f.dict
	java -cp /data/tools/gridss/2.13.2-modified/gridss.jar \
		picard.cmdline.PicardCommandLine \
		CreateSequenceDictionary \
		R=$f \
		O=$f.dict
done
