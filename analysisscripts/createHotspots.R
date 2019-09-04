#!/usr/bin/Rscript

#### This file used to create the hotspot list for our HMF pipeline.

library(tidyr)
library(dplyr)
library(RMySQL)

pasteNA <- function(...) {
  source = paste(..., sep =",")
  source = gsub("NA,", "", source)
  source = gsub(",NA", "", source)
  return (source)
}

chromosomeLevels = c(1:22, "X","Y")

cgi = read.csv("/data/experiments/compare_hotspots/prod_hotspots/cgi_variant_list", sep = "\t", stringsAsFactors = F)
cgi = cgi %>% 
  select(chromosome = CHROMOSOME, position = POSITION, ref = REF, alt = ALT, cgiGene = GENE) %>%
  group_by(chromosome, position, ref, alt) %>% distinct(chromosome, position, ref, alt)

onco = read.csv("/data/experiments/compare_hotspots/prod_hotspots/oncoKb_variant_list", sep = "\t", stringsAsFactors = F)
onco = onco %>% 
  filter(INFO %in% c("Likely Oncogenic", "Oncogenic")) %>%
  select(chromosome = CHROMOSOME, position = POSITION, ref = REF, alt = ALT, oncoGene = GENE) %>%
  group_by(chromosome, position, ref, alt) %>% distinct(chromosome, position, ref, alt)

civic = read.csv("/data/experiments/compare_hotspots/prod_hotspots/civic_variant_list", sep = "\t", stringsAsFactors = F)
civic = civic %>% 
  filter(INFO == "true") %>%
  select(chromosome = CHROMOSOME, position = POSITION, ref = REF, alt = ALT, civicGene = GENE) %>%
    group_by(chromosome, position, ref, alt) %>%
    distinct(chromosome, position, ref, alt)

dbProd = dbConnect(MySQL(), user = db_user, password = db_password, dbname = db_name, groups = "RAnalysis")
query = "SELECT distinct chromosome, position, ref, alt FROM somaticVariant where gene = 'TERT' and position in (1295242,1295228,1295250) AND canonicalEffect = 'upstream gene variant'"
tert = dbGetQuery(dbProd, query)
dbDisconnect(dbProd)

hotspots = merge(cgi, onco, by = c("chromosome", "position", "ref", "alt"), all = T)
hotspots = merge(hotspots, civic, by = c("chromosome", "position", "ref", "alt"), all = T)
hotspots = merge(hotspots, tert, by = c("chromosome", "position", "ref", "alt"), all = T)

hotspots = hotspots %>% distinct(chromosome, position, ref, alt) %>%
  mutate(chromosome = factor(chromosome, chromosomeLevels)) %>%
  arrange(chromosome, position)

write.table(hotspots, file = "/data/experiments/compare_hotspots/prod_hotspots/Hotspot_prod.tsv", sep = "\t", col.names = F, quote = F, row.names = F)

#bgzip Hotspot.tsv
#tabix -f -s 1 -b 2 -e 2 -S 0 Hotspot.tsv.gz
#tabix Hotspot.tsv.gz 1:9780851-9780853

#Apply hotspots to somatic and germine pon
#/data/common/tools/bcftools_v1.3/bcftools annotate -a ./Hotspot.tsv.gz -c CHROM,POS,REF,ALT -m +HOTSPOT ./SOMATIC_PON.vcf.gz -o ./SOMATIC_PON_HOTSPOT.vcf
#/data/common/tools/bcftools_v1.3/bcftools filter -i 'HOTSPOT=1' ./SOMATIC_PON_HOTSPOT.vcf -o SOMATIC_PON_HOTSPOT_FILTERED.vcf
#/data/common/tools/bcftools_v1.3/bcftools annotate -a ./Hotspot.tsv.gz -c CHROM,POS,REF,ALT -m +HOTSPOT ./GERMLINE_PON.vcf.gz -o ./GERMLINE_PON_HOTSPOT.vcf
#/data/common/tools/bcftools_v1.3/bcftools filter -i 'HOTSPOT=1' ./GERMLINE_PON_HOTSPOT.vcf -o GERMLINE_PON_HOTSPOT_FILTERED.vcf

#germline = read.table("~/hmf/analysis/hotspot/GERMLINE_PON_HOTSPOT_FILTERED.vcf",header=TRUE,sep='\t', stringsAsFactors = F)
#somatic = read.table("~/hmf/analysis/hotspot/SOMATIC_PON_HOTSPOT_FILTERED.vcf",header=TRUE,sep='\t', stringsAsFactors = F)
#combinedPon = merge(germline,somatic,by=c('CHROM','POS','ID','REF','ALT','FILTER','QUAL'),all=TRUE)  %>%
#  separate(INFO.x,into=c('GTotal','a','b','c'),sep=';') %>% separate(GTotal,into=c('d','GPONCount'),sep='=') %>%
#  separate(INFO.y,into=c('STotal','g','e','f'),sep=';') %>% separate(STotal,into=c('h','SPONCount'),sep='=') %>%
#  mutate(GPONCount=as.numeric(GPONCount)) %>%
#  mutate(SPONCount=as.numeric(SPONCount)) %>%
#  select(chromosome = CHROM,position = POS,ref = REF,alt = ALT,GPONCount,SPONCount)
#combinedPon$chromosome <- as.character(combinedPon$chromosome)
#rm(germline, somatic)

#hotspots = merge(hotspots, combinedPon, c("chromosome", "position", "ref", "alt"), all = T)
#hotspots[is.na(hotspots)] <- 0
#hotspots = hotspots %>% filter(GPONCount <= 2) %>% select(chromosome, position, ref, alt)
#write.table(hotspots, file = "~/hmf/resources/Hotspot.tsv", sep = "\t", col.names = F, quote = F, row.names = F)





