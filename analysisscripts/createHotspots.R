library(tidyr)
library(dplyr)

pasteNA <- function(...) {
  source = paste(..., sep =",")
  source = gsub("NA,", "", source)
  source = gsub(",NA", "", source)
  return (source)
}

mnv_to_svn <- function(mnvs) {
  snvs = data.frame(chromosome = character(), position = integer(), ref = character(), alt = character(), stringsAsFactors = F)
  for (i in 1:nrow(mnvs)) {
    
    mnvEntry = mnvs[i, ]
    j = nchar(mnvEntry$ref)
    
    chromosome = rep(mnvEntry$chromosome, j)
    position = mnvEntry$position + c(0:(j-1))
    ref = strsplit(mnvEntry$ref, "")[[1]]
    alt = strsplit(mnvEntry$alt, "")[[1]]
    
    snvs = rbind(snvs, data.frame(chromosome, position, alt, ref, stringsAsFactors = F))
  }
  
  return (snvs)
}

chromosomeLevels = c(1:22, "X","Y")

cgi = read.csv("~/Downloads/cgi_variant_list", sep ="\t", stringsAsFactors = F)
cgi = cgi %>% 
  select(chromosome = CHROMOSOME, position = POSITION, ref = REF, alt = ALT, cgiGene = GENE) %>%
  group_by(chromosome, position, ref, alt) %>% distinct(chromosome, position, ref, alt)

onco = read.csv("~/Downloads/oncoKb_variant_list", sep ="\t", stringsAsFactors = F)
onco = onco %>% 
  filter(INFO %in% c("Likely Oncogenic", "Oncogenic")) %>%
  select(chromosome = CHROMOSOME, position = POSITION, ref = REF, alt = ALT, oncoGene = GENE) %>%
  group_by(chromosome, position, ref, alt) %>% distinct(chromosome, position, ref, alt)

civic = read.csv("~/Downloads/civic_variant_list", sep ="\t", stringsAsFactors = F)
civic = civic %>% 
  filter(INFO == "true") %>%
  select(chromosome = CHROMOSOME, position = POSITION, ref = REF, alt = ALT, civicGene = GENE) %>%
  group_by(chromosome, position, ref, alt) %>% distinct(chromosome, position, ref, alt) 

tert = data.frame(chromosome = "5", position = 1295250, ref = "G", alt = "A", stringsAsFactors = F)

hotspots = merge(cgi, onco, by = c("chromosome", "position", "ref", "alt"), all = T)
hotspots = merge(hotspots, civic, by = c("chromosome", "position", "ref", "alt"), all = T)
hotspots = merge(hotspots, tert, by = c("chromosome", "position", "ref", "alt"), all = T)

mnvInd = nchar(hotspots$ref) > 1 & nchar(hotspots$ref) == nchar(hotspots$alt)
mnvs = hotspots[mnvInd, ]
snvs = mnv_to_svn(mnvs)

hotspots = hotspots[!mnvInd, ]
hotspots = rbind(hotspots, snvs)
hotspots = hotspots %>% distinct(chromosome, position, ref, alt) %>%
  mutate(chromosome = factor(chromosome, chromosomeLevels)) %>%
  arrange(chromosome, position)

rm(cgi,civic,onco,tert, mnvs, mnvInd, snvs)
write.table(hotspots, file = "~/hmf/resources/Hotspot.tsv", sep = "\t", col.names = F, quote = F, row.names = F)

#bgzip Hotspot.tsv.gz
#tabix -f -s 1 -b 2 -e 2 -S 0 Hotspot.tsv.gz
#tabix Hotspot.tsv.gz 1:9780851-9780853

#Apply hotspots to somatic and germine pon
#/data/common/tools/bcftools_v1.3/bcftools annotate -a ./Hotspot.tsv.gz -c CHROM,POS,REF,ALT -m +HOTSPOT ./SOMATIC_PON.vcf.gz -o ./SOMATIC_PON_HOTSPOT.vcf
#/data/common/tools/bcftools_v1.3/bcftools filter -i 'HOTSPOT=1' ./SOMATIC_PON_HOTSPOT.vcf -o SOMATIC_PON_HOTSPOT_FILTERED.vcf
#/data/common/tools/bcftools_v1.3/bcftools annotate -a ./Hotspot.tsv.gz -c CHROM,POS,REF,ALT -m +HOTSPOT ./GERMLINE_PON.vcf.gz -o ./GERMLINE_PON_HOTSPOT.vcf
#/data/common/tools/bcftools_v1.3/bcftools filter -i 'HOTSPOT=1' ./GERMLINE_PON_HOTSPOT.vcf -o GERMLINE_PON_HOTSPOT_FILTERED.vcf

germline = read.table("~/hmf/analysis/hotspot/GERMLINE_PON_HOTSPOT_FILTERED.vcf",header=TRUE,sep='\t', stringsAsFactors = F)
somatic = read.table("~/hmf/analysis/hotspot/SOMATIC_PON_HOTSPOT_FILTERED.vcf",header=TRUE,sep='\t', stringsAsFactors = F)
combinedPon = merge(germline,somatic,by=c('CHROM','POS','ID','REF','ALT','FILTER','QUAL'),all=TRUE)  %>% 
  separate(INFO.x,into=c('GTotal','a','b','c'),sep=';') %>% separate(GTotal,into=c('d','GPONCount'),sep='=') %>% 
  separate(INFO.y,into=c('STotal','g','e','f'),sep=';') %>% separate(STotal,into=c('h','SPONCount'),sep='=') %>% 
  mutate(GPONCount=as.numeric(GPONCount)) %>% 
  mutate(SPONCount=as.numeric(SPONCount)) %>%
  select(chromosome = CHROM,position = POS,ref = REF,alt = ALT,GPONCount,SPONCount)
combinedPon$chromosome <- as.character(combinedPon$chromosome)
rm(germline, somatic)

hotspots = merge(hotspots, combinedPon, c("chromosome", "position", "ref", "alt"), all = T)
hotspots[is.na(hotspots)] <- 0
hotspots = hotspots %>% filter(GPONCount <= 2) %>% select(chromosome, position, ref, alt)
write.table(hotspots, file = "~/hmf/resources/Hotspot.tsv", sep = "\t", col.names = F, quote = F, row.names = F)





