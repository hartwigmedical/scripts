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


write.table(hotspots, file = "~/hmf/resources/Hotspot.tsv", sep = "\t", col.names = F, quote = F, row.names = F)

#bgzip Hotspot.tsv.gz
#tabix -f -s 1 -b 2 -e 2 -S 0 Hotspot.tsv.gz
#tabix Hotspot.tsv.gz 1:9780851-9780853






