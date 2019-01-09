library(vcfR)
library(dplyr)
library(RMySQL)


read_sage_unfiltered <- function(filename) {
  vcf <- read.vcfR(filename, verbose = F) 
  Z <- vcfR2tidy(vcf, single_frame = T)$dat %>% filter(HOTSPOT=='known', FILTER %in% c('PASS', 'GERMLINE_INDEL'))
  
  Z %>% spread(Indiv, gt_AD)
  
  normal = Z %>% filter(substr(Indiv, 13, 13) == 'R') %>% 
    dplyr::select(-Indiv, -gt_GT, -gt_GT_alleles, -HOTSPOT, -GHBL, -ID, -QUAL, -FILTER, -AF)
  colnames(normal) <- c("chromosome","position", "ref", "alt", "AD", "DP")
  
  tumor = Z %>% filter(substr(Indiv, 13, 13) == 'T')  %>% 
    dplyr::select(-gt_GT, -gt_GT_alleles, -HOTSPOT, -GHBL, -ID)
  colnames(tumor) <- c("chromosome","position", "ref", "alt", "qual", "filter", "AF", "sampleId", "AD", "DP")
  
  combined = normal %>% left_join(tumor, by = c("chromosome", "position","ref","alt"), suffix = c(".normal",".tumor")) %>%
    dplyr::select(sampleId, everything())
  
  return (combined)
}

read_sage_filtered <- function(filename, filter = 'SAGE_PON') {
  vcf <- read.vcfR(filename, verbose = F)
  Z <- vcfR2tidy(vcf, single_frame = T)$dat %>% filter(FILTER == filter)
  
  tumor = Z %>% filter(substr(Indiv, 13, 13) == 'T')  %>% dplyr::select(Indiv, CHROM, POS, REF, ALT, FILTER)
  colnames(tumor) <- c("sampleId", "chromosome","position", "ref", "alt", "filter")
  return (tumor)
}

####### PREPARE SAGE DATA
sageUnfiltered = data.frame()
for (file in list.files(path = "/Users/jon/hmf/analysis/sage/", pattern = "*unfiltered.vcf")) {
  sageUnfiltered = bind_rows(sageUnfiltered, read_sage_unfiltered(paste0("/Users/jon/hmf/analysis/sage/", file)))
}

ponFiltered = data.frame()
for (file in list.files(path = "/Users/jon/hmf/analysis/sage/", pattern = "*hotspot.vcf")) {
  ponFiltered = bind_rows(ponFiltered, read_sage_filtered(paste0("/Users/jon/hmf/analysis/sage/", file)))
}

load(file = "hmf/RData/Processed/highestPurityCohortSummary.RData")
sageSamples =  gsub("\\.hotspot\\.unfiltered\\.vcf", "", list.files(path = "/Users/jon/hmf/analysis/sage/", pattern = "*unfiltered.vcf"))
sageSamplesInCohort = data.frame(sampleId = sageSamples, stringsAsFactors = F) %>% filter(sampleId %in% highestPurityCohortSummary$sampleId)
germlineIndel = sageUnfiltered %>% filter(filter == 'GERMLINE_INDEL', sampleId %in% highestPurityCohortSummary$sampleId)
sageHotspots = sageUnfiltered %>% filter(filter == 'PASS',sampleId %in% highestPurityCohortSummary$sampleId)
save(sageUnfiltered, germlineIndel, sageHotspots, sageSamplesInCohort, ponFiltered, file = "/Users/jon/hmf/analysis/sage/sageHotspots.RData")


##### PREPARE DATABASE DATA
load(file = "/Users/jon/hmf/analysis/sage/sageHotspots.RData")

dbProd = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis")
sampleIdString = paste("'", sageSamplesInCohort$sampleId, "'", collapse = ",", sep = "")
query = paste0("SELECT * from somaticVariant where hotspot='HOTSPOT' and sampleId in (", sampleIdString, ")")
databaseHotspots = dbGetQuery(dbProd, query)
save(databaseHotspots, file = "/Users/jon/hmf/analysis/sage/databaseHotspots.RData")
dbDisconnect(dbProd)
rm(dbProd)


##### PROCESS
load(file = "/Users/jon/hmf/analysis/sage/sageHotspots.RData")
load(file = "/Users/jon/hmf/analysis/sage/databaseHotspots.RData")

germlineIndel = germlineIndel %>% dplyr::select(sampleId, chromosome, position, ref, alt, filter)

sage = left_join(sageHotspots, ponFiltered, by = c("sampleId", "chromosome", "position", "ref", "alt")) %>% 
  filter(is.na(filter.y)) %>%
  dplyr::select(-filter.y)

combined = databaseHotspots %>% filter(filter == 'PASS') %>% left_join(sage, by = c("sampleId", "chromosome","position", "ref", "alt"))

missingFromDatabase = right_join(databaseHotspots, sage, by = c("sampleId", "chromosome","position", "ref", "alt")) %>% filter(is.na(id))
missingFromSage = combined %>% filter(is.na(AF)) %>% 
  left_join(ponFiltered, by = c("sampleId", "chromosome","position", "ref", "alt")) %>%
  left_join(germlineIndel, by = c("sampleId", "chromosome","position", "ref", "alt"))

combinedDatabaseAndSage = combined %>% filter(!is.na(AF)) 

load(file = "~/hmf/RData/Processed/hpcDriversByGene.RData")
analysisSet = merge(combinedDatabaseAndSage,hpcDriversByGene %>% distinct(gene), by='gene') %>% 
  filter(repeatCount<8,DP.tumor<300,canonicalCodingEffect %in% c('MISSENSE','NONSENSE_OR_FRAMESHIFT')) %>% 
  separate(AD.tumor,c('AD.tumor.REF','AD.tumor.ALT'),sep=',') %>% 
  mutate(reasonMissed=ifelse(recovered==1,ifelse(as.numeric(AD.tumor.ALT)<=6,'LowTumorAltCount',ifelse(DP.normal<=17,'LowCoverageRef','Unknown')),'Found'))

View(analysisSet %>% group_by(recovered) %>% dplyr::count())
View(analysisSet %>% group_by(reasonMissed) %>% dplyr::count())

