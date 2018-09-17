#install.packages("vcfR")
# scp -C jon@hmf-datastore:/data/experiments/180815_pcawg_on_HMF_pipeline_validation/PCAWG_HMF/*/*somaticSV_bpi.vcf.gz /Users/jon/hmf/analysis/pcawg/ours/
# scp -C jon@hmf-datastore:/data/experiments/180815_pcawg_on_HMF_pipeline_validation/PCAWG_HMF_set/*.vcf.gz /Users/jon/hmf/analysis/pcawg/theirs/

######### MAPPING FILES
library(vcfR)

hmfSvVcfs = c(
  "EGAZ00001225530_EGAZ00001225588_somaticSV_bpi.vcf.gz",
  "f7f3e156-0ddb-72b9-e040-11ac0d48542c_f7f3e156-0dde-72b9-e040-11ac0d48542c_somaticSV_bpi.vcf.gz",
  "fc8130df-1cda-cade-e040-11ac0d485dec_fc8130df-1cdd-cade-e040-11ac0d485dec_somaticSV_bpi.vcf.gz",
  "EGAZ00001225554_EGAZ00001225511_somaticSV_bpi.vcf.gz",
  "6ff52f43-aae0-432d-ae31-374c22c45fd9_d392ded3-afc8-4c79-b278-40245f18f2f8_somaticSV_bpi.vcf.gz",
  "fc8130df-1f1e-c8f9-e040-11ac0d485dfc_fc8130df-1f21-c8f9-e040-11ac0d485dfc_somaticSV_bpi.vcf.gz",
  "EGAZ00001222819_EGAZ00001222815_somaticSV_bpi.vcf.gz",
  "EGAZ00001225550_EGAZ00001225531_somaticSV_bpi.vcf.gz",
  "1abec81b-9e0f-40a6-9840-967ee0ab8b54_e84debc4-b47d-48ed-a0d0-2859f0ebf987_somaticSV_bpi.vcf.gz",
  "EGAZ00001225631_EGAZ00001225528_somaticSV_bpi.vcf.gz",
  "EGAZ00001222797_EGAZ00001224455_somaticSV_bpi.vcf.gz",
  "EGAZ00001222825_EGAZ00001222794_somaticSV_bpi.vcf.gz",
  "d1b992ab-41b0-4d39-992f-c4f8da898576_b75b2663-dcc6-411c-bfcc-574aa33cf388_somaticSV_bpi.vcf.gz",
  "EGAZ00001222785_EGAZ00001222839_somaticSV_bpi.vcf.gz",
  "b5c7e4d8-19bb-4253-adeb-abaf1b31aa5d_2a8d63eb-0174-4213-9214-413f391f512c_somaticSV_bpi.vcf.gz",
  "6b38c35a-7b39-4e96-9237-b02206c5bcc6_f26b1f44-12de-43ba-85bb-bc61741a5a88_somaticSV_bpi.vcf.gz",
  "dd48a2fc-7c5b-44b5-bf3c-ffc8350445da_8b28f6d2-4b7d-493b-826e-b119a4fb0cb4_somaticSV_bpi.vcf.gz",
  "EGAZ00001224456_EGAZ00001222832_somaticSV_bpi.vcf.gz",
  "fc8130df-685d-7677-e040-11ac0d485ddc_fc8130df-6860-7677-e040-11ac0d485ddc_somaticSV_bpi.vcf.gz",
  "1319bfce-ea65-4181-874d-b0f4c4d4a789_0d0793c1-df1b-4db1-ba36-adcb960cc0f5_somaticSV_bpi.vcf.gz",
  "6a7fdabe-fd87-481c-b9a3-1f46ad3097ef_c2ec7f57-8510-4bbf-a2e9-dbd9ce8dcad1_somaticSV_bpi.vcf.gz",
  "EGAZ00001225621_EGAZ00001225598_somaticSV_bpi.vcf.gz",
  "EGAZ00001222791_EGAZ00001222814_somaticSV_bpi.vcf.gz",
  "b589d2a2-5bf9-4080-bc99-0ae984dbae43_bbe59385-5f83-43f6-a485-517c860bef6f_somaticSV_bpi.vcf.gz")

pcawgSvVcfs = c(
  "fc8130df-18fe-c74d-e040-11ac0d485df2.pcawg_consensus_1.6.161116.somatic.sv.vcf.gz",
  "f7f3e156-0dde-72b9-e040-11ac0d48542c.pcawg_consensus_1.6.161116.somatic.sv.vcf.gz",
  "fc8130df-1cdd-cade-e040-11ac0d485dec.pcawg_consensus_1.6.161116.somatic.sv.vcf.gz",
  "fc8130df-e399-e34d-e040-11ac0c483279.pcawg_consensus_1.6.161116.somatic.sv.vcf.gz",
  "d392ded3-afc8-4c79-b278-40245f18f2f8.pcawg_consensus_1.6.161116.somatic.sv.vcf.gz",
  "fc8130df-1f21-c8f9-e040-11ac0d485dfc.pcawg_consensus_1.6.161116.somatic.sv.vcf.gz",
  "fc68c24d-47ad-7961-e040-11ac0c48595c.pcawg_consensus_1.6.161116.somatic.sv.vcf.gz",
  "fc8130df-6bec-7627-e040-11ac0d485e04.pcawg_consensus_1.6.161116.somatic.sv.vcf.gz",
  "e84debc4-b47d-48ed-a0d0-2859f0ebf987.pcawg_consensus_1.6.161116.somatic.sv.vcf.gz",
  "fc8130e0-0f1a-b6eb-e040-11ac0c48328f.pcawg_consensus_1.6.161116.somatic.sv.vcf.gz",
  "f221cbb5-eefa-187f-e040-11ac0c481708.pcawg_consensus_1.6.161116.somatic.sv.vcf.gz",
  "fc63cbab-d27a-5ebb-e040-11ac0c48724f.pcawg_consensus_1.6.161116.somatic.sv.vcf.gz",
  "b75b2663-dcc6-411c-bfcc-574aa33cf388.pcawg_consensus_1.6.161116.somatic.sv.vcf.gz",
  "fc447d4f-2532-c8ea-e040-11ac0c48469f.pcawg_consensus_1.6.161116.somatic.sv.vcf.gz",
  "2a8d63eb-0174-4213-9214-413f391f512c.pcawg_consensus_1.6.161116.somatic.sv.vcf.gz",
  "f26b1f44-12de-43ba-85bb-bc61741a5a88.pcawg_consensus_1.6.161116.somatic.sv.vcf.gz",
  "8b28f6d2-4b7d-493b-826e-b119a4fb0cb4.pcawg_consensus_1.6.161116.somatic.sv.vcf.gz",
  "f393bb07-270c-2c93-e040-11ac0d484533.pcawg_consensus_1.6.161116.somatic.sv.vcf.gz",
  "fc8130df-6860-7677-e040-11ac0d485ddc.pcawg_consensus_1.6.161116.somatic.sv.vcf.gz",
  "0d0793c1-df1b-4db1-ba36-adcb960cc0f5.pcawg_consensus_1.6.161116.somatic.sv.vcf.gz",
  "c2ec7f57-8510-4bbf-a2e9-dbd9ce8dcad1.pcawg_consensus_1.6.161116.somatic.sv.vcf.gz",
  "fc8130df-8ec8-5b1e-e040-11ac0d485e06.pcawg_consensus_1.6.161116.somatic.sv.vcf.gz",
  "f393bb08-4121-cad8-e040-11ac0d484535.pcawg_consensus_1.6.161116.somatic.sv.vcf.gz",
  "bbe59385-5f83-43f6-a485-517c860bef6f.pcawg_consensus_1.6.161116.somatic.sv.vcf.gz")

vcfFilenames = data.frame(hmfSv = hmfSvVcfs, pcawgSv =pcawgSvVcfs, stringsAsFactors = F)

load_sv_file <- function(directory, fileName) {
  chromFactors = c(1:22, "X","Y")
  
  vcf <- read.vcfR(paste0(directory, fileName), verbose = F)
  vcf  = vcfR2tidy(vcf, single_frame = F)$fix %>% filter(!is.na(MATEID), FILTER %in% c(".", "PASS")) %>% 
    mutate(CHROM = factor(CHROM, chromFactors)) %>%
    arrange(CHROM, POS)
  
  vcf.id = vcf %>% select(ID) %>% mutate(IDRank = row_number())
  vcf.mateid = vcf %>% select(ID = MATEID) %>% mutate(MATEIDRank = row_number())  
  
  vcf = left_join(vcf, vcf.id, by = "ID") %>% 
    left_join(vcf.mateid, by = "ID") %>% 
    mutate(primary = IDRank < MATEIDRank, CHROM = as.character(CHROM)) %>% 
    select(-IDRank, -MATEIDRank) 
  
  vcf.start = vcf
  vcf.end = vcf
  vcf.joined = inner_join(vcf.start, vcf.end, by = c("ID" = "MATEID"), suffix = c("_start", '_end')) %>% 
    filter(primary_start) %>% 
    select(-primary_start, -primary_end) %>%
    select(startChromosome = CHROM_start, endChromosome = CHROM_end, startPosition = POS_start, endPosition = POS_end, everything()) 
  
  return (vcf.joined)
}

sv_overlaps <- function(query, subject, maxgap = -1) {
  require(tidyr)
  require(dplyr)
  require(GenomicRanges)
  
  queryStartRange <- GRanges(query$startChromosome, IRanges(query$startPosition, query$startPosition))
  subjectStartRange <- GRanges(subject$startChromosome, IRanges(subject$startPosition, subject$startPosition))
  startOverlaps = data.frame(findOverlaps(queryStartRange, subjectStartRange, type="any", select="all", maxgap = maxgap))
  
  queryEndRange <- GRanges(query$endChromosome, IRanges(query$endPosition, query$endPosition))
  subjectEndRange <- GRanges(subject$endChromosome, IRanges(subject$endPosition, subject$endPosition))
  endOverlaps = data.frame(findOverlaps(queryEndRange, subjectEndRange, type="any", select="all", maxgap = maxgap))
  
  overlaps = inner_join(startOverlaps, endOverlaps, by = c("queryHits", "subjectHits"))
  
  overlapQueryData = query[overlaps$queryHits, ] %>%
    mutate(queryHits = overlaps$queryHits) %>%
    select(queryHits, sampleId, startChromosome, endChromosome, startPosition, endPosition, startOrientation, endOrientation)
  
  overlapSubjectData = subject[overlaps$subjectHits, ] %>%
    mutate(subjectHits = overlaps$subjectHits) %>%
    select(subjectHits, subjectStartPosition = startPosition, subjectEndPosition = endPosition, subjectStartOrientation = startOrientation, subjectEndOrientation = endOrientation)
  
  overlapsData = bind_cols(overlapQueryData, overlapSubjectData) %>%
    filter(startOrientation == subjectStartOrientation, endOrientation == subjectEndOrientation) %>%
    select(-subjectStartOrientation, -subjectEndOrientation) %>%
    mutate(startPositionDiff = abs(startPosition - subjectStartPosition), endPositionDiff = abs(endPosition - subjectEndPosition), positionDiff = startPositionDiff + endPositionDiff) %>%
    group_by(startChromosome, endChromosome, startPosition, endPosition,startOrientation,endOrientation) %>%
    top_n(1, -positionDiff) %>%
    group_by(queryHits) %>%
    top_n(1, -subjectHits)
  
  return (overlapsData %>% select(queryHits, subjectHits))
}


##### EXECUTE ALGORITM
hmfDirectory = "/Users/jon/hmf/analysis/pcawg/ours/"
pcawgDirectory = "/Users/jon/hmf/analysis/pcawg/theirs/"

query = "select * from structuralVariant where sampleId not like 'CPCT%'"
dbProd = dbConnect(MySQL(), dbname='landscape_paper_validation', groups="RAnalysis")
hmfSVs = dbGetQuery(dbProd, query)
dbDisconnect(dbProd)
rm(dbProd, query)

combinedSV = data.frame()

for (i in c(1:nrow(vcfFilenames))) {
  cat("Processing", i, "of 24", "\n")

  sample = gsub("_somaticSV_bpi.vcf.gz", "", vcfFilenames[i, "hmfSv"])
  sample = gsub(".*_", "", sample)
  
  hmfSv = hmfSVs %>% 
    filter(
      sampleId == sample, 
      filter == 'PASS') %>% 
    mutate(source = "hmf", scope = "Private") 
  
  pcawgVCFFile = vcfFilenames[i, "pcawgSv"]
  pcawgSv = load_sv_file(pcawgDirectory, pcawgVCFFile) %>%
    mutate(
      startOrientation = ifelse(STRAND_start == "+", 1, -1),
      endOrientation = ifelse(STRAND_end == "+", 1, -1)) %>% 
    mutate(source = "pcawg", scope = "Private")
  
  overlaps = sv_overlaps(hmfSv, pcawgSv, maxgap = 100)
  hmfSv[overlaps$queryHits, "scope"] <- "Shared"
  pcawgSv[overlaps$subjectHits, "scope"] <- "Shared"
  
  combined = bind_rows(pcawgSv, hmfSv) %>%
    mutate(sampleId = sample)
  
  combinedSV = bind_rows(combinedSV, combined)
}

#################### PRE FILTER STATS
combinedSV %>% filter(source == 'hmf') %>% mutate(exclude = type == 'BND' & ploidy < 0.2) %>% group_by(exclude) %>% count() %>% ungroup() %>% mutate(percent = n / sum(n))
combinedSV %>% filter(source == 'hmf', ploidy < 0.2) %>% group_by(type) %>% count() %>% ungroup() %>% mutate(percent = n / sum(n))

### FILTER
combinedSV = combinedSV %>% 
  filter(
    !(source == "hmf" & type == 'BND' & ploidy < 0.2), 
    startChromosome != endChromosome | endPosition-startPosition>1000)

#################### POST FILTER STATS
combinedSV %>% group_by(source) %>% count()
combinedSV %>% filter(source == 'pcawg') %>% group_by(scope) %>% count() %>% ungroup() %>% mutate(percent = n / sum(n))
combinedSV %>% filter(source == 'hmf') %>% group_by(scope) %>% count() %>% ungroup() %>% mutate(percent = n / sum(n))



save(combinedSV, file = "/Users/jon/hmf/analysis/pcawg/combinedSV.RData")


