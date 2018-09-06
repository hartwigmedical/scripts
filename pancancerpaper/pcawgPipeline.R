#install.packages("vcfR")
# scp -C jon@hmf-datastore:/data/experiments/180815_pcawg_on_HMF_pipeline_validation/PCAWG_HMF/*/*_post_processed.vcf.gz* /Users/jon/hmf/analysis/pcawg/ours/
# scp -C jon@hmf-datastore:/data/experiments/180815_pcawg_on_HMF_pipeline_validation/PCAWG_HMF/*/*somaticSV_bpi.vcf.gz /Users/jon/hmf/analysis/pcawg/ours/
# scp -C jon@hmf-datastore:/data/experiments/180815_pcawg_on_HMF_pipeline_validation/PCAWG_HMF_set/*.vcf.gz /Users/jon/hmf/analysis/pcawg/theirs/

######### MAPPING FILES
library(vcfR)

hmfVcfs = c(
"EGAZ00001225530_EGAZ00001225588_post_processed.vcf.gz",
"f7f3e156-0ddb-72b9-e040-11ac0d48542c_f7f3e156-0dde-72b9-e040-11ac0d48542c_post_processed.vcf.gz",
"fc8130df-1cda-cade-e040-11ac0d485dec_fc8130df-1cdd-cade-e040-11ac0d485dec_post_processed.vcf.gz",
"EGAZ00001225554_EGAZ00001225511_post_processed.vcf.gz",
"6ff52f43-aae0-432d-ae31-374c22c45fd9_d392ded3-afc8-4c79-b278-40245f18f2f8_post_processed.vcf.gz",
"fc8130df-1f1e-c8f9-e040-11ac0d485dfc_fc8130df-1f21-c8f9-e040-11ac0d485dfc_post_processed.vcf.gz",
"EGAZ00001222819_EGAZ00001222815_post_processed.vcf.gz",
"EGAZ00001225550_EGAZ00001225531_post_processed.vcf.gz",
"1abec81b-9e0f-40a6-9840-967ee0ab8b54_e84debc4-b47d-48ed-a0d0-2859f0ebf987_post_processed.vcf.gz",
"EGAZ00001225631_EGAZ00001225528_post_processed.vcf.gz",
"EGAZ00001222797_EGAZ00001224455_post_processed.vcf.gz",
"EGAZ00001222825_EGAZ00001222794_post_processed.vcf.gz",
"d1b992ab-41b0-4d39-992f-c4f8da898576_b75b2663-dcc6-411c-bfcc-574aa33cf388_post_processed.vcf.gz",
"EGAZ00001222785_EGAZ00001222839_post_processed.vcf.gz",
"b5c7e4d8-19bb-4253-adeb-abaf1b31aa5d_2a8d63eb-0174-4213-9214-413f391f512c_post_processed.vcf.gz",
"6b38c35a-7b39-4e96-9237-b02206c5bcc6_f26b1f44-12de-43ba-85bb-bc61741a5a88_post_processed.vcf.gz",
"dd48a2fc-7c5b-44b5-bf3c-ffc8350445da_8b28f6d2-4b7d-493b-826e-b119a4fb0cb4_post_processed.vcf.gz",
"EGAZ00001224456_EGAZ00001222832_post_processed.vcf.gz",
"fc8130df-685d-7677-e040-11ac0d485ddc_fc8130df-6860-7677-e040-11ac0d485ddc_post_processed.vcf.gz",
"1319bfce-ea65-4181-874d-b0f4c4d4a789_0d0793c1-df1b-4db1-ba36-adcb960cc0f5_post_processed.vcf.gz",
"6a7fdabe-fd87-481c-b9a3-1f46ad3097ef_c2ec7f57-8510-4bbf-a2e9-dbd9ce8dcad1_post_processed.vcf.gz",
"EGAZ00001225621_EGAZ00001225598_post_processed.vcf.gz",
"EGAZ00001222791_EGAZ00001222814_post_processed.vcf.gz",
"b589d2a2-5bf9-4080-bc99-0ae984dbae43_bbe59385-5f83-43f6-a485-517c860bef6f_post_processed.vcf.gz")

hmfSampleIds = gsub("_.*", "", hmfVcfs)

pcawgSnvVcfs = c(
"fc8130df-18fe-c74d-e040-11ac0d485df2.consensus.20160830.somatic.snv_mnv.vcf.gz",
"f7f3e156-0dde-72b9-e040-11ac0d48542c.consensus.20160830.somatic.snv_mnv.vcf.gz",
"fc8130df-1cdd-cade-e040-11ac0d485dec.consensus.20160830.somatic.snv_mnv.vcf.gz",
"fc8130df-e399-e34d-e040-11ac0c483279.consensus.20160830.somatic.snv_mnv.vcf.gz",
"d392ded3-afc8-4c79-b278-40245f18f2f8.consensus.20160830.somatic.snv_mnv.vcf.gz",
"fc8130df-1f21-c8f9-e040-11ac0d485dfc.consensus.20160830.somatic.snv_mnv.vcf.gz",
"fc68c24d-47ad-7961-e040-11ac0c48595c.consensus.20160830.somatic.snv_mnv.vcf.gz",
"fc8130df-6bec-7627-e040-11ac0d485e04.consensus.20160830.somatic.snv_mnv.vcf.gz",
"e84debc4-b47d-48ed-a0d0-2859f0ebf987.consensus.20160830.somatic.snv_mnv.vcf.gz",
"fc8130e0-0f1a-b6eb-e040-11ac0c48328f.consensus.20160830.somatic.snv_mnv.vcf.gz",
"f221cbb5-eefa-187f-e040-11ac0c481708.consensus.20160830.somatic.snv_mnv.vcf.gz",
"fc63cbab-d27a-5ebb-e040-11ac0c48724f.consensus.20160830.somatic.snv_mnv.vcf.gz",
"b75b2663-dcc6-411c-bfcc-574aa33cf388.consensus.20160830.somatic.snv_mnv.vcf.gz",
"fc447d4f-2532-c8ea-e040-11ac0c48469f.consensus.20160830.somatic.snv_mnv.vcf.gz",
"2a8d63eb-0174-4213-9214-413f391f512c.consensus.20160830.somatic.snv_mnv.vcf.gz",
"f26b1f44-12de-43ba-85bb-bc61741a5a88.consensus.20160830.somatic.snv_mnv.vcf.gz",
"8b28f6d2-4b7d-493b-826e-b119a4fb0cb4.consensus.20160830.somatic.snv_mnv.vcf.gz",
"f393bb07-270c-2c93-e040-11ac0d484533.consensus.20160830.somatic.snv_mnv.vcf.gz",
"fc8130df-6860-7677-e040-11ac0d485ddc.consensus.20160830.somatic.snv_mnv.vcf.gz",
"0d0793c1-df1b-4db1-ba36-adcb960cc0f5.consensus.20160830.somatic.snv_mnv.vcf.gz",
"c2ec7f57-8510-4bbf-a2e9-dbd9ce8dcad1.consensus.20160830.somatic.snv_mnv.vcf.gz",
"fc8130df-8ec8-5b1e-e040-11ac0d485e06.consensus.20160830.somatic.snv_mnv.vcf.gz",
"f393bb08-4121-cad8-e040-11ac0d484535.consensus.20160830.somatic.snv_mnv.vcf.gz",
"bbe59385-5f83-43f6-a485-517c860bef6f.consensus.20160830.somatic.snv_mnv.vcf.gz")

pcawgIndelVcfs = c(
"fc8130df-18fe-c74d-e040-11ac0d485df2.consensus.20161006.somatic.indel.vcf.gz",
"f7f3e156-0dde-72b9-e040-11ac0d48542c.consensus.20161006.somatic.indel.vcf.gz",
"fc8130df-1cdd-cade-e040-11ac0d485dec.consensus.20161006.somatic.indel.vcf.gz",
"fc8130df-e399-e34d-e040-11ac0c483279.consensus.20161006.somatic.indel.vcf.gz",
"d392ded3-afc8-4c79-b278-40245f18f2f8.consensus.20161006.somatic.indel.vcf.gz",
"fc8130df-1f21-c8f9-e040-11ac0d485dfc.consensus.20161006.somatic.indel.vcf.gz",
"fc68c24d-47ad-7961-e040-11ac0c48595c.consensus.20161006.somatic.indel.vcf.gz",
"fc8130df-6bec-7627-e040-11ac0d485e04.consensus.20161006.somatic.indel.vcf.gz",
"e84debc4-b47d-48ed-a0d0-2859f0ebf987.consensus.20161006.somatic.indel.vcf.gz",
"fc8130e0-0f1a-b6eb-e040-11ac0c48328f.consensus.20161006.somatic.indel.vcf.gz",
"f221cbb5-eefa-187f-e040-11ac0c481708.consensus.20161006.somatic.indel.vcf.gz",
"fc63cbab-d27a-5ebb-e040-11ac0c48724f.consensus.20161006.somatic.indel.vcf.gz",
"b75b2663-dcc6-411c-bfcc-574aa33cf388.consensus.20161006.somatic.indel.vcf.gz",
"fc447d4f-2532-c8ea-e040-11ac0c48469f.consensus.20161006.somatic.indel.vcf.gz",
"2a8d63eb-0174-4213-9214-413f391f512c.consensus.20161006.somatic.indel.vcf.gz",
"f26b1f44-12de-43ba-85bb-bc61741a5a88.consensus.20161006.somatic.indel.vcf.gz",
"8b28f6d2-4b7d-493b-826e-b119a4fb0cb4.consensus.20161006.somatic.indel.vcf.gz",
"f393bb07-270c-2c93-e040-11ac0d484533.consensus.20161006.somatic.indel.vcf.gz",
"fc8130df-6860-7677-e040-11ac0d485ddc.consensus.20161006.somatic.indel.vcf.gz",
"0d0793c1-df1b-4db1-ba36-adcb960cc0f5.consensus.20161006.somatic.indel.vcf.gz",
"c2ec7f57-8510-4bbf-a2e9-dbd9ce8dcad1.consensus.20161006.somatic.indel.vcf.gz",
"fc8130df-8ec8-5b1e-e040-11ac0d485e06.consensus.20161006.somatic.indel.vcf.gz",
"f393bb08-4121-cad8-e040-11ac0d484535.consensus.20161006.somatic.indel.vcf.gz",
"bbe59385-5f83-43f6-a485-517c860bef6f.consensus.20161006.somatic.indel.vcf.gz")

vcfFilenames = data.frame(hmf = hmfVcfs, pcawfSnv = pcawgSnvVcfs, pcawgIndel = pcawgIndelVcfs, stringsAsFactors = F)
rm(hmfVcfs, pcawgSnvVcfs, pcawgIndelVcfs, hmfSvVcfs, pcawgSvVcfs)

load_sv_file <- function(directory, fileName) {
  vcf <- read.vcfR(paste0(directory, fileName), verbose = F)
  vcf  = vcfR2tidy(vcf, single_frame = F)$fix %>% filter(!is.na(MATEID)) 
  
  vcf.id = vcf %>% select(ID) %>% mutate(IDRank = row_number())
  vcf.mateid = vcf %>% select(ID = MATEID) %>% mutate(MATEIDRank = row_number())  
  
  vcf = left_join(vcf, vcf.id, by = "ID") %>% left_join(vcf.mateid, by = "ID") %>% mutate(primary = IDRank < MATEIDRank) %>% select(-IDRank, -MATEIDRank)
  
  vcf.start = hmf.vcf
  vcf.end = hmf.vcf
  vcf.joined = inner_join(vcf.start, vcf.end, by = c("ID" = "MATEID"), suffix = c("_start", '_end')) %>% filter(primary_start) %>% select(-primary_start, -primary_end) %>%
    select(startChromosome = CHROM_start, endChromosome = CHROM_end, startPosition = POS_start, endPosition = POS_end)
}


hmf_data_from_file <- function(hmfFile) {
  hmf.vcf <- read.vcfR(paste0(hmfDirectory, hmfFile), verbose = F)
  hmf.vcf <- vcfR2tidy(hmf.vcf, single_frame = T)$dat
  hmf.vcf$source = "HMF"
  hmf.vcf = hmf.vcf %>%
    mutate(
      type = ifelse(nchar(REF) > 1 | nchar(ALT) > 1, "MNV", "SNV"),
      type = ifelse(nchar(REF)!=nchar(ALT), "INDEL", type))
  
  return (hmf.vcf)
}

hmf_data_from_db <- function(sample) {
  hmf.vcf = hmfSomaticsDB %>% filter(sampleId == sample) %>%
    select(CHROM = chromosome, POS = position, REF = ref, ALT = alt, FILTER = filter, everything()) %>%
    mutate(
      type = ifelse(type == "MNP","MNV", type), 
      type = ifelse(type == "SNP","SNV", type),
      source = "HMF")
  
  return (hmf.vcf)
}

#i = 5
#pcawgSnvFile = vcfFilenames[i, "pcawfSnv"]

combined_data_frame <- function(hmf.vcf, pcawgSnvFile, pcawgIndelFile, hmfDirectory, pcawgDirectory) {
    chromosomeFactors = c(1:22, 'X', 'Y')

    pcawg.vcf.snv <- read.vcfR(paste0(pcawgDirectory, pcawgSnvFile), verbose = F)
    pcawg.vcf.snv <- vcfR2tidy(pcawg.vcf.snv, info_only = T)$fix
    pcawg.vcf.snv = pcawg.vcf.snv %>%
        mutate(CHROM = factor(CHROM, chromosomeFactors, ordered = T)) %>%
        arrange(CHROM, POS) %>%
        mutate(
        prevPOS = lag(POS),
        prevCHROM = lag(CHROM),
        nextPOS = lead(POS),
        nextCHROM = lead(CHROM),
        nextPOS2 = lead(POS, 2),
        nextCHROM2 = lead(CHROM, 2),
        prevDistance = ifelse(CHROM == prevCHROM, prevPOS - POS, 0),
        nextDistance = ifelse(CHROM == nextCHROM, nextPOS - POS, 0),
        nextDistance2 = ifelse(CHROM == nextCHROM2, nextPOS2 - POS, 0),
        type = ifelse(prevDistance != -1 & nextDistance == 1, "MNV", "SNV"),
        type = ifelse(prevDistance != -1 & nextDistance == 1 & nextDistance2 == 2, "MNV_LONG", type),
        type = ifelse(prevDistance == -1, "MNV_BODY", type),
        type = ifelse(is.na(type), "SNV", type),
        CHROM = as.character(CHROM)
        ) %>%
        #select(prevCHROM, prevPOS, CHROM, POS,  nextCHROM, nextPOS, nextCHROM2, nextPOS2, prevDistance, nextDistance, nextDistance2, type)
        select(-nextPOS, -nextCHROM, -prevPOS, -prevCHROM, -nextDistance, -prevDistance, -nextCHROM2, -nextPOS2, -nextDistance2)


    pcawg.vcf.indel <- read.vcfR(paste0(pcawgDirectory, pcawgIndelFile), verbose = F)
    pcawg.vcf.indel <- vcfR2tidy(pcawg.vcf.indel, info_only = T)$fix
    pcawg.vcf.indel$type = "INDEL"

    pcawg.vcf <- bind_rows(pcawg.vcf.indel, pcawg.vcf.snv)
    pcawg.vcf$source = "PCAWG"

    combined = bind_rows(hmf.vcf, pcawg.vcf) %>%
        filter(CHROM %in% chromosomeFactors) %>%
        mutate(
        CHROM = factor(CHROM, chromosomeFactors, ordered = T),
        FILTER = ifelse(is.na(FILTER),"PASS", FILTER)) %>%
        filter(FILTER == 'PASS') %>%
        select(source, everything())

    return (combined)
}

add_scope <- function(combined) {
    combined = combined %>%
        mutate(groupColumn = ifelse(type == "MNV", "MNV", paste0(REF,"|",ALT))) %>%
        group_by(CHROM, POS, groupColumn) %>%
        mutate(scope = ifelse(n() > 1, "Shared", source)) %>%
        select(source, scope, everything())

    return (combined)
}

filter_count <- function(combined) {
    filterCount  = combined %>%
        select(CHROM, POS, groupColumn, FILTER, source, type) %>%
        spread(source, FILTER) %>%
        group_by(type, HMF, PCAWG)%>%   #Group by type also
        count() %>%
        select(sampleId, everything()) %>%
        arrange(-n)
}

######### EXECUTE ALGORITHM
query = "select * from somaticVariant where sampleId not like 'CPCT%' AND filter = 'PASS'"
dbProd = dbConnect(MySQL(), dbname='landscape_paper_validation', groups="RAnalysis")
hmfSomaticsDB = dbGetQuery(dbProd, query)
dbDisconnect(dbProd)
rm(dbProd, query)

hmfDirectory = "/Users/jon/hmf/analysis/pcawg/ours/"
pcawgDirectory = "/Users/jon/hmf/analysis/pcawg/theirs/"

combinedSnvMnvIndel = data.frame()
for (i in c(1:nrow(vcfFilenames))) {
    cat("Processing", i, "of 24", "\n")

    sample = vcfFilenames[i, "hmf"]
    sample = gsub("_post_processed.vcf.gz", "", sample)
    sample = gsub(".*_", "", sample)
    
    hmf.vcf = hmf_data_from_db(sample)
    #hmf.vcf = hmf_data_from_file(vcfFilenames[i, "hmf"])
    
    combined = combined_data_frame(hmf.vcf, vcfFilenames[i, "pcawfSnv"], vcfFilenames[i, "pcawgIndel"], hmfDirectory, pcawgDirectory)
    combinedWithScope = add_scope(combined) %>% mutate(sampleId = sample)

    combinedSnvMnvIndel = bind_rows(combinedWithScope, combinedSnvMnvIndel)
}

save(combinedSnvMnvIndel, file = "/Users/jon/hmf/analysis/pcawg/combinedSnvMnvIndel.RData")


######### STATS!!
load(file = "/Users/jon/hmf/analysis/pcawg/combinedSnvMnvIndel.RData")
combinedSnvMnvIndel = combinedSnvMnvIndel %>% filter(!type %in% c('MNV_BODY', 'MNV_LONG'), type != 'MNV' | nchar(REF) <= 2) 
combinedSnvMnvIndel %>% group_by(source, type) %>% count() %>% spread(type, n)

combinedSnvMnvIndel %>% filter(source == 'PCAWG') %>% group_by(type, scope) %>% count() %>% group_by(type) %>% mutate(percent = n / sum(n))
combinedSnvMnvIndel %>% filter(source == 'PCAWG', t_alt_count > 5) %>% group_by(type, scope) %>% count() %>% group_by(type) %>% mutate(percent = n / sum(n))
combinedSnvMnvIndel %>% filter(source == 'PCAWG', t_alt_count <= 5) %>% group_by(type, scope) %>% count() %>% group_by(type) %>% mutate(percent = n / sum(n))
combinedSnvMnvIndel %>% filter(source == 'HMF') %>% group_by(type, scope) %>% count() %>% group_by(type) %>% mutate(percent = n / sum(n))

######### PLOT!!!
load(file = "/Users/jon/hmf/analysis/pcawg/combinedSnvMnvIndel.RData")
combinedSnvMnvIndelSummary = combinedSnvMnvIndel %>%
  filter(!type %in% c('MNV_BODY', 'MNV_LONG'), type != 'MNV' | nchar(REF) <= 2)  %>%
  group_by(source,sampleId, type) %>% 
  count() %>% 
  spread(source,n, fill = 0)

load(file = "/Users/jon/hmf/analysis/pcawg/combinedSV.RData")
combinedSVSummary = combinedSV %>% 
  mutate(source = ifelse(source == "hmf", "HMF", "PCAWG")) %>%
  filter(startChromosome!=endChromosome|endPosition-startPosition>1000) %>%  
  group_by(source,sampleId) %>% 
  count() %>% 
  spread(source,n, fill = 0) %>%
  mutate(type = 'SV')

combinedSummary = bind_rows(combinedSnvMnvIndelSummary, combinedSVSummary)

lm_eqn = function(m) {
  
  l <- list(a = format(coef(m)[1], digits = 3),
            b = format(abs(coef(m)[2]), digits = 3),
            r2 = format(summary(m)$r.squared, digits = 3));
  
  if (coef(m)[2] >= 0)  {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
  }
  
  as.character(as.expression(eq));                 
}

singleBlue = "#6baed6"
singleRed = "#d94701"

create_plot <- function(data, xLab, yLab) {
  ggplot(data=data,aes(PCAWG,HMF)) + 
    geom_smooth(method = 'lm', colour = singleBlue, se = F, alpha = 0.2) +
    geom_point(color = "black") +
    geom_text(aes(x = xLab, y = yLab, label = lm_eqn(lm(PCAWG ~ HMF, data))), parse = TRUE, size = 3, hjust = 0) +
    xlab("PCAWG Pipeline") + ylab("HMF Pipeline") +# ggtitle("Structural Variants") +
    scale_x_continuous(expand = c(0.02,0.02)) + scale_y_continuous(expand = c(0.02,0.02)) +
    theme(panel.border = element_blank(), axis.ticks = element_blank(), panel.grid.minor = element_blank())
}


p1 = create_plot(combinedSummary %>% filter(type == 'SNV'), 10000, 30000) + xlab("") + ggtitle("SNV")
p2 = create_plot(combinedSummary %>% filter(type == 'MNV'), 25, 100) + xlab("") + ylab("") + ggtitle("MNV")
p3 = create_plot(combinedSummary %>% filter(type == 'INDEL'), 5000, 30000) + ggtitle("INDEL")
p4 = create_plot(combinedSummary %>% filter(type == 'SV'), 100, 400) + ylab("") + ggtitle("SV") + ylim(0, 500)

pTotal = plot_grid(p1,p2,p3,p4, nrow = 2, labels = "AUTO")
pTotal
save_plot("~/hmf/RPlot/Extended Figure 4 - Pipeline Comparison.png", pTotal, base_width = 8, base_height = 8)
