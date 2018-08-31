#install.packages("vcfR")
# scp -C jon@hmf-datastore:/data/experiments/180815_pcawg_on_HMF_pipeline_validation/PCAWG_HMF/*/*_post_processed.vcf.gz* /Users/jon/hmf/analysis/pcawg/ours/
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
rm(hmfVcfs, pcawgSnvVcfs, pcawgIndelVcfs)


combined_data_frame <- function(hmfFile, pcawgSnvFile, pcawgIndelFile, hmfDirectory, pcawgDirectory) {
    chromosomeFactors = c(1:22, 'X', 'Y')

    hmf.vcf <- read.vcfR(paste0(hmfDirectory, hmfFile), verbose = F)
    hmf.vcf <- vcfR2tidy(hmf.vcf, single_frame = T)$dat
    hmf.vcf$source = "HMF"
    hmf.vcf = hmf.vcf %>%
    mutate(
    type = ifelse(nchar(REF) > 1 | nchar(ALT) > 1, "MNV", "SNV"),
    type = ifelse(nchar(REF)!=nchar(ALT), "INDEL", type))

    pcawg.vcf.snv <- read.vcfR(paste0(pcawgDirectory, pcawgSnvFile), verbose = F)
    pcawg.vcf.snv <- vcfR2tidy(pcawg.vcf.snv, info_only = T)$fix
    pcawg.vcf.snv = pcawg.vcf.snv %>%
        mutate(CHROM = factor(CHROM, chromosomeFactors, ordered = T)) %>%
        arrange(CHROM, POS) %>%
        mutate(
        nextPOS = lead(POS),
        nextCHROM = lead(CHROM),
        prevPOS = lag(POS),
        prevCHROM = lag(CHROM),
        nextDistance = ifelse(CHROM == nextCHROM, nextPOS - POS, 0),
        prevDistance = ifelse(CHROM == prevCHROM, prevPOS - POS, 0),
        type = ifelse(nextDistance == 1 & prevDistance != -1, "MNV", "SNV"),
        type = ifelse(prevDistance == -1, "MNV_BODY", type),
        CHROM = as.character(CHROM)
        ) %>%
        select(-nextPOS, -nextCHROM, -prevPOS, -prevCHROM, -nextDistance, -prevDistance)


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
    sample = first(combined %>% filter(!is.na(Indiv)) %>% pull(Indiv))

    filterCount  = combined %>%
        select(CHROM, POS, groupColumn, FILTER, source, type) %>%
        spread(source, FILTER) %>%
        group_by(type, HMF, PCAWG)%>%   #Group by type also
        count() %>%
        mutate(sampleId = sample) %>%
        select(sampleId, everything()) %>%
        arrange(-n)
}

######### EXECUTE ALGORITHM
hmfDirectory = "/Users/jon/hmf/analysis/pcawg/ours/"
pcawgDirectory = "/Users/jon/hmf/analysis/pcawg/theirs/"

result = data.frame()
allSomatics = data.frame()
for (i in c(1:nrow(vcfFilenames))) {
    cat("Processing", i, "of 24")

    combined = combined_data_frame(vcfFilenames[i, "hmf"], vcfFilenames[i, "pcawfSnv"], vcfFilenames[i, "pcawgIndel"], hmfDirectory, pcawgDirectory)
    sample = first(combined %>% filter(!is.na(Indiv)) %>% pull(Indiv))
    combinedWithScope = add_scope(combined) %>% mutate(sampleId = sample)

    allSomatics = bind_rows(combinedWithScope, allSomatics)

    #filterCount = filter_count(combinedWithScope %>% filter(type != "MNV_BODY"))
    #result = bind_rows(filterCount, result)
}

save(allSomatics, file ="~/hmf/allSomatics.RData")
save(result, file ="~/hmf/result.RData")


result$PCAWGSummarised = ifelse(grepl('NORMAL',result$PCAWG),'PON',ifelse(grepl('PASS',result$PCAWG)|is.na(result$PCAWG),result$PCAWG,'OTHER'))
result$HMFSummarised = ifelse(grepl('PON',result$HMF),'PON',ifelse(grepl('PASS',result$HMF)|is.na(result$HMF),result$HMF,'OTHER'))
summary=(result %>% filter(HMFSummarised=='PASS'|PCAWGSummarised=='PASS') %>% unite(HMF_PCAWG,HMFSummarised,PCAWGSummarised) %>% group_by(HMF_PCAWG,type,sampleId) %>% summarise(sum=sum(n)) %>% spread(HMF_PCAWG,sum) )
library(ggplot2)
ggplot()+geom_point(data=summary %>% filter(isIndel==F),aes(PASS_PASS,PASS_OTHER),color = "red")+
    geom_point(data=summary %>% filter(isIndel==F),aes(PASS_PASS,PASS_NA),color = "blue") +
    geom_point(data=summary %>% filter(isIndel==F),aes(PASS_PASS,PASS_PON),color = "green")+
    geom_point(data=summary %>% filter(isIndel==F),aes(PASS_PASS,NA_PASS),color = "black")# +scale_x_log10()+scale_y_log10()

