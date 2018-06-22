library(purple)
detach("package:purple", unload=TRUE)
library(purple)
library(RMySQL)
library(StructuralVariantAnnotation)
#dbmanta = dbConnect(MySQL(), dbname = "hmfpatients_pilot")
dbcon = dbConnect(MySQL(), dbname = "gridss_test")

grcohort = query_structural_variants_as_GRanges(dbcon, cohort=data.frame(sampleId=c("CPCT02100013T", "CPCT02100013TII")))
grbp = grcohort[!is.na(grcohort$partner)]
grbp$remote=paste0(seqnames(partner(grbp)), ":", start(partner(grbp)), "-", end(partner(grbp)))
grbp$svlen = abs(start(partner(grbp)) - start(grbp))
grbp$type=simpleEventType(grbp)

grbp1 = grbp[grbp$sampleId=="CPCT02100013T"]
grbp2 = grbp[grbp$sampleId=="CPCT02100013TII"]

# calculate overlaps

hits = findBreakpointOverlaps(grbp1, grbp2) # require exact/overlapping call positions
grbp1$hit = NA_character_
grbp2$hit = NA_character_
grbp1$hit[hits$queryHits] = names(grbp2)[hits$subjectHits]
grbp2$hit[hits$subjectHits] = names(grbp1)[hits$queryHits]

ggplot(as.data.frame(c(grbp1, grbp2))) +
  aes(x=ploidy, fill=type) +
  geom_histogram() +
  facet_wrap(ifelse(is.na(hit), "No match", "Match") ~ sampleId)

as.data.frame(c(grbp1, grbp2)) %>% filter(is.na(hit)) %>% arrange(seqnames, start) %>% select(28, 1:27) %>% View()


# Is the issue the filtering? How many of these variants are found in the full GRIDSS?
full_vcf = readVcf("CPCT02100013R_CPCT02100013TII.gridss.vcf", "hg19")
bpgr = breakpointRanges(full_vcf, unpartneredBreakends=FALSE) # allow 300bp margin

hits = findBreakpointOverlaps(grbp1, bpgr, maxgap=300) %>%
  mutate(QUAL=bpgr$QUAL[subjectHits]) %>%
  arrange(QUAL)
grbp1$hits_all = NA_character_
grbp1$hits_all[hits$queryHits] = names(bpgr)[hits$subjectHits]
grbp1$hit_status=ifelse(is.na(grbp1$hits_all), "No match", ifelse(is.na(grbp1$hit), "filtered", "match" ))
table(grbp1$hit_status)

as.data.frame(grbp2) %>% filter(is.na(hit)) %>% arrange(seqnames, start) %>% select(28, 1:27) %>% View()
info(full_vcf[grbp2$vcfid[is.na(grbp2$hit)]]) %>% as.data.frame() %>% View()

ggplot(info(full_vcf[grbp2$vcfid]) %>%
    as.data.frame() %>%
    mutate(hasMatch=!is.na(grbp2$hit)) %>%
    mutate(FILTER=rowRanges(full_vcf[grbp2$vcfid])$FILTER)) +
  aes(x=RP+SR, y=ASSR+ASRP, color=hasMatch, shape=FILTER) +
  geom_point() +
  geom_jitter() +
  scale_x_log10() +
  scale_y_log10()


########
# Breakend comparison
grbe = grcohort[is.na(grcohort$partner)]
grbe1 = grbe[grbe$sampleId=="CPCT02100013T"]
grbe2 = grbe[grbe$sampleId=="CPCT02100013TII"]
hits = findOverlaps(grbe1, grbe2, maxgap=8) %>% as.data.frame()
grbe1$hit = "Private"
grbe2$hit = "Private"
grbe1$hit[hits$queryHits] = "Shared breakend"
grbe2$hit[hits$subjectHits] = "Shared breakend"

grbe1$hit[queryHits(findOverlaps(grbe1, grbp2, maxgap=8))] = "Breakend to breakpoint"
grbe2$hit[queryHits(findOverlaps(grbe2, grbp1, maxgap=8))] = "Breakend to breakpoint"

ggplot(as.data.frame(c(grbe1, grbe2))) +
  aes(x=ploidy, y=QUAL, color=hit) +
  geom_point() +
  facet_wrap( ~ sampleId) +
  labs(title="Shared vs private single breakends")

ggplot(info(full_vcf[grbe2$vcfid]) %>%
         as.data.frame() %>%
         mutate(hit=grbe2$hit) %>%
         mutate(FILTER=rowRanges(full_vcf[grbe2$vcfid])$FILTER)) +
  aes(x=BSC+BUM, y=BASSR+BASRP, color=hit, shape=FILTER) +
  geom_point() +
  geom_jitter() +
  scale_x_log10() +
  scale_y_log10()


