library("libgridss.R")
sample1 = "C:/hartwig/down/CPCT02100013R_CPCT02100013T.gridss.vcf"
sample2 = "C:/hartwig/down/CPCT02100013R_CPCT02100013TII.gridss.vcf"
vcf1 = readVcf(sample1, "hg19")
vcf2 = readVcf(sample2, "hg19")

#########
# Breakpoint comparison
bpgr1 = breakpointRanges(vcf1)
bpgr2 = breakpointRanges(vcf2)

bpgr1$type = simpleEventType(bpgr1)
bpgr2$type = simpleEventType(bpgr2)

bpgr1$af = gridss_somatic_bp_af(bpgr1, vcf1)
bpgr2$af = gridss_somatic_bp_af(bpgr2, vcf2)

bpgr1$FILTER = gridss_breakpoint_filter(bpgr1, vcf1, pon_dir=pon_dir)
bpgr2$FILTER = gridss_breakpoint_filter(bpgr2, vcf2, pon_dir=pon_dir)

bpgr1$FILTER = .addFilter(bpgr1$FILTER, "LOW_Qual", bpgr1$QUAL < gridss.min_qual)
bpgr2$FILTER = .addFilter(bpgr2$FILTER, "LOW_Qual", bpgr2$QUAL < gridss.min_qual)

hits = findBreakpointOverlaps(bpgr1, bpgr2, maxgap=8) %>%
  mutate(filters1 = str_count(bpgr1$FILTER[queryHits], stringr::fixed(";")),
         filters2 = str_count(bpgr2$FILTER[queryHits], stringr::fixed(";")))
bestHits1 = hits %>% group_by(queryHits) %>% top_n(1, filters1)
bpgr1$matchFilter = "NoMatch"
bpgr1$matchFilter[bestHits1$queryHits] = bpgr2$FILTER[bestHits1$subjectHits]
bestHits2 = hits %>% group_by(subjectHits) %>% top_n(1, filters2)
bpgr2$matchFilter = "NoMatch"
bpgr2$matchFilter[bestHits2$subjectHits] = bpgr1$FILTER[bestHits2$queryHits]

# How many calls are missed by one
ggplot(bpgr1[bpgr1$FILTER == "" & bpgr1$matchFilter == "NoMatch"] %>% as.data.frame()) +
  aes(x=af, fill=type) +
  geom_histogram()
ggplot(bpgr2[bpgr2$FILTER == "" & bpgr2$matchFilter == "NoMatch"] %>% as.data.frame()) +
  aes(x=af, fill=type) +
  geom_histogram()

# What are the filters that are applied in one and not the other?
table(bpgr1$matchFilter[bpgr1$FILTER == "" & bpgr1$matchFilter != "" & !is.na(bpgr1$matchFilter)])
table(bpgr2$matchFilter[bpgr2$FILTER == "" & bpgr2$matchFilter != "" & !is.na(bpgr2$matchFilter)])

# Should Be Unfiltered bpgr1
sbuf1 = bpgr1[bpgr1$FILTER != "" & bpgr1$matchFilter == ""]
table(sbuf1$FILTER)


ggplot(info(vcf1[bpgr1$FILTER == "" & bpgr1$matchFilter == "NoMatch"]) %>% as.data.frame()) +
  aes(x=SR+RP, y=ASQ) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10()

# TODO: contamination-based somatic model
# Probability that we observed this many variant reads assuming X% tumour contamination and Y AF
# vs
# Probability that we observed this many variant reads assuming het/hom variant





