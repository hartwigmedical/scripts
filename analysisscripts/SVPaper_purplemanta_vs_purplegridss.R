source("libSVPaper.R")

dbconn1 = dbConnect(MySQL(), dbname='hmfpatients_20180418', groups="RAnalysis")
dbconn2 = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
data1 = list(
  sv = DBI::dbGetQuery(dbconn1, "
      SELECT id,
        sampleId,
        startChromosome,
        endChromosome,
        startPosition,
        endPosition,
        startOrientation,
        endOrientation,
        ploidy,
        adjustedStartAF AS adjustedAFStart,
        adjustedEndAF AS adjustedAFEnd
      FROM structuralVariant
      WHERE filter = 'PASS' OR filter = 'INFERRED'"),
  cn=DBI::dbGetQuery(dbconn1, "
      SELECT id,
        sampleId,
        chromosome,
        start,
        end,
        copyNumber
      FROM copyNumber"))
data2 = list(
  sv = DBI::dbGetQuery(dbconn2, "
      SELECT id,
        sampleId,
        startChromosome,
        endChromosome,
        startPosition,
        endPosition,
        startOrientation,
        endOrientation,
        ploidy,
        adjustedAFStart,
        adjustedAFEnd
      FROM structuralVariant
      WHERE filter = 'PASS' OR filter = 'INFERRED'"),
  cn=DBI::dbGetQuery(dbconn2, "
      SELECT id,
        sampleId,
        chromosome,
        start,
        end,
        copyNumber
      FROM copyNumber"))
save(data1, data2, file = "D:/hartwig/svtoolkit/purplemanta_vs_purplegridss.RData")

load("D:/hartwig/svtoolkit/purplemanta_vs_purplegridss.RData")
data1$sv = data1$sv %>% filter(sampleId %in% data2$sv$sampleId)
data2$sv = data2$sv %>% filter(sampleId %in% data1$sv$sampleId)
data1$cn = data1$cn %>% filter(sampleId %in% data2$cn$sampleId)
data2$cn = data2$cn %>% filter(sampleId %in% data1$cn$sampleId)


annotate_data = function(hmftables, endadjustment = 0) {
  svdf = hmftables$sv %>%
    mutate(isBreakend=endChromosome==0)
  cndf = hmftables$cn
  bpdf = bind_rows(
    svdf %>%
      mutate(posid = paste0(sampleId, "_", startChromosome, "_", startPosition, "_", startOrientation)) %>%
      dplyr::select(svid=id, posid, svploidy=ploidy, svaf=adjustedAFStart, isBreakend),
    svdf %>%
      mutate(posid = paste0(sampleId, "_", endChromosome, "_", endPosition, "_", endOrientation)) %>%
      dplyr::select(svid=id, posid, svploidy=ploidy, svaf=adjustedAFEnd, isBreakend)) %>%
    mutate(posid = ifelse(duplicated(posid), paste0("dup_", posid, "_", svid), posid)) %>%
    # filter out foldback inversions
    filter(!duplicated(posid))
  rownames(bpdf) = bpdf$posid
  cndf$startsv = paste0(cndf$sampleId, "_", cndf$chromosome, "_", cndf$start, "_", "-1")
  cndf$endsv = paste0(cndf$sampleId, "_", cndf$chromosome, "_", cndf$end + endadjustment, "_", "1")
  cndf$startsv[!(cndf$startsv %in% bpdf$posid)] = NA
  cndf$endsv[!(cndf$endsv %in% bpdf$posid)] = NA

  cndf = cndf %>%
    group_by(sampleId, chromosome) %>%
    arrange(start) %>%
    mutate(
      hasPrev=!is.na(lag(start)),
      hasNext=!is.na(lead(start)),
      prevTransitionSvs=as.numeric(!is.na(lag(endsv))) + as.numeric(!is.na(startsv)),
      nextTransitionSvs=as.numeric(!is.na(lead(startsv))) + as.numeric(!is.na(endsv)))
  cntranstion = bind_rows(
    cndf %>%
      mutate(
        cnadj=lead(copyNumber)) %>%
      filter(
        hasNext,
        !is.na(endsv),
        nextTransitionSvs == 1) %>%
      ungroup() %>%
      dplyr::select(cnid=id, cn=copyNumber, cnadj, posid=endsv),
    cndf %>%
      mutate(
        cnadj=lag(copyNumber)) %>%
      filter(
        hasPrev,
        !is.na(startsv),
        prevTransitionSvs == 1) %>%
      ungroup() %>%
      dplyr::select(cnid=id, cn=copyNumber, cnadj, posid=startsv)) %>%
    inner_join(bpdf, by="posid") %>%
    mutate(cndelta=cn - cnadj)
  pairedsvcn = cntranstion %>%
    group_by(svid) %>%
    filter(n() == 2)
  pairedsvcn = pairedsvcn %>% inner_join(pairedsvcn, by=c("svid"), suffix=c("", ".2")) %>%
    filter(posid != posid.2) %>%
    group_by(svid) %>%
    top_n(1)
  cndf = cndf %>% left_join()
  return(list(cn=cndf, bp=pairedsvcn))
}

mantadf = annotate_data(data1, endadjustment=1)
gridssdf = annotate_data(data2, endadjustment=0)

purplecndf = bind_rows(
    gridssdf$cn %>% mutate(caller="purple/gridss2"),
    mantadf$cn %>% mutate(caller="purple/manta")) %>%
  ungroup()
purplebpdf = bind_rows(
  gridssdf$bp %>% mutate(caller="purple/gridss2"),
    mantadf$bp %>% mutate(caller="purple/manta")) %>%
  ungroup() %>%
  mutate(cndelta_percentage_inconsistency=abs(cndelta-cndelta.2)/pmax(cn, cn.2, cnadj, cnadj.2))

colour_scheme=c("#1b9e77", "#7570b3")

ggplot(purplebpdf) +
  aes(x=cndelta_percentage_inconsistency, fill=caller) +
  geom_histogram(aes(y=..density..), position="identity", bins=50, alpha=0.5, boundary=0) +
  #geom_density(cut=0) +
  scale_x_continuous(limits=c(0, 0.3), expand=c(0,0), labels = scales::percent) +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_manual(values=colour_scheme) +
  labs(title="Distribution of copy number\ninconsistency across breakpoints", x="Copy number inconsistency (%)", y="") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
figsave("purple_manta_vs_gridss_cn_inconsistency", height=4, width=5)

ggplot(purplebpdf %>% ungroup() %>% sample_n(100000)) +
  aes(x=svploidy, y=abs(cndelta), colour=caller) +
  geom_point(size=0.2, alpha=0.03) +
  #geom_density2d() +
  scale_x_continuous(limits=c(0, 4)) +
  scale_y_continuous(limits=c(0, 4)) +
  scale_colour_manual(values=colour_scheme) +
  facet_grid(~ caller) +
  labs(title="Consistency of structural variant and flanking segment copy number", x="Copy number delta across breakend", y="Structural variant copy number")
figsave("purple_manta_vs_gridss_sv_segment_scatter", height=5, width=6)

sv_cn_diffdf = purplebpdf %>% mutate(
    delta=abs(svploidy-cndelta),
    delta_portion=abs(svploidy-cndelta)/pmax(abs(svploidy), abs(cndelta)))
ggplot(sv_cn_diffdf) +
  aes(x=abs(delta), fill=caller) +
  geom_histogram(position="identity", alpha=0.5, bins=100) +
  scale_x_continuous(limits=c(0, 1)) +
  scale_fill_manual(values=c("#1b9e77", "#7570b3")) +
  labs(title="Difference of SV copy number and copy number segment delta", x="copy number difference")
ggplot(sv_cn_diffdf) +
  aes(x=abs(delta_portion), fill=caller) +
  geom_histogram(position="identity", alpha=0.6, bins=100) +
  scale_x_continuous(limits=c(0, 1), labels = scales::percent) +
  scale_fill_manual(values=c("#1b9e77", "#7570b3")) +
  labs(title="Relative difference of SV copy number\nand copy number segment delta", x="Relative difference")
figsave("purple_manta_vs_gridss_sv_segment_relative_error", height=4, width=5)
purplebpdf %>%
  group_by(caller) %>%
  filter(!is.na(cndelta_percentage_inconsistency) & !is.infinite(cndelta_percentage_inconsistency)) %>%
  summarise(median(cndelta_percentage_inconsistency), mean(cndelta_percentage_inconsistency))

sv_cn_diffdf %>% group_by(caller) %>%
  summarise(
    delta_mean=mean(abs(delta), na.rm=TRUE),
    delta_median=median(abs(delta), na.rm=TRUE),
    delta_sd=sd(delta, na.rm=TRUE),
    mean_delta_portion=mean(abs(delta_portion), na.rm=TRUE),
    median_delta_portion=median(abs(delta_portion), na.rm=TRUE))
t.test(
  sv_cn_diffdf %>% filter(caller == "purple/gridss2") %>% pull(delta_portion),
  sv_cn_diffdf %>% filter(caller != "purple/gridss2") %>% pull(delta_portion))

bind_rows(
  purplecndf %>% filter(hasNext) %>% dplyr::select(caller, transitionsvs=nextTransitionSvs),
  purplecndf %>% filter(hasPrev) %>% dplyr::select(caller, transitionsvs=prevTransitionSvs)) %>%
  group_by(caller) %>%
  summarise(unexplained_cn_transition=1 - sum(transitionsvs > 0) / n())
# fisher's exact test










