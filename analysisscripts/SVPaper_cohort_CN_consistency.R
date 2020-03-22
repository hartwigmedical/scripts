source("libSVPaper.R")
# TODO: export CN from LINX
# ggplot(sva_svs) +
#   aes(x=(CNChgEnd-Ploidy)/max(CNChgEnd-Ploidy)) +
#   geom_histogram() +
#   facet_wrap(~ResolvedType, scales="free")
# ggplot(sva_svs) +
#   aes(
#     y=(CNChgStart-Ploidy)/max(CNChgStart-Ploidy),
#     x=(CNChgEnd-Ploidy)/max(CNChgEnd-Ploidy)) +
#   facet_wrap(~ ResolvedType)
#   geom_point()


library(RMySQL)
dbConn = dbConnect(MySQL(), dbname = "hmfpatients", port=3307)
# copyNumberdbTable = dbGetQuery(dbConn, "SELECT * FROM copyNumber")
#
#
# bedf = sva_svs %>%
#   dplyr::select(SampleId, Id, Type, Chr=ChrStart, Pos=PosStart, Orient=OrientStart, AF=AFStart, CN=CNStart, CNChg=CNChgStart, Ploidy, PloidyMin, PloidyMax, ResolvedType) %>%
#   bind_rows(sva_svs %>% filter(ChrEnd != 0) %>%
#             dplyr::select(SampleId, Id, Type, Chr=ChrEnd, Pos=PosEnd, Orient=OrientEnd, AF=AFEnd, CN=CNEnd, CNChg=CNChgEnd, Ploidy, PloidyMin, PloidyMax, ResolvedType)) %>%
#   mutate(
#     cnLeftEnd=Pos + ifelse(Orient==-1, -1, 0),
#     cnRightStart=Pos + ifelse(Orient==-1, 0, 1)
#   )
# cnNames = names(copyNumberdbTable)
# names(copyNumberdbTable) = paste0("left.", cnNames)
# becndf = left_join(bedf, copyNumberdbTable, by=c("SampleId"="left.sampleId", "Chr"="left.chromosome","cnLeftEnd"="left.end"))
# names(copyNumberdbTable) = paste0("right.", cnNames)
# becndf = left_join(becndf, copyNumberdbTable, by=c("SampleId"="right.sampleId", "Chr"="right.chromosome","cnRightStart"="right.start"))
# names(copyNumberdbTable) = cnNames
#
# becndf = becndf %>%
#   filter(
#     !is.na(becndf$left.id) &
#     !is.na(becndf$right.id) &
#     Type != "NONE" &
#     left.copyNumberMethod %in% c("BAF_WEIGHTED", "STRUCTURAL_VARIANT") &
#     right.copyNumberMethod %in% c("BAF_WEIGHTED", "STRUCTURAL_VARIANT")) %>%
#   mutate(
#     # TODO: adjust to +ve means SV ploidy is higher than expected
#     cn_error_abs=abs(left.copyNumber - right.copyNumber) - Ploidy,
#     cn_error_rel=cn_error_abs / pmax(left.copyNumber, right.copyNumber))
# # save.image("cohort_CN_consistency.RData")
#
# ggplot(becndf) +
#   aes(x=cn_error_abs, fill=factor(Orient)) +
#   geom_histogram(bins=200) +
#   scale_x_continuous(limits=c(-1, 1)) +
#   facet_wrap(right.copyNumberMethod ~ left.copyNumberMethod, scales="free")
#
# # TODO: add expected theoretical distribution based on flanking CN and read depth
# ggplot(becndf) +
#   aes(x=cn_error_rel, fill=left.copyNumberMethod) +
#   geom_histogram(bins=200) +
#   scale_x_continuous(limits=c(-0.25, 0.25), labels=percent_format()) +
#   facet_wrap(right.copyNumberMethod ~ left.copyNumberMethod)
#
#
#
cndf = dbGetQuery(dbConn, "SELECT * FROM copyNumber")
svdf = dbGetQuery(dbConn, "SELECT * FROM structuralVariant")

svcndf = dbGetQuery(dbConn, "SELECT ploidy, adjustedCopyNumberStart, adjustedCopyNumberEnd, adjustedCopyNumberChangeStart, adjustedCopyNumberChangeEnd FROM structuralVariant")


ggplot(svcndf) +
  aes(x=(adjustedCopyNumberChangeEnd - adjustedCopyNumberChangeStart)) +
  scale_x_continuous(limits=c(-0.5, 0.5)) +
  geom_histogram(bins=100) +
  labs(title="PURPLE copy number consistency across breakpoints", x="copy number error", y="count")

ggplot(rbind(
    svcndf %>% filter(!is.na(adjustedCopyNumberChangeEnd)) %>% mutate(delta=adjustedCopyNumberChangeEnd-ploidy) %>% dplyr::select(delta),
    svcndf %>% mutate(delta=adjustedCopyNumberChangeStart-ploidy) %>% dplyr::select(delta)
  ) %>% filter(delta != 0)) +
  aes(x=delta) +
  scale_x_continuous(limits=c(-5, 5)) +
  geom_histogram(bins=100) +
  labs(title="PURPLE copy number consistency with\nGRIDSS breakpoint copy number inferred from VAF\n/supporting fragment count.", x="copy number error", y="count")



