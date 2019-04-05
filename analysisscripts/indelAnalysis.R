library(dplyr)
library(ggplot2)


### SUMMARISE INDELS
indel_summary <- function(somatics) {
  
  lengthLabels = c("1", "2", "4", "8", "16", "32+")
  lengthBreaks = c(1, 2, 4, 8, 16, 32, 100000000)
  
  repeatCountLabels = c("0", "1", "2", "3", "4+")
  repeatCountBreaks = c(0, 1, 2, 3, 4, 10000000)
  
  result = somatics %>% filter(filter == 'PASS', type == 'INDEL') %>%
    mutate(
      isMicrohomology = microhomology != "",
      indelType = ifelse(nchar(ref) > nchar(alt), "DEL", "INS"),
      length = abs(nchar(ref) - nchar(alt)), lengthFactor = cut(length, breaks = lengthBreaks, labels = lengthLabels, include.lowest = T, right = F),
      repeatCountFactor = cut(repeatCount, breaks = repeatCountBreaks, labels = repeatCountLabels, include.lowest = T, right = F)
    ) %>%
    select(sampleId, chromosome, position, ref, alt, indelType, isMicrohomology, length, lengthFactor, repeatCount, repeatCountFactor) %>%
    group_by(sampleId, indelType, isMicrohomology,lengthFactor, repeatCountFactor)  %>% 
    count()
}


#load(file = "~/hmf/RData/reference/allSomatics_p1.RData")
#allIndels_p1 = indel_summary(allSomatics_p1)
#rm(allSomatics_p1)
#load(file = "~/hmf/RData/reference/allSomatics_p2.RData")
#allIndels_p2 = indel_summary(allSomatics_p2)
#rm(allSomatics_p2)
#allIndels = bind_rows(allIndels_p1, allIndels_p2)
#save(allIndels, file = "~/hmf/analysis/sv/RData/allIndels.RData")


#svs = read.csv(file = "~/hmf/analysis/sv/input/SVA_VIS_SVS.csv", stringsAsFactors = F)
#svLengthLabels = c("NA", "<1000", "1k-10k", "10k-100k", "100k+")
#svLengthBreaks = c(0, 1, 1000, 10000, 100000, 10000000000)
#allSvs = svs %>%
#  mutate(svLength = ifelse(ChrStart == ChrEnd, abs(PosStart - PosEnd), 0), svLengthFactor = cut(svLength, breaks = svLengthBreaks, labels = svLengthLabels, include.lowest = T, right = F)) %>%
#  group_by(SampleId, Type, ResolvedType, svLengthFactor) %>% count() %>%
#  select(sampleId = SampleId, type = Type, resolvedType = ResolvedType, svLengthFactor, n)
#save(allSvs, file = "~/hmf/analysis/sv/RData/allSvs.RData")


plot_indel_v_sv <- function(title, indel, sv) {
  df = inner_join(indel, sv, by = "sampleId") %>% left_join(allPurity %>% select(sampleId, cancerType), by = "sampleId")
  
  ggplot(data=df, aes(x = svCount, y = indelCount)) +
    geom_segment(aes(x = 1e1, xend=1e4, y = 12400, yend=12380), linetype = "dashed") + annotate("text", x = 2000, y = 20000, label = "MSI Threshold", size = 3, hjust = 0) +
    geom_point(aes(color = cancerType)) + 
    scale_color_manual(values = cancerTypeColours) + 
    scale_x_continuous(trans="log10", limits = c(1e1, 1e4)) + 
    scale_y_continuous(trans="log10", limits = c(1e2, 1e6)) + 
    theme(legend.position = "none", legend.title = element_blank()) + guides(colour = guide_legend(ncol = 1)) +
    xlab("SVs") + ylab("Indels") + ggtitle(title)
}

load(file = "~/hmf/analysis/sv/RData/allSvs.RData")
load(file = "~/hmf/analysis/sv/RData/allIndels.RData")
load(file = "~/hmf/RData/Reference/cancerTypeColours.RData")
load(file = "~/hmf/RData/Reference/allPurity.RData")

## DISPLAY AVAILABLE LEVELS
levels(allSvs$svLengthFactor)
levels(allIndels$lengthFactor)
levels(allIndels$repeatCountFactor)

dfSV = allSvs %>% 
  filter(type == 'DUP', resolvedType == 'SimpleSV', svLengthFactor %in% c('<1000', '1k-10k')) %>%   # MODIFY AS YOU SEE FIT
  group_by(sampleId) %>% summarise(svCount = sum(n))

dfIndel = allIndels %>%
  filter(lengthFactor %in% c('1','2'), repeatCountFactor == "4+") %>%   # MODIFY AS YOU SEE FIT
  group_by(sampleId) %>% summarise(indelCount = sum(n))

plot_indel_v_sv("Anything you want here", dfIndel, dfSV)


dfSV = allSvs %>% 
  filter(type == 'DEL', resolvedType == 'SimpleSV', svLengthFactor %in% c('<1000')) %>%   # MODIFY AS YOU SEE FIT
  group_by(sampleId) %>% summarise(svCount = sum(n))

dfIndel = allIndels  %>%
  filter(lengthFactor=='8+', repeatCountFactor != "4+",isMicrohomology==T,indelType=='DEL') %>%   # MODIFY AS YOU SEE FIT
  group_by(sampleId) %>% summarise(indelCount = sum(n))

plot_indel_v_sv("Anything you want here", dfIndel, dfSV)




nonRepeatIndels = allIndels %>% filter(!isRepeat) %>% group_by(sampleId) %>% summarise(indelCount = sum(n))
simpleSVs = svs %>% group_by(SampleId, ResolvedType) %>% summarise(svCount = n()) %>% filter(ResolvedType == 'SimpleSV') %>% select(sampleId = SampleId, svCount)
plot_indel_v_sv("Non-repeat Indels vs Simple SVs", nonRepeatIndels, simpleSVs, allPurity)

nonRepeatInserts = allIndels %>% filter(!isRepeat, indelType == 'INS') %>% group_by(sampleId) %>% summarise(indelCount = sum(n))
simpleSVInserts = svs %>% filter(ResolvedType == 'SimpleSV', Type == 'INS') %>% group_by(SampleId, ResolvedType) %>% summarise(svCount = n())  %>% select(sampleId = SampleId, svCount)
plot_indel_v_sv("Non-repeat Inserts vs Simple SV Inserts", nonRepeatInserts, simpleSVInserts, allPurity)






