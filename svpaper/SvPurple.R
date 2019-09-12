library(tidyr)
library(dplyr)
library(ggplot2)
detach("package:purple", unload=TRUE)
library(purple)

singleBlue = "#6baed6"

hmfTheme = theme_bw() +
  theme(
    axis.title = element_text(size=7),
    axis.text = element_text(size=5),
    axis.ticks = element_blank(),
    legend.title = element_text(size=5), 
    legend.text = element_text(size=5),
    legend.key.size = unit(0.2, "cm"),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.text = element_text(size=6)
  )
theme_set(hmfTheme)  

load(file = "~/hmf/analysis/cohort/cohort.RData")
load(file = "~/hmf/analysis/cohort/hpcStructuralVariants.RData")
#load(file = "~/hmf/analysis/cohort/hpcInferredStructuralVariants.RData")

hpcStructuralVariants = hpcStructuralVariants %>% filter(sampleId %in% highestPurityCohort$sampleId, !grepl("PON", filter))
hpcStructuralVariantsSummary = hpcStructuralVariants %>% mutate(single = is.na(endChromosome) | endChromosome == 0) %>% group_by(filter, single, recovered) %>% count() %>%
  mutate(weight = ifelse(single, 1, 2), n = n * weight) %>%
  group_by(filter, recovered) %>% summarise(n = sum(n)) %>%
  ungroup() %>%
  mutate(relPercentage = n / sum(n))

svInferredDF = hpcStructuralVariants %>% filter(filter =='INFERRED') %>% group_by(sampleId) %>% summarise(inferred = n())
svPassingDF = hpcStructuralVariants %>% 
  filter(filter  == 'PASS') %>% mutate(type = ifelse(recovered, "recovered", "supported")) %>%
  group_by(sampleId, type) %>% 
  count() %>%
  spread(type, n, fill = 0)

svDF = full_join(svInferredDF, svPassingDF, by = 'sampleId')
svDF[is.na(svDF)] <- 0  
svDF = svDF %>% gather(type, value, inferred, recovered, supported) %>% group_by(sampleId) %>% 
  mutate(sampleTotal = sum(value), relValue = value / sampleTotal) 

svOrdering = svDF %>% group_by(sampleId) %>% summarise(n = sum(value)) %>% arrange(-n)
svDF = svDF %>% ungroup() %>% mutate(sampleId = factor(sampleId, levels = svOrdering$sampleId, ordered = T))

p1 = ggplot(svDF) +
  geom_bar(aes(x = sampleId, y = value, fill = type), stat = "identity", width = 1) + 
  theme(axis.text.x = element_blank(), legend.title = element_blank(), plot.background = element_blank(), legend.text = element_text(size=7)) + xlab("Samples") + ylab("SV Counts") + 
  coord_cartesian(ylim = c(0, 2000))

#ggplot(svDF) +
#  geom_bar(aes(x = sampleId, y = relValue, fill = type), stat = "identity", width = 1) + 
#  theme(axis.text.x = element_blank()) + xlab("Samples") + ylab("Relative SV Counts") +
#  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0%", "25%", "50%", "75%", "100%"))


sv2DF = svDF %>% group_by(type) %>% summarise(value = sum(value)) %>% mutate(relValue = value / sum(value))
p2 = ggplot(sv2DF) +
  geom_bar(aes(x = 1, y = relValue, fill = type), stat = "identity", width = 1) + 
  theme(axis.text.x = element_blank()) + xlab("") + ylab("Total Proportion") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0%", "25%", "50%", "75%", "100%")) + 
  theme(legend.position = "none", panel.grid = element_blank()) +
  coord_cartesian(ylim = c(0, 1))
  
p3 = plot_grid(p1, p2, nrow = 1, rel_widths = c(8, 1))
p3
ggplot2::ggsave("~/hmf/analysis/svPaper/plot/PurpleSV.pdf", p3, width = 200, height = 100, units = "mm", dpi = 300)
ggplot2::ggsave("~/hmf/analysis/svPaper/plot/PurpleSV.png", p3, width = 200, height = 100, units = "mm", dpi = 300)



load(file = '~/hmf/analysis/cohort/hpcCopyNumbers.RData')
names(hpcCopyNumbers)
centromeres = hpcCopyNumbers %>% ungroup() %>%
  filter(segmentStartSupport == "CENTROMERE" | segmentEndSupport == "CENTROMERE") %>% 
  filter(!chromosome %in% c("13","14","15", "21", "22")) %>%
  mutate(chromosome = factor(chromosome, levels = c(1:22, "X", "Y"), ordered = T)) %>%
  group_by(sampleId, chromosome) %>% 
  mutate(type = ifelse(segmentStartSupport == "CENTROMERE", "start", "end")) %>%
  select(sampleId, chromosome, type, copyNumber) %>%
  spread(type, copyNumber) %>%
  mutate(changeAcrossCentromere = abs(end - start)) 


ggplot(centromeres) +
  stat_ecdf(aes(x = changeAcrossCentromere),  geom = "point", pad = F, size = 0.05) +
  coord_flip() +
  facet_grid(~chromosome) + theme(axis.text.x = element_blank()) + ylab("Proportion of samples") + xlab("Change Across Centromere") 
  scale_x_log10(breaks=c(2,11,101,1001), labels=c(1,10,100,1000))
  


  
  #Base level accuracy in copy number predictions
  #has direct clinical relevance in accurately calling disrupted tumor suppressor genes as we found in the same cohort that over
  #Z%
  #of homozygous deletions of TSG drivers [CITE] are caused by deletion of 10k bases or less 
  #[FIGURE 3B].
  
  load(file = "~/hmf/analysis/cohort/reference/hpcGeneCopyNumberDeletes.RData")
  load(file = "~/hmf/analysis/cohort/reference/hpcDriverCatalog.RData")


hpcDriverDels = hpcDriverCatalog %>% filter(driver == 'DEL', category == 'TSG') %>% select(gene, sampleId) %>% left_join(hpcGeneCopyNumberDeletes, by = c("gene", "sampleId")) %>%
  mutate(delLength = minRegionEnd - minRegionStart + 1, bucket = cut(delLength, breaks = c(0, 1000, 10000, 100000, 1e100), labels = c("1-1k", "1k-10k", "10k-100k", "100k+"))) %>%
  group_by(bucket) %>% count() %>% ungroup() %>% mutate(percentage = n / sum(n))

ggplot(hpcDriverDels, aes(bucket, y = percentage)) +
  geom_bar(stat = "identity") + xlab("Del Length") + ylab("Del Drivers") +
  scale_y_continuous(breaks = 0.1*c(0:7), labels = paste0(10*c(0:7), "%"))
  
