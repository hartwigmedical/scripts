pcawg_structural_variant_summary <- function(dbConnect, cohort) {
  sampleIdString = paste("'", cohort$sampleId, "'", collapse = ",", sep = "")
  query = paste(
    "SELECT sampleId, count(*) as SV",
    " FROM structuralVariant ",
    "WHERE sampleId in (",sampleIdString, ") AND (type <> 'BND' OR ploidy >= 0.2) AND (startChromosome <> endChromosome OR endPosition-startPosition>1000)",
    " GROUP BY sampleId",
    sep = "")
  result = dbGetQuery(dbConnect, query)
  return (result)
}

load(file = "~/hmf/RData/reference/highestPurityCohort.RData")
dbProd = dbConnect(MySQL(), dbname='hmfpatients_20180418', groups="RAnalysis")
pcawgSV = pcawg_structural_variant_summary(dbProd, highestPurityCohort)
dbDisconnect(dbProd)
rm(dbProd)


load(file = "~/hmf/RData/reference/allSomatics_p1.RData")
allSomatics_p1 = allSomatics_p1 %>% filter(type != 'MNP' | nchar(ref) == 2)
somatics_summary_p1 = cohort_somatic_summary(allSomatics_p1)
rm(allSomatics_p1)

load(file = "~/hmf/RData/reference/allSomatics_p2.RData")
allSomatics_p2 = allSomatics_p2 %>% filter(type != 'MNP' | nchar(ref) == 2)
somatics_summary_p2 = cohort_somatic_summary(allSomatics_p2)
rm(allSomatics_p2)

pcawgSomaticsSummary = rbind(somatics_summary_p1, somatics_summary_p2) %>% filter(sampleId %in% highestPurityCohort$sampleId) %>%
  mutate(
    INDEL = TOTAL_INDEL,
    MNV = TOTAL_MNV,
    SNV = TOTAL_SNV) %>%
  select(sampleId, INDEL, MNV, SNV)

filteredHmfMutationalLoad = highestPurityCohort %>% select(sampleId, cancerType) %>%
  left_join(pcawgSomaticsSummary, by = "sampleId") %>%
  left_join(pcawgSV, by = "sampleId")

save(filteredHmfMutationalLoad, file = "~/hmf/RData/processed/filteredHmfMutationalLoad.RData")


########### PLOT FUNCTIONS
display_cancer_types <- function(cancerTypes) {
  
  for (i in 1:length(cancerTypes)) {
    n = paste0('(n=',hpcCancerTypeCounts[hpcCancerTypeCounts$cancerType == cancerTypes[i], ] %>% pull(N), ')')
    
    if (cancerTypes[i] == "Mesothelioma") {
      cancerTypes[i]= "Meso- thelioma"
    }
    if (cancerTypes[i] == "Esophagus") {
      cancerTypes[i]= "Esoph- agus"
    }
    if (cancerTypes[i] == "Colon/Rectum") {
      cancerTypes[i]= "Colon/ Rectum"
    }
    if (cancerTypes[i] == "Head and neck") {
      cancerTypes[i]= "Head & Neck"
    }
    if (cancerTypes[i] == "Bone/Soft tissue") {
      cancerTypes[i]= "Bone/Soft Tissue"
    }
    if (cancerTypes[i] == "Urinary tract") {
      cancerTypes[i]= "Urinary Tract"
    }
    
    cancerTypes[i] = paste0(cancerTypes[i], " ", n)
  }
  
  return (label_wrap_gen(10)(cancerTypes))
}

somaticColours = c("#a6611a","#dfc27d","#80cdc1","#018571")
somaticColours = setNames(somaticColours, c("Adjusted HMF SNV","PCAWG SNV", "PCAWG MNV", "Adjusted HMF MNV"))
somaticLinetypes = c("solid","longdash","longdash","solid")
somaticLinetypes = setNames(somaticLinetypes, c("Adjusted HMF SNV","PCAWG SNV", "PCAWG MNV", "Adjusted HMF MNV"))

indelSVColours = c("#d01c8b","#f1b6da","#b8e186","#4dac26")
indelSVColours = setNames(indelSVColours, c("Adjusted HMF INDEL","PCAWG INDEL", "PCAWG SV", "Adjusted HMF SV"))
indelSVLinetypes = c("solid","longdash","longdash","solid")
indelSVLinetypes = setNames(indelSVLinetypes, c("Adjusted HMF INDEL","PCAWG INDEL", "PCAWG SV", "Adjusted HMF SV"))

########### PLOT
load(file = "~/hmf/RData/processed/highestPurityCohortSummary.RData")
oldMutationalLoad = highestPurityCohortSummary %>% 
  select(sampleId, cancerType, ends_with("INDEL"), ends_with("SNV"), ends_with("MNV"), TRL, DEL, INS, INV, DUP) %>%
  mutate(cancerType = factor(cancerType, levels = cancerTypeFactors)) %>%
  mutate(
    INDEL = TOTAL_INDEL,
    MNV = TOTAL_MNV,
    SNV = TOTAL_SNV,
    SV = TRL + DEL + INS + INV + DUP) %>% 
  select(sampleId, cancerType, INDEL, SNV, MNV, SV)

load(file = "~/hmf/RData/processed/filteredHmfMutationalLoad.RData")
load(file = '~/hmf/RData/Reference/hpcCancerTypeCounts.RData')
cancerTypeFactors =  factor(hpcCancerTypeCounts$cancerType, levels = hpcCancerTypeCounts$cancerType)

filteredHmfMutationalLoad = filteredHmfMutationalLoad %>%
  mutate(
    SV = (1 - 0.37) * SV,
    SNV = (1 - 0.07) * SNV,
    MNV = (1 - 0.21) * MNV,
    INDEL = (1 - 0.45) * INDEL)


hmfMutationalLoad = filteredHmfMutationalLoad
#hmfMutationalLoad = oldMutationalLoad

pcawgRaw = read.csv("~/hmf/resources/PCAWG_counts.txt", sep = '\t', stringsAsFactors = F)
pcawg_histology_tier2 = sort(unique(pcawgRaw$histology_tier2))

pcawgCancerTypeMapping = data.frame(histology_tier2 = pcawg_histology_tier2, cancerType = pcawg_histology_tier2, stringsAsFactors = F)
pcawgCancerTypeMapping[pcawgCancerTypeMapping$histology_tier2 == "Bladder", "cancerType"] = "Urinary tract"
pcawgCancerTypeMapping[pcawgCancerTypeMapping$histology_tier2 == "Bone/SoftTissue", "cancerType"] = "Bone/Soft tissue"
pcawgCancerTypeMapping[pcawgCancerTypeMapping$histology_tier2 == "Cervix", "cancerType"] = NA
pcawgCancerTypeMapping[pcawgCancerTypeMapping$histology_tier2 == "Head/Neck", "cancerType"] = "Head and neck"
pcawgCancerTypeMapping[pcawgCancerTypeMapping$histology_tier2 == "Myeloid", "cancerType"] = "Other"
pcawgCancerTypeMapping[pcawgCancerTypeMapping$histology_tier2 == "Lymphoid", "cancerType"] = "Other"
pcawgCancerTypeMapping[pcawgCancerTypeMapping$histology_tier2 == "Thyroid", "cancerType"] = "Other"
pcawgCancerTypeMapping = pcawgCancerTypeMapping[!is.na(pcawgCancerTypeMapping$cancerType), ]

pcawgMutationalLoad = pcawgRaw %>% left_join(pcawgCancerTypeMapping, by = "histology_tier2") %>%
  filter(!is.na(cancerType)) %>% select(PCAWG_SNV = all.SNVs, PCAWG_INDEL = all.Indels, PCAWG_SV = SV.events, age,PCAWG_MNV = all.MNVs,  cancerType) %>%
  mutate(source = "PCAWG")

combinedMutationalLoad =  hmfMutationalLoad %>% select(sampleId, cancerType, INDEL, SNV, MNV, SV) %>%
  mutate(source = "HMF") %>%
  bind_rows(pcawgMutationalLoad) %>%
  mutate(cancerType = factor(cancerType, levels = cancerTypeFactors)) %>%
  filter(cancerType != "Other")

combinedMutationalLoad = combinedMutationalLoad %>% 
  group_by(cancerType) %>% 
  mutate(
    medianSNV = median(SNV, na.rm = T), 
    medianMNV = median(MNV, na.rm = T), 
    medianPCAWG_SNV = median(PCAWG_SNV, na.rm = T),
    medianPCAWG_MNV = median(PCAWG_MNV, na.rm = T),
    medianINDEL = median(INDEL, na.rm = T), 
    medianSV = median(SV, na.rm = T), 
    medianPCAWG_INDEL = median(PCAWG_INDEL, na.rm = T),
    medianPCAWG_SV = median(PCAWG_SV, na.rm = T)
  ) %>% ungroup()



theme_set(theme_bw() +  theme(axis.text=element_text(size=5), axis.title=element_text(size=7),legend.text = element_text(size=5), strip.text.x = element_text(size = 5, face = "plain"),
                              legend.background=element_blank(), legend.key=element_blank(), panel.spacing = unit(1, "pt"), plot.title = element_text(size=7)))

p3 = ggplot(data=combinedMutationalLoad) +
  stat_ecdf(size = 0.3, aes(SNV,color='Adjusted HMF SNV',linetype='Adjusted HMF SNV'), geom = "step", pad = FALSE) + geom_segment(size = 0.3,aes(x = medianSNV, xend = medianSNV, y = 0.25, yend = 0.75, color='Adjusted HMF SNV'), show.legend = F) + 
  stat_ecdf(size = 0.3,aes(MNV,color='Adjusted HMF MNV',linetype='Adjusted HMF MNV') ,geom = "step", pad = FALSE) + geom_segment(size = 0.3,aes(x = medianMNV, xend = medianMNV, y = 0.25, yend = 0.75, color='Adjusted HMF MNV'), show.legend = F) + 
  stat_ecdf(size = 0.4,aes(PCAWG_SNV,color='PCAWG SNV',linetype='PCAWG SNV'), geom = "step", pad = FALSE) + geom_segment(size = 0.3,aes(x = medianPCAWG_SNV, xend = medianPCAWG_SNV, y = 0.25, yend = 0.75, color='PCAWG SNV'), show.legend = F) + 
  stat_ecdf(size = 0.4,aes(PCAWG_MNV,color='PCAWG MNV',linetype='PCAWG MNV'), geom = "step", pad = FALSE) + geom_segment(size = 0.3,aes(x = medianPCAWG_MNV, xend = medianPCAWG_MNV, y = 0.25, yend = 0.75, color='PCAWG MNV'), show.legend = F) + 
  scale_x_log10() + 
  facet_grid(~cancerType, scales = "free_x", labeller = labeller(cancerType = display_cancer_types)) + theme(panel.spacing = unit(1, "pt")) +  
  scale_colour_manual(name = "Combined", values=somaticColours) + 
  scale_linetype_manual(name = "Combined", values = somaticLinetypes) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(),
        panel.grid.minor.x = element_blank(),
         legend.position="top", legend.title = element_blank(), 
        legend.background=element_blank(), legend.key=element_blank()) +
  
  xlab("Somatic Variants") +
  guides(colour = guide_legend(nrow = 1)) +
  coord_flip()
p3

p4 = ggplot(data=combinedMutationalLoad) +
  stat_ecdf(size = 0.3,aes(INDEL, color='Adjusted HMF INDEL', linetype = 'Adjusted HMF INDEL'),geom = "step", pad = FALSE) + geom_segment(size = 0.3,aes(x = medianINDEL, xend = medianINDEL, y = 0.25, yend = 0.75, color='Adjusted HMF INDEL'), show.legend = F) + 
  stat_ecdf(size = 0.3,aes(SV,color='Adjusted HMF SV',linetype='Adjusted HMF SV'),geom = "step", pad = FALSE) + geom_segment(size = 0.3,aes(x = medianSV, xend = medianSV, y = 0.25, yend = 0.75, color='Adjusted HMF SV'), show.legend = F) +
  stat_ecdf(size = 0.4,aes(PCAWG_INDEL,color='PCAWG INDEL',linetype='PCAWG INDEL'),geom = "step", pad = FALSE) + geom_segment(size = 0.3,aes(x = medianPCAWG_INDEL, xend = medianPCAWG_INDEL, y = 0.25, yend = 0.75, color='PCAWG INDEL'), show.legend = F) +
  stat_ecdf(size = 0.4,aes(PCAWG_SV,color='PCAWG SV', linetype='PCAWG SV'),geom = "step",  pad = FALSE) + geom_segment(size = 0.3,aes(x = medianPCAWG_SV, xend = medianPCAWG_SV, y = 0.25, yend = 0.75, color='PCAWG SV'), show.legend = F) +
  scale_x_log10() +  facet_grid(~cancerType) + theme(panel.spacing = unit(1, "pt")) +
  scale_colour_manual(name = "Combined", values=indelSVColours) + 
  scale_linetype_manual(name = "Combined", values = indelSVLinetypes) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background = element_blank(), 
        strip.text.x = element_blank(), 
        legend.position="bottom",
        legend.title = element_blank(), 
        legend.background=element_blank(), legend.key=element_blank()
        ) +
  xlab("Somatic Variants") +
  xlab("INDELs & SVs") + 
  guides(colour = guide_legend(nrow = 1)) +
  coord_flip()

p4
pFigure1Revisited = plot_grid(p3, p4, ncol=1, align="v", rel_heights = c(1.2, 1), labels = c("e", "f"), label_size = 8)
pFigure1Revisited
dev.off()
save_plot("~/hmf/RPlot/Extended Figure 5 - Adjusted TMB Comparison.png", pFigure1Revisited, base_width = 14, base_height = 5)
#save_plot("~/hmf/RPlot/Extended Figure 10 - OLD.png", pFigure1Revisited, base_width = 14, base_height = 5)


########### WILSON COX
cancerTypes = unique(pcawgMutationalLoad$cancerType)
result = data.frame(cancerType = cancerTypes, stringsAsFactors = F)
variantTypes = c("SNV","INDEL","SV","MNV")
selectedVariant = "SV"

hmfMutationalLoad = filteredHmfMutationalLoad
#hmfMutationalLoad = oldMutationalLoad
for (selectedCancerType in cancerTypes) {
  for (selectedVariant in variantTypes) {
    pcawgSelected = pcawgMutationalLoad %>% filter(cancerType == selectedCancerType) %>% pull(paste0("PCAWG_", selectedVariant))
    pcawgSelected = pcawgSelected[!is.na(pcawgSelected)]
    
    hmfSelected = hmfMutationalLoad %>% filter(cancerType == selectedCancerType) %>% pull(selectedVariant)
    hmfSelected = hmfSelected[!is.na(hmfSelected)]
    
    w = wilcox.test(pcawgSelected,hmfSelected,conf.int = T)
    result[result$cancerType == selectedCancerType, selectedVariant] <- w[["p.value"]]
    result[result$cancerType == selectedCancerType, paste('medianHMF',selectedVariant,sep='_')] <- median(hmfSelected)
    result[result$cancerType == selectedCancerType, paste('medianPCAWG',selectedVariant,sep='_')] <- median(pcawgSelected)
    result[result$cancerType == selectedCancerType, paste('diff',selectedVariant,sep='_')] <- w[["estimate"]]
  }
}

adjustedHmfWilcoxResult = result
save(adjustedHmfWilcoxResult, file = "~/hmf/RData/Processed/adjustedHmfWilcoxResult.RData")

