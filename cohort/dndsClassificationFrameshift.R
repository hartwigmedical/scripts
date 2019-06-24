##### BRINGING IN FRAMESHIFT ETC EXPERIMENT


load(file = "~/hmf/analysis/cohort/processed/genePanelCv.RData")
oldGenePanel = genePanelCv %>% select(gene, oldClassification = classification) %>% distinct()


load(file = "~/hmf/analysis/cohort/reference/hpcExonicSomatics.RData")
load("~/hmf/analysis/cohort/processed/dndsFilteredAnnotatedMutations.RData")
combinedSomatics = dnds_annotate_somatics(dndsFilteredAnnotatedMutations, hpcExonicSomatics, distinguishIndels = T)
somaticCountsPerGene = combinedSomatics %>% group_by(gene.x, impact) %>% count() %>% filter(impact != "") %>% ungroup() %>%
  mutate(
  impact = ifelse(impact == "Synonymous", "n_syn", impact),
  impact = ifelse(impact == "Missense", "n_mis", impact),
  impact = ifelse(impact == "Splice", "n_spl", impact),
  impact = ifelse(impact == "Nonsense", "n_non", impact),
  impact = ifelse(impact == "Frameshift", "n_frameshift", impact),
  impact = ifelse(impact == "Inframe", "n_inframe", impact)) %>%
  spread(impact, n, fill = 0) %>%
  select(gene = gene.x, everything())

load("~/hmf/analysis/cohort/processed/completeGenePanel.RData")
dndsGenePanel = completeGenePanel %>% filter(martincorenaDnds | hmfDnds | cosmicCurated | actionableVariant | actionableDrup) %>% distinct()

load("~/hmf/analysis/cohort/processed/HmfRefCDSCv.RData")
genePanelCv = HmfRefCDSCv %>% filter(cancerType == "All", gene_name %in% dndsGenePanel$gene) %>%
  select(-n_ind, -n_mis, -n_spl, -n_non) %>%
  left_join(somaticCountsPerGene %>% select(gene, n_mis, n_spl, n_non, n_frameshift, n_inframe), by = c("gene_name" = "gene")) %>%
  mutate(
  d_ind =  ifelse(n_frameshift>0, n_frameshift * pmax(0,(wind_cv-1)/wind_cv), 0),
  d_mis =  ifelse(n_mis>0, n_mis * pmax(0,(wmis_cv-1)/wmis_cv), 0),
  d_spl =  ifelse(n_spl>0, n_spl * pmax(0,(wspl_cv-1)/wspl_cv), 0),
  d_non =  ifelse(n_non>0, n_non * pmax(0,(wnon_cv-1)/wnon_cv), 0)
  ) %>%
  select(gene = gene_name, wmis_cv, wnon_cv, wspl_cv, wnon_cv, d_ind, d_mis, d_spl, d_non, n_frameshift, n_inframe, n_mis, n_spl, n_non)
genePanelCv = left_join(genePanelCv, dndsGenePanel[, c("gene", "cosmicTsg", "cosmicOncogene", "hmfDnds", "martincorenaDnds", "actionableVariant", "actionableDrup", "actionableDrupCategory", "actionableAmplification", "actionableDeletion", "hmfAmplification", "hmfDeletion")], by = "gene")

trainTsg = genePanelCv %>% filter(cosmicTsg, !cosmicOncogene, gene != "PHOX2B") %>% mutate(tsg = 1) %>% select(wmis_cv, wnon_cv, n_frameshift, n_inframe, tsg)
trainOnco = genePanelCv %>% filter(cosmicOncogene, !cosmicTsg, gene != "PHOX2B") %>% mutate(tsg = 0) %>% select(wmis_cv, wnon_cv, n_frameshift, n_inframe,  tsg)
trainData = rbind(trainTsg, trainOnco)
model <- glm(tsg ~.,family=binomial(link='logit'),data=trainData)
summary(model)
anova(model, test="Chisq")
genePanelCv$continuousClassification <- predict(model,newdata = genePanelCv %>% select(wmis_cv, wnon_cv, n_frameshift, n_inframe),type='response')
genePanelCv[is.na(genePanelCv)] <- F
genePanelCv$inDnds <- genePanelCv$hmfDnds | genePanelCv$martincorenaDnds
genePanelCv$unambigiousCosmic = xor(genePanelCv$cosmicTsg, genePanelCv$cosmicOncogene)

genePanelCv$classification = ifelse(genePanelCv$continuousClassification > 0.5,"tsg","onco")
genePanelCv$classification = ifelse(!genePanelCv$inDnds & !genePanelCv$cosmicTsg & genePanelCv$cosmicOncogene, "onco", genePanelCv$classification)
genePanelCv$classification = ifelse(!genePanelCv$inDnds & genePanelCv$cosmicTsg & !genePanelCv$cosmicOncogene, "tsg", genePanelCv$classification)
genePanelCv$classification = ifelse(!genePanelCv$inDnds & !genePanelCv$unambigiousCosmic & (genePanelCv$hmfAmplification | genePanelCv$actionableAmplification), "onco", genePanelCv$classification)
genePanelCv$classification = ifelse(!genePanelCv$inDnds & !genePanelCv$unambigiousCosmic & (genePanelCv$hmfDeletion | genePanelCv$actionableDeletion), "tsg", genePanelCv$classification)

genePanelCvWithFrameshift = genePanelCv %>% left_join(oldGenePanel, by = "gene")
save(genePanelCvWithFrameshift, file = "~/hmf/analysis/cohort/processed/genePanelCvWithFrameshift.RData")
