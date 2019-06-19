
load(file = "~/hmf/analysis/cohort/processed/genePanelCv.RData")
newPanel = genePanelCv %>%
  mutate(
    isActionable = actionableVariant | actionableAmplification | actionableDeletion | actionableDrup,
    source = ifelse(isActionable, "actionable", "cosmic"),
    source = ifelse(inDnds, "dnds", source)) %>%
  select(gene, n_ind, n_mis, n_spl, n_non, d_ind, d_mis, d_spl, d_non, classification, source)

load(file = "~/hmf/analysis/cohort/processed/genePanelCvPaper.RData")
oldPanel = genePanelCvPaper %>%
  mutate(
    source = ifelse(hmf | martincorena, "dnds", "cosmic")) %>%
  select(gene = gene_name, wind_cv, wmis_cv, wspl_cv, wnon_cv, classification, source)

driverCatalogDifferences = full_join(newPanel, oldPanel, by = "gene", suffix = c(".new", ".old")) %>%
  mutate(
    d_ind.old =  ifelse(n_ind>0, n_ind * pmax(0,(wind_cv-1)/wind_cv), 0),
    d_mis.old =  ifelse(n_mis>0, n_mis * pmax(0,(wmis_cv-1)/wmis_cv), 0),
    d_spl.old =  ifelse(n_spl>0, n_spl * pmax(0,(wspl_cv-1)/wspl_cv), 0),
    d_non.old =  ifelse(n_non>0, n_non * pmax(0,(wnon_cv-1)/wnon_cv), 0),
    changedClassification = classification.old != classification.new,
    new = is.na(classification.old),
    removed = is.na(classification.new)
    ) %>%
  select(-wind_cv, -wmis_cv, -wspl_cv, -wnon_cv)

save(driverCatalogDifferences, file = "~/hmf/analysis/cohort/processed/driverCatalogDifferences.RData")

View(driverCatalogDifferences %>% mutate(chmis=d_mis-d_mis.old,chind=d_ind-d_ind.old,chspl=d_spl-d_spl.old,chnon=d_non-d_non.old) %>% 
       select(gene,chmis,chind,chspl,chnon,changedClassification,removed,new,everything()))

load(file = "~/hmf/RData/Processed/geneCopyNumberAmplificationTargets.RData")
load(file = "~/hmf/RData/Processed/geneCopyNumberDeleteTargets.RData")
oldAmps = geneCopyNumberAmplificationTargets %>% ungroup() %>% select(gene = target) %>% mutate(old = T)
oldDels = geneCopyNumberDeleteTargets %>% ungroup() %>% select(gene = target) %>% mutate(old = T)

load(file = "~/hmf/analysis/cohort/processed/geneCopyNumberAmplificationTargets.RData")
load(file = "~/hmf/analysis/cohort/Processed/geneCopyNumberDeleteTargets.RData")
newAmps = geneCopyNumberAmplificationTargets %>% ungroup() %>% select(gene = target) %>% mutate(new = T)
newDels = geneCopyNumberDeleteTargets %>% ungroup() %>% select(gene = target) %>% mutate(new = T)

ampDifferences = full_join(newAmps, oldAmps, by = "gene")
delDifferences = full_join(newDels, oldDels, by = "gene")
save(ampDifferences, delDifferences, file = "~/hmf/analysis/cohort/processed/ampDelDifferences.RData")

