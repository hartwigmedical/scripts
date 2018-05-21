dnds_excess <- function(HmfRefCDSCv) {
  HmfRefCDSCv$prob_mis = ifelse(HmfRefCDSCv$n_mis>0,pmax(0,(HmfRefCDSCv$wmis_cv-1)/HmfRefCDSCv$wmis_cv),0)
  HmfRefCDSCv$prob_non = ifelse(HmfRefCDSCv$n_non,pmax(0,(HmfRefCDSCv$wnon_cv-1)/HmfRefCDSCv$wnon_cv),0)
  HmfRefCDSCv$prob_spl = ifelse(HmfRefCDSCv$n_spl>0,pmax(0,(HmfRefCDSCv$wspl_cv-1)/HmfRefCDSCv$wspl_cv),0)
  HmfRefCDSCv$prob_ind = ifelse(HmfRefCDSCv$n_ind>0,pmax(0,(HmfRefCDSCv$wind_cv-1)/HmfRefCDSCv$wind_cv),0)
  HmfRefCDSCv$excess_mis = HmfRefCDSCv$prob_mis*HmfRefCDSCv$n_mis
  HmfRefCDSCv$excess_non = HmfRefCDSCv$prob_non*HmfRefCDSCv$n_non
  HmfRefCDSCv$excess_spl = HmfRefCDSCv$prob_spl*HmfRefCDSCv$n_spl
  HmfRefCDSCv$excess_ind = HmfRefCDSCv$prob_ind*HmfRefCDSCv$n_ind

  return (HmfRefCDSCv)
}


dnds_annotate_somatics <- function(annotmuts, somatics) {
  result = annotmuts %>% select(sampleId = sampleID, chromosome = chr, position = pos, ref = ref, alt = mut, gene, impact) %>%  left_join(somatics, by = c("sampleId","chromosome","position","ref", "alt"))

  result = result %>%
    mutate(impact = ifelse(impact == "Essential_Splice", "Splice", impact)) %>%
    mutate(impact = ifelse(impact == "no-SNV", paste0(substr(canonicalCodingEffect, 1, 1), tolower(substring(canonicalCodingEffect, 2))), impact)) %>%
    mutate(impact = ifelse(impact == "None", paste0(substr(worstCodingEffect, 1, 1), tolower(substring(worstCodingEffect, 2))), impact)) %>%
    mutate(impact = ifelse(impact %in% c("Stop_loss","Nonsense_or_frameshift"), "Nonsense", impact)) %>%
    mutate(impact = ifelse(type == "INDEL" & impact == "Missense", "Inframe", impact)) %>%
    mutate(impact = ifelse(type == "INDEL" & impact == "Nonsense", "Frameshift", impact)) %>%
    filter(!is.na(impact), impact != "None")

  return (result)
}
