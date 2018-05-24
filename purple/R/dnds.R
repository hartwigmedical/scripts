



onco_mutations <- function(mutations) {
  mutations$nearHotspot <- nearHotspot(mutations)

  result = mutations %>%
    #result = mutations[mutations$gene %in% c("APC", "NRAS","CDKN2A", "KRAS"), ] %>%
    filter(impact %in% c("Inframe", "Missense", "Frameshift")) %>%
    mutate(impact = factor(impact, levels = c("Inframe", "Missense", "Frameshift"))) %>%
    group_by(sampleId, gene) %>%
    mutate(
      row_number = row_number(),
      n = n(),
      driverType = oncoGeneStatus(hotspot, nearHotspot, n),
      redundant = row_number != oncoGeneStatusPrimaryRowNumber(row_number, hotspot, nearHotspot, biallelic, impact)) %>%
    ungroup() %>%
    mutate(
      hotspot = hotspot > 0,
      biallelic = biallelic > 0,
      impact = as.character(impact),
      driverType = as.character(driverType),
      driverType = ifelse(redundant, "Redundant", driverType)) %>%
    filter()

  result = result %>%
    mutate(knownDriver = hotspot | nearHotspot | impact == "Inframe") %>%
    mutate(knownDriver = ifelse(redundant, F, knownDriver))

  return (result)
}


tsg_mutations <- function(mutations) {
  result = mutations %>%
    #filter(gene %in% c("MET", "NRAS","CDKN2A", "KRAS")) %>%
    filter(impact != "Synonymous") %>%
    group_by(sampleId, gene) %>%
    mutate(
      n = n(),
      impact = factor(impact, levels = c("MNV", "Frameshift", "Nonsense", "Splice", "Missense", "Inframe")),
      driverType = factor(tsGeneStatus(hotspot, biallelic, n),levels = c("Biallelic", "MultiHit", "SingleHit")),
      geneStatusPrimaryPositions = tsGeneStatusPrimaryPositions(position, driverType, hotspot, impact),
      redundant = tsGeneStatusRedundant(position, geneStatusPrimaryPositions)) %>%
    select(-geneStatusPrimaryPositions) %>%
    ungroup() %>%
    mutate(
      hotspot = hotspot > 0,
      biallelic = biallelic > 0,
      impact = as.character(impact),
      driverType = as.character(driverType))

  result$driverType <- ifelse(result$redundant, "Redundant", result$driverType)

  # Add knownDriver field: Hotspot, Biallelic, or Multihit (excluding 2xMissense & 2xIndel)
  result = result %>%
    mutate(impact = ifelse(impact %in% c("Inframe", "Frameshift"), "Indel", impact)) %>%
    group_by(sampleId, gene, impact) %>%
    mutate(knownDriver = ifelse(n() > 1 & impact %in% c("Missense", "Indel"), T, F))%>%
    ungroup() %>%
    mutate(knownDriver = ifelse(biallelic | hotspot | driverType == "MultiHit", T, knownDriver)) %>%
    mutate(knownDriver = ifelse(redundant, F, knownDriver))

  return (result)
}


dnds_driver_likelihood <- function(mutations, expectedDriversPerGene) {
  result = mutations %>%
    filter(redundant == F) %>%
    group_by(gene, impact, knownDriver) %>%
    summarise(n = n()) %>%
    mutate(knownDriver = ifelse(knownDriver, "knownDrivers", "unknownDrivers")) %>%
    spread(knownDriver, n, fill = 0) %>%
    left_join(expectedDriversPerGene, by = c("gene","impact")) %>%
    mutate(driverLikelihood = pmin(1, pmax(0, (expectedDrivers - knownDrivers)/(unknownDrivers))))

  result[is.na(result)] <- 0

  return (result)
}


dnds_expected_drivers <- function(HmfRefCDSCv, annotmuts, somatics) {

  geneImpactCount = dnds_annotate_somatics(annotmuts, somatics, distinguishIndels = F) %>%
    group_by(gene, impact) %>%
    summarise(n = n()) %>%
    filter(!impact %in% c("", "Synonymous")) %>%
    mutate(impact = tolower(paste("n", impact, sep = "_"))) %>%
    spread(impact, n, fill = 0)

  result = left_join(HmfRefCDSCv %>% filter(cancerType == 'All'), geneImpactCount, by = c("gene_name" = "gene"))

  result$prob_mis = ifelse(result$n_mis>0,pmax(0,(result$wmis_cv-1)/result$wmis_cv),0)
  result$prob_non = ifelse(result$n_non,pmax(0,(result$wnon_cv-1)/result$wnon_cv),0)
  result$prob_spl = ifelse(result$n_spl>0,pmax(0,(result$wspl_cv-1)/result$wspl_cv),0)
  result$prob_ind = ifelse(result$n_ind>0,pmax(0,(result$wind_cv-1)/result$wind_cv),0)

  result$Missense = result$prob_mis*result$n_missense
  result$Nonsense = result$prob_non*result$n_nonsense
  result$Splice = result$prob_spl*result$n_splice
  result$Indel = result$prob_ind*result$n_indel

  result = result %>%
    select(gene = gene_name, Missense, Nonsense, Splice, Indel) %>%
    gather("impact", "expectedDrivers", Missense, Nonsense, Splice, Indel)

  return (result)
}


#dnds_excess <- function(HmfRefCDSCv) {
#
#
#  HmfRefCDSCv$excess_mis = HmfRefCDSCv$prob_mis*HmfRefCDSCv$n_mis
#  HmfRefCDSCv$excess_non = HmfRefCDSCv$prob_non*HmfRefCDSCv$n_non
#  HmfRefCDSCv$excess_spl = HmfRefCDSCv$prob_spl*HmfRefCDSCv$n_spl
#  HmfRefCDSCv$excess_ind = HmfRefCDSCv$prob_ind*HmfRefCDSCv$n_ind
#
#  return (HmfRefCDSCv)
#}


dnds_annotate_somatics <- function(annotmuts, somatics, distinguishIndels = T) {
  result = annotmuts %>% select(sampleId = sampleID, chromosome = chr, position = pos, ref = ref, alt = mut, gene, impact) %>%  left_join(somatics, by = c("sampleId","chromosome","position","ref", "alt"))

  result = result %>%
    mutate(impact = ifelse(impact == "Essential_Splice", "Splice", impact)) %>%
    mutate(impact = ifelse(impact == "no-SNV", paste0(substr(canonicalCodingEffect, 1, 1), tolower(substring(canonicalCodingEffect, 2))), impact)) %>%
    mutate(impact = ifelse(impact == "None", paste0(substr(worstCodingEffect, 1, 1), tolower(substring(worstCodingEffect, 2))), impact)) %>%
    mutate(impact = ifelse(impact %in% c("Stop_loss","Nonsense_or_frameshift"), "Nonsense", impact)) %>%
    mutate(impact = ifelse(type == "INDEL" & impact == "Missense", "Inframe", impact)) %>%
    mutate(impact = ifelse(type == "INDEL" & impact == "Nonsense", "Frameshift", impact)) %>%
    filter(!is.na(impact), impact != "None")

  if (!distinguishIndels) {
    result[result$impact %in% c("Inframe", "Frameshift"), "impact"] <- "Indel"
  }

  return (result)
}
