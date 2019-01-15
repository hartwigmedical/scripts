collapsePHGVS <- function(x) {
  x = x[!is.na(x) & x != ""]
  if (length(x) == 0) {
    return ("")
  }
  return (paste(x, collapse = ","))
}

hotspot_category <- function(hotspot, nearhotspot) {
  result = ifelse(hotspot, "OnHotspot", "OffHotspot")
  result = ifelse(nearHotspot, "NearHotspot", result)
  return (result)
}

variant_not_polymorphism <- function(type) {
  type = ifelse(type == 'SNP', 'SNV', type)
  type = ifelse(type == 'MNP', 'MNV', type)
  return (type)
}

genomic_coordinates <- function(chromosome, position, ref, alt) {
  return (paste0(chromosome, ":", position, ref, ">", alt))
}


nearHotspot <-function(mutations, distance = 5) {
  hotspots = mutations %>% filter(hotspot > 0) %>% select(chromosome, position) %>% distinct
  hrange <- GRanges(hotspots$chromosome, IRanges(hotspots$position, hotspots$position + distance))
  mrange <- GRanges(mutations$chromosome, IRanges(mutations$position, mutations$position + nchar(mutations$ref) - 1 + distance))

  ol = as.matrix(findOverlaps(hrange, mrange, type="any", select="all"))
  mutations$nearHotspot <- FALSE
  mutations[ol[,2], c("nearHotspot")] <- TRUE
  return (mutations$nearHotspot & !mutations$hotspot)
}

tsGeneStatus <-function(hotspot, biallelic, n) {
  levels = c("Biallelic", "MultiHit", "SingleHit")
  result = ifelse(n > 1, "MultiHit" , "SingleHit")
  result = ifelse(biallelic, "Biallelic", result)
  return (factor(result, levels))
}

tsGeneStatusPrimaryPositions <- function(pos, geneStatus, hotspot, impact) {
  df = data.frame(pos = pos, geneStatus = geneStatus, hotspot = hotspot, impact = impact) %>%
    arrange(geneStatus, -hotspot, impact)

  if (df[1, "geneStatus"] == "MultiHit") {
    return (paste(df[1:2, c("pos")], collapse =","))
  }

  return (as.character(df[1, c("pos")]))
}

tsGeneStatusRedundant <- function(position, driverPositions) {
  result = vector(mode = "character", length = length(position))
  for (i in 1:length(position)) {
    drivers = as.integer(strsplit(driverPositions[i], ",")[[1]])
    if (position[i] %in% drivers) {
      result[i] <- FALSE
    } else {
      result[i] <- TRUE
    }
  }
  return (result)
}

oncoGeneStatus <-function(hotspot, nearHotspot, n) {
  result = ifelse(nearHotspot, "NearHotspot", "Hit")
  result = ifelse(hotspot, "Hotspot", result)
  return (result)
}

oncoGeneStatusPrimaryRowNumber <- function(pos, hotspot, nearHotspot, biallelic, impact) {
  df = data.frame(pos = pos, hotspot = hotspot, nearHotspot = nearHotspot, biallelic = biallelic, impact = impact)
  ordered = df %>% arrange(-hotspot, -nearHotspot, -biallelic, impact)
  if (nrow(ordered) > 0) {
    return (ordered[1, c("pos")])
  }

  return (NA)
}


onco_mutations <- function(mutations) {
  mutations$nearHotspot <- nearHotspot(mutations)

  result = mutations %>%
    #result = mutations[mutations$gene %in% c("APC", "NRAS","CDKN2A", "KRAS"), ] %>%
    filter(impact %in% c("Inframe", "Missense")) %>%
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
    mutate(knownDriver = hotspot | nearHotspot | (impact == "Inframe" & repeatCount < 8)) %>%
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
    mutate(knownDriver = ifelse(biallelic | hotspot, T, F)) %>%
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

dnds_tsg_drivers <- function(sampleSomatics, mutations, expectedDriversPerGene) {
  result = list()
  cohortSize = nrow(sampleSomatics)
  totalSomatics = sampleSomatics %>% ungroup() %>% summarise(total_SNV = sum(sample_SNV), total_INDEL = sum(sample_INDEL))

  tsgMutations = tsg_mutations(mutations)
  result[["tsgMutations"]] <- tsgMutations

  tsgDriverRates = dnds_driver_likelihood(tsgMutations, expectedDriversPerGene)
  result[["tsgDriverRates"]] <- tsgDriverRates

  tsgDrivers = tsgMutations %>%
    filter(redundant == F) %>%
    left_join(tsgDriverRates %>% select(gene, impact, driverLikelihood), by = c("gene","impact")) %>%
    mutate(driverLikelihood = ifelse(knownDriver, 1, driverLikelihood))

  tsgUnknownDriversTotals = tsgDrivers %>%
    group_by(gene, impact) %>%
    filter(!knownDriver) %>%
    summarise(gene_drivers = sum(driverLikelihood), gene_non_drivers = sum(1 - driverLikelihood))
  result[["tsgUnknownDriversTotals"]] <- tsgUnknownDriversTotals

  tsgDrivers = tsgDrivers %>%
    left_join(sampleSomatics, by = "sampleId") %>%
    left_join(tsgUnknownDriversTotals, by = c("gene","impact")) %>%
    mutate(
      coordinate = genomic_coordinates(chromosome, position, ref, alt),
      variant = variant_not_polymorphism(type),
      p_variant_nondriver_snv_single = 1 - ppois(0, sample_SNV / totalSomatics$total_SNV  * gene_non_drivers),
      p_variant_nondriver_snv_multi = 1 - ppois(1, sample_SNV / totalSomatics$total_SNV  * gene_non_drivers),
      p_variant_nondriver_indel_single = 1 - ppois(0, sample_INDEL / totalSomatics$total_INDEL  * gene_non_drivers),
      p_variant_nondriver_indel_multi = 1 - ppois(1, sample_INDEL / totalSomatics$total_INDEL  * gene_non_drivers),
      p_variant_nondriver_single = ifelse(impact == "Indel", p_variant_nondriver_indel_single, p_variant_nondriver_snv_single),
      p_variant_nondriver_multi = ifelse(impact == "Indel", p_variant_nondriver_indel_multi, p_variant_nondriver_snv_multi)
    ) %>%
    group_by(sampleId, gene, impact) %>%
    mutate(
      sameImpact = n() > 1,
      p_variant_nondriver = ifelse(sameImpact, p_variant_nondriver_multi, p_variant_nondriver_single)
    )%>%
    ungroup() %>%
    group_by(sampleId, gene) %>%
    summarise(
      coordinate = paste0(coordinate, collapse = ","),
      variant = paste0(variant, collapse = ","),
      multihit = n() > 1,
      biallelic = any(biallelic),
      hotspot = any(hotspot),
      impact = paste(impact, collapse = ","),
      knownDriver = any(knownDriver),
      driverLikelihood = max(driverLikelihood),
      sameImpact = any(sameImpact),
      p_variant_nondriver = ifelse(sameImpact, max(p_variant_nondriver), prod(p_variant_nondriver)),
      gene_drivers = max(gene_drivers),
      sample_SNV = max(sample_SNV),
      subclonalLikelihood = min(subclonalLikelihood),
      shared = any(shared),
      pHGVS = collapsePHGVS(pHGVS)
    ) %>%
    mutate(
      p_driver = gene_drivers / cohortSize,
      p_driver_variant = p_driver / (p_driver + p_variant_nondriver * (1-p_driver)),
      driverLikelihoodAdjusted = ifelse(driverLikelihood > 0 & driverLikelihood < 1, p_driver_variant, driverLikelihood),
      driver = ifelse(multihit, "Multihit", impact),
      type = 'TSG') %>%
    select(sampleId, gene, coordinate, variant, driver, impact, type, multihit, biallelic, hotspot, subclonalLikelihood, shared, knownDriver, driverLikelihood, driverLikelihoodAdjusted, sample_SNV, pHGVS) %>%
    ungroup()
  result[["tsgDrivers"]] <- tsgDrivers

  return (result)
}

dnds_onco_drivers <- function(sampleSomatics, mutations, expectedDriversPerGene) {
  result = list()
  cohortSize = nrow(sampleSomatics)
  totalSomatics = sampleSomatics %>% ungroup() %>% summarise(total_SNV = sum(sample_SNV), total_INDEL = sum(sample_INDEL))

  oncoMutations = onco_mutations(mutations)
  result[["oncoMutations"]] <- oncoMutations

  oncoDriverRates = dnds_driver_likelihood(oncoMutations, expectedDriversPerGene)
  result[["oncoDriverRates"]] <- oncoDriverRates

  oncoDrivers = oncoMutations %>%
    filter(redundant == F) %>%
    left_join(oncoDriverRates %>% select(gene, impact, driverLikelihood), by = c("gene","impact")) %>%
    mutate(driverLikelihood = ifelse(knownDriver, 1, driverLikelihood)) %>%
    mutate(driverLikelihood = ifelse(impact == "Frameshift", 0, driverLikelihood))

  oncoUnknownDriversTotals = oncoDrivers %>%
    group_by(gene, impact) %>%
    filter(!knownDriver) %>%
    summarise(gene_drivers = sum(driverLikelihood), gene_non_drivers = sum(1 - driverLikelihood))
  result[["oncoUnknownDriversTotals"]] <- oncoUnknownDriversTotals

  oncoDrivers = oncoDrivers %>%
    left_join(sampleSomatics, by = "sampleId") %>%
    left_join(oncoUnknownDriversTotals, by = c("gene","impact")) %>%
    mutate(
      coordinate = genomic_coordinates(chromosome, position, ref, alt),
      variant = variant_not_polymorphism(type),
      p_variant_nondriver = 1 - ppois(0, sample_SNV / totalSomatics$total_SNV  * gene_non_drivers),
      p_driver = gene_drivers / cohortSize,
      p_driver_variant = p_driver / (p_driver + p_variant_nondriver * (1-p_driver)),
      driverLikelihoodAdjusted = ifelse(driverLikelihood > 0 & driverLikelihood < 1, p_driver_variant, driverLikelihood)
    ) %>%
    group_by(sampleId, gene) %>%
    summarise(
      chromosome = first(chromosome),
      position = first(position),
      ref = first(ref),
      alt = first(alt),
      coordinate = paste0(coordinate, collapse = ","),
      variant = paste0(variant, collapse = ","),
      driver = impact,
      nearHotspot = any(nearHotspot),
      hotspot = sum(hotspot) > 0,
      impact = paste(impact, collapse = ","),
      knownDriver = any(knownDriver),
      driverLikelihood = max(driverLikelihood),
      driverLikelihoodAdjusted = max(driverLikelihoodAdjusted),
      sample_SNV = max(sample_SNV),
      subclonalLikelihood = min(subclonalLikelihood),
      shared = any(shared),
      pHGVS = collapsePHGVS(pHGVS)
    ) %>%
    mutate(driver = impact, type = 'ONCO') %>%
    select(sampleId, gene, coordinate, chromosome, position, ref, alt, variant, driver, impact, type, hotspot, nearHotspot, subclonalLikelihood, shared, knownDriver, driverLikelihood, driverLikelihoodAdjusted, sample_SNV, pHGVS) %>%
    ungroup()
  result[["oncoDrivers"]] <- oncoDrivers

  return (result)
}






