nearHotspot <-function(mutations, distance = 10) {
  hotspots = mutations %>% filter(hotspot > 0) %>% select(chromosome, position) %>% distinct
  hrange <- GRanges(hotspots$chromosome, IRanges(hotspots$position, hotspots$position + distance))
  mrange <- GRanges(mutations$chromosome, IRanges(mutations$position, mutations$position + distance))

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



#tsg_driver_likelihood <- function(tsgAnnotatedMutations, excessTsgRates) {
#  tsgAnnotatedMutations$geneStatus <- ifelse(tsgAnnotatedMutations$redundant, "Redundant", tsgAnnotatedMutations$geneStatus)
#  tsgAnnotatedMutations = left_join(tsgAnnotatedMutations, excessTsgRates %>% select(gene, impact, MultiHitDriverRate, SingleHitDriverRate), by = c("gene", "impact"))
#  tsgAnnotatedMutations$driver = ifelse(tsgAnnotatedMutations$geneStatus == "SingleHit", tsgAnnotatedMutations$SingleHitDriverRate, NA )
#  tsgAnnotatedMutations$driver = ifelse(tsgAnnotatedMutations$geneStatus == "MultiHit", tsgAnnotatedMutations$MultiHitDriverRate, tsgAnnotatedMutations$driver )
#  tsgAnnotatedMutations$driver = ifelse(tsgAnnotatedMutations$geneStatus == "Biallelic", 1, tsgAnnotatedMutations$driver)
#  tsgAnnotatedMutations$driver = ifelse(tsgAnnotatedMutations$hotspot > 0, 1, tsgAnnotatedMutations$driver)
#  tsgAnnotatedMutations$driver = ifelse(tsgAnnotatedMutations$redundant, 0, tsgAnnotatedMutations$driver)
#  return (tsgAnnotatedMutations$driver)
#}


