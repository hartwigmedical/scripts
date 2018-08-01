svTypes = c("DUP","DEL","BND","INS","INV")
svColours = c("#33a02c","#e31a1c","#1f78b4","#ffff33","#060809")
svColours = setNames(svColours, svTypes)

redGradient = c("#fee5d9", "#fcae91", "#fb6a4a", "#de2d26", "#a50f15")
greenGradient = c("#edf8e9", "#bae4b3", "#74c476", "#31a354", "#006d2c")
blackGradient = c("#f7f7f7", "#cccccc", "#969696", "#636363", "#252525")

sv_length_buckets<-function() {
  return (c("0", "<1k", "1k-10k", "10k-100k", "100k-1M", ">1M"))
}

sv_type_length_buckets<-function() {
  buckets = sv_length_buckets()[-1]
  labels = c("INS_0","BND_0", paste("DEL", buckets, sep = "_"), paste("DUP", buckets, sep = "_"), paste("INV", buckets, sep = "_"))
  return (labels)
}

create_empty_sv_signature<-function() {
  buckets = sv_type_length_buckets()
  type = factor(buckets,levels=buckets,ordered=TRUE)
  return (data.frame(type))
}

sv_overlaps <- function(query, subject, maxgap = -1) {

  require(tidyr)
  require(dplyr)
  require(GenomicRanges)

  queryStartRange <- GRanges(query$startChromosome, IRanges(query$startPosition, query$startPosition))
  subjectStartRange <- GRanges(subject$startChromosome, IRanges(subject$startPosition, subject$startPosition))
  startOverlaps = data.frame(findOverlaps(queryStartRange, subjectStartRange, type="any", select="all", maxgap = maxgap))

  queryEndRange <- GRanges(query$endChromosome, IRanges(query$endPosition, query$endPosition))
  subjectEndRange <- GRanges(subject$endChromosome, IRanges(subject$endPosition, subject$endPosition))
  endOverlaps = data.frame(findOverlaps(queryEndRange, subjectEndRange, type="any", select="all", maxgap = maxgap))

  overlaps = inner_join(startOverlaps, endOverlaps, by = c("queryHits", "subjectHits"))

  overlapQueryData = query[overlaps$queryHits, ] %>%
    mutate(queryHits = overlaps$queryHits) %>%
    select(queryHits, sampleId, startChromosome, endChromosome, startPosition, endPosition, startOrientation, endOrientation, type)

  overlapSubjectData = subject[overlaps$subjectHits, ] %>%
    mutate(subjectHits = overlaps$subjectHits) %>%
    select(subjectHits, subjectStartPosition = startPosition, subjectEndPosition = endPosition, subjectStartOrientation = startOrientation, subjectEndOrientation = endOrientation, subjectType = type)

  overlapsData = bind_cols(overlapQueryData, overlapSubjectData) %>%
    filter(type == subjectType, startOrientation == subjectStartOrientation, endOrientation == subjectEndOrientation) %>%
    select(-subjectType, -subjectStartOrientation, -subjectEndOrientation) %>%
    mutate(startPositionDiff = abs(startPosition - subjectStartPosition), endPositionDiff = abs(endPosition - subjectEndPosition), positionDiff = startPositionDiff + endPositionDiff) %>%
    group_by(startChromosome, endChromosome, startPosition, endPosition,startOrientation,endOrientation,type) %>%
    top_n(1, -positionDiff) %>%
    group_by(queryHits) %>%
    top_n(1, -subjectHits)

  return (overlapsData %>% select(queryHits, subjectHits))
}

sv_signature <- function(variants) {
  require(tidyr)
  require(dplyr)

  empty = create_empty_sv_signature()

  signature = variants %>%
    mutate(
      length = pmax(startPosition, endPosition) - pmin(startPosition, endPosition),
      length = ifelse(type %in% c("BND","INS"), 0, length),
      length = cut(length, right = FALSE, breaks = c(-Inf, 1, 1000, 10000, 1e+05, 1e+06, Inf), labels = sv_length_buckets()),
      type = factor(paste(type, length, sep="_"), levels=sv_type_length_buckets(), ordered = TRUE)
    ) %>%
    group_by(type, sampleId) %>%
    count() %>%
    ungroup() %>%
    spread(sampleId, n, fill = 0)

  result = merge(empty, signature, all.x=TRUE)
  result[is.na(result)] <- 0

  return (result)
}

plot_sv_signature <- function(svSignature) {
  require(ggplot2)

  svSignatureColours = c(svColours[["INS"]], svColours[["BND"]], redGradient, greenGradient, blackGradient)
  svSignatureColours = setNames(svSignatureColours, sv_type_length_buckets())

  p1 <- plot_absolute_contribution(svSignature) +
    theme(axis.text.x = element_text(size=9),legend.text=element_text(size=5),axis.title.y = element_text(size=10))+
    ggtitle("Structural Variant Signatures")+
    scale_fill_manual(name = "", values = svSignatureColours, drop = FALSE) +
    theme(panel.border = element_blank(), axis.ticks.x = element_blank())

  return (p1)
}


# Deprecated
sv_signature_by_scope <- function(variants) {
  DT = sv_signature_data(variants)

  result = dcast(DT, type ~ scope, value.var = "sampleId", fun.aggregate = length)
  result = merge(create_empty_sv_signature(), result, all.x = TRUE)
  result[is.na(result)] <- 0
  return (result)
}

# Deprecated
sv_signature_by_clonality <- function(variants) {
  DT = sv_signature_data(variants)

  result = DT[, .(sample=.N), by=type]
  result = merge(create_empty_sv_signature(), result, all.x = TRUE)
  result[is.na(result)] <- 0
  return (result)
}

# Deprecated
sv_signature_data <- function(variants) {
  DT = data.table(variants)

  ## Length
  DT$length <- pmax(DT$startPosition, DT$endPosition) - pmin(DT$startPosition, DT$endPosition)
  DT$length <- ifelse(DT$type == "BND", 0, DT$length)
  DT$length <- ifelse(DT$type == "INS", 0, DT$length)
  DT$length <- cut(DT$length, right = FALSE, breaks = c(-Inf, 1, 1000, 10000, 1e+05, 1e+06, Inf), labels = sv_length_buckets())
  DT$type <- factor(paste(DT$type, DT$length, sep="_"), levels=sv_type_length_buckets(), ordered = TRUE)

  return (DT)
}



