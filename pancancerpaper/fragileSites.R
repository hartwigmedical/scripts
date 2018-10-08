load(file = "~/hmf/RData/processed/highestPurityCohortSummary.RData")

query_sv<-function(dbConnect, cohort) {
  sampleIdString = paste("'", cohort$sampleId, "'", collapse = ",", sep = "")
  query = paste(
    "SELECT * FROM structuralVariant where filter = 'PASS'",
    "   AND sampleId in (",sampleIdString, ")",
    sep = " ")
  return (dbGetQuery(dbConnect, query))
}

dbProd = dbConnect(MySQL(), dbname='hmfpatients_20180418', groups="RAnalysis")
hpcStructuralVariants = query_sv(dbProd, highestPurityCohortSummary)
dbDisconnect(dbProd)
rm(dbProd)
#save(hpcStructuralVariants, file =  "~/hmf/RData/reference/hpcStructuralVariants.RData")


load(file =  "~/hmf/RData/reference/canonicalTranscripts.RData")
load(file =  "~/hmf/RData/reference/hpcStructuralVariants.RData")
load(file = "~/hmf/RData/processed/highestPurityCohortSummary.RData")
load(file = "~/hmf/RData/processed/geneCopyNumberDeleteTargets.RData")

deletedGenes = canonicalTranscripts %>%
  filter(gene %in% geneCopyNumberDeleteTargets$target) %>%
  select(gene, chromosome, start = geneStart, end = geneEnd, exonBases, codingBases) %>%
  mutate(length = end - start)

mediumSvDels = hpcStructuralVariants %>% 
  filter(type == 'DEL') %>%
  mutate(length = endPosition - startPosition) %>%
  filter(length > 20000, length < 1000000) %>%
  select(sampleId, chromosome = startChromosome, startPosition, endPosition)

subjectGene = GRanges(deletedGenes$chromosome, IRanges(deletedGenes$start,deletedGenes$end))

queryStart = GRanges(mediumSvDels$chromosome, IRanges(mediumSvDels$startPosition, mediumSvDels$startPosition))
queryEnd = GRanges(mediumSvDels$chromosome, IRanges(mediumSvDels$endPosition, mediumSvDels$endPosition))
olStart = as.matrix(findOverlaps(queryStart, subjectGene, type="any", select="all", maxgap = -1))
olEnd = as.matrix(findOverlaps(queryEnd, subjectGene, type="any", select="all", maxgap = -1))
olAll = bind_rows(data.frame(olStart), data.frame(olEnd)) %>% distinct(queryHits, subjectHits)
olDelsSummary = olAll %>% group_by(subjectHits) %>% count()
rm(queryStart, queryEnd, olStart, olEnd, olAll)

queryStart = GRanges(hpcStructuralVariants$startChromosome, IRanges(hpcStructuralVariants$startPosition, hpcStructuralVariants$startPosition))
queryEnd = GRanges(hpcStructuralVariants$endChromosome, IRanges(hpcStructuralVariants$endPosition, hpcStructuralVariants$endPosition))
olStart = as.matrix(findOverlaps(queryStart, subjectGene, type="any", select="all", maxgap = -1))
olEnd = as.matrix(findOverlaps(queryEnd, subjectGene, type="any", select="all", maxgap = -1))
olAll = bind_rows(data.frame(olStart), data.frame(olEnd)) %>% distinct(queryHits, subjectHits)
olAllSummary = olAll %>% group_by(subjectHits) %>% count()
rm(queryStart, queryEnd, olStart, olEnd, olAll)

deletedGenes[olDelsSummary$subjectHits, 'svDels'] <- olDelsSummary$n
deletedGenes[olAllSummary$subjectHits, 'svAll'] <- olAllSummary$n
deletedGenes[is.na(deletedGenes)] <- 0

fragileGenes = deletedGenes %>%
  mutate(svDelProportion = ifelse(svAll > 0, svDels / svAll, 0)) %>%
  mutate(fragile = (svDelProportion > 0.3 & length > 500000) | svDelProportion > 0.6) %>%
  filter(fragile) %>%
  select(gene_name = gene, fragile)
save(fragileGenes, file = "~/hmf/RData/Processed/fragileGenes.RData")

