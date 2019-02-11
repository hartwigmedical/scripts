library(dplyr)
library(RMySQL)
library(GenomicRanges)


dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
dbDisconnect(dbProd)
rm(dbProd)

query_indels <- function(dbConnect, sampleId) {
  query = paste0("select * from somaticVariant where type = 'INDEL' and sampleId ='", sampleId, "'")
  return (dbGetQuery(dbConnect, query))
}

query_svs <- function(dbConnect, sampleId) {
  query = paste0("select * from structuralVariant where  type in ('DEL', 'DUP', 'INS') and sampleId ='", sampleId, "'")
  return (dbGetQuery(dbConnect, query))
}

CPCT02030516T
CPCT02010497T


svs = query_svs(dbProd, "CPCT02010497T")
indels = query_indels(dbProd, "CPCT02010497T")
save(svs, indels, file= "~/hmf/analysis/sv/RData/CPCT02010497T.RData")


load(file= "~/hmf/analysis/sv/RData/indelGridssComparisonSampleData.RData")

indels = indels %>% 
  mutate(length = abs(nchar(ref) - nchar(alt))) %>%
  filter(length >= 32, length <= 50)

svs = svs %>%
  mutate(length = endPosition - startPosition) %>%
  filter(length >= 32, length <= 50)

indelRanges = GRanges(indels$chromosome, IRanges(indels$position, indels$position + indels$length))
svsRanges = GRanges(svs$startChromosome, IRanges(svs$startPosition, svs$endPosition))
ol = as.matrix(findOverlaps(indelRanges, svsRanges, type="any", select="all", maxgap = 50))

indels$scope <- "Private"
indels[ol[,1], "scope"]<- "Shared"
indels[ol[, 1], "svId"] <- svs[ol[,2], "id"]

svs$scope <- "Private"
svs[ol[, 2], "scope"]<- "Shared"
svs[ol[, 2], "indelId"] <- indels[ol[, 1], "id"]

exampleJoin = indels %>% full_join(svs %>% select(id, startChromosome, endChromosome, startPosition, endPosition), by = c("svId"="id"))



                                   