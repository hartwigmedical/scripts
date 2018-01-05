library(RMySQL)
library(data.table)
detach("package:purple", unload=TRUE)
library(purple)

####### LOAD DATA FROM DATABASE #######
#prodDB = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
#clinicalData = purple::query_clinical_data(prodDB)
#geneCopyNumberDeletes = purple::query_gene_copy_number_deletes(prodDB)
#geneCopyNumberAmplifactions = purple::query_gene_copy_number_amplifications(prodDB)
#allGenes = query_gene_copy_number(prodDB, allGeneDeletes[1, c('sampleId')])
#genes = allGenes[allGenes$chromosome != 'Y', c("chromosome","start","end", "gene")]
#geneCopyNumberDeletes$cancerType <- sapply(geneCopyNumberDeletes$sampleId, function(x) {clinicalData[match(x, clinicalData$sampleId), c("cancerType")] })
#geneCopyNumberAmplifactions$cancerType <- sapply(geneCopyNumberAmplifactions$sampleId, function(x) {clinicalData[match(x, clinicalData$sampleId), c("cancerType")] })
#save(geneCopyNumberAmplifactions, geneCopyNumberDeletes, genes, file = "~/hmf/geneCopyNumber.RData")
#dbDisconnect(prodDB)
#rm(prodDB)


####### LOAD DATA FROM FILE #######
load("~/hmf/geneCopyNumber.RData")
allGenes = genes
allGeneCopyNumbers = geneCopyNumberDeletes
allGeneCopyNumbers = geneCopyNumberAmplifactions
copyNumberDrivers = copy_number_drivers(allGenes, allGeneCopyNumbers, maxDriversPerChromosome = 3, chromosomes = c(1:2))

summary = copyNumberDrivers$summary
topSummary = summary[summary[, N==max(N), by=.(chromosome)]$V1]

#allRemovedWithLowSD = allRemoved[sd < 0.5]
#topRemovedWithLowSD = allRemovedWithLowSD[allRemovedWithLowSD[, N==max(N), by=.(chromosome)]$V1]
