library(RMySQL)
library(data.table)
detach("package:purple", unload=TRUE)
library(purple)

####### LOAD DATA FROM DATABASE #######
#prodDB = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
#clinicalData = purple::query_clinical_data(prodDB)
#allDeletes = purple::query_gene_copy_number_deletes(prodDB)
#allDeletes$cancerType <- sapply(allDeletes$sampleId, function(x) {clinicalData[match(x, clinicalData$sampleId), c("cancerType")] })
#allGenes = query_all_genes_for_sample(prodDB, allDeletes[1, c('sampleId')])
#allGenes = allGenes[allGenes$chromosome != 'Y', c("chromosome","start","end", "gene")]
#save(allDeletes, allGenes, file = "~/hmf/copyNumberDeletions.RData")
#dbDisconnect(prodDB)
#rm(prodDB)


####### LOAD DATA FROM FILE #######
load("~/hmf/copyNumberDeletions.RData")
copyNumberDeletions = copy_number_deletions(data.table(allGenes), data.table(allDeletes))

str(data.table(allGenes))

#DT = data.table(allDeletes[!is.na(allDeletes$cancerType), ])
#cancerTypes = unique(DT$cancerType)
#cancerTypes = cancerTypes[!is.na(cancerTypes)]

#copyNumberDeletions = copy_number_deletions(allGenes, allDeletes)

#allRemoved = copyNumberDeletions$summary
#topRemoved = allRemoved[allRemoved[, N==max(N), by=.(chromosome)]$V1]

#allRemovedWithLowSD = allRemoved[sd < 0.5]
#topRemovedWithLowSD = allRemovedWithLowSD[allRemovedWithLowSD[, N==max(N), by=.(chromosome)]$V1]

#chrom9 = copyNumberDeletions[["9"]]
#chrom9[[1]]

#chrom1 = copyNumberDeletions[["1"]]
#chrom1[[1]]

#chromX = copyNumberDeletions[["X"]]
#chromX[[1]]

#inside = allGenes[allGenes$start < shift(allGenes$end, type = 'lag') & allGenes$end > shift(allGenes$start, type = 'lag'), ]


#overlapping = allGenes[allGenes$start < shift(allGenes$end, type = 'lag') | allGenes$end > shift(allGenes$start, type = 'lead'), ]
#overlapping[overlapping$chromosome == 'X',]

#shift(allGenes$start, type = 'lag')
#?shift