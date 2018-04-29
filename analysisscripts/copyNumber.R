library(RMySQL)
detach("package:purple", unload=TRUE)
library(purple)
library(dplyr)


####### LOAD DATA FROM DATABASE #######
dbPaper = dbConnect(MySQL(), dbname='hmfpatients_20180418', groups="RAnalysis")
highestPurityCohort = purple::query_highest_purity_cohort(dbPaper)
genes = purple::query_canonical_transcript(dbPaper) %>% select(gene, chromosome, start = geneStart, end = geneEnd)
geneCopyNumberDeletes = purple::query_gene_copy_number_deletes(dbPaper, highestPurityCohort)
geneCopyNumberAmplifications = purple::query_gene_copy_number_amplifications(dbPaper, highestPurityCohort)
dbDisconnect(dbPaper)
rm(dbPaper)

geneCopyNumberDeletes$cancerType <- sapply(geneCopyNumberDeletes$sampleId, function(x) {highestPurityCohort[match(x, highestPurityCohort$sampleId), c("primaryTumorLocation")] })
geneCopyNumberAmplifications$cancerType <- sapply(geneCopyNumberAmplifications$sampleId, function(x) {highestPurityCohort[match(x, highestPurityCohort$sampleId), c("primaryTumorLocation")] })

#save(genes, geneCopyNumberDeletes, file = "~/hmf/RData/geneCopyNumberDeletes.RData")
#save(genes, geneCopyNumberAmplifications, file = "~/hmf/RData/geneCopyNumberAmplifications.RData")


####### LOAD DATA FROM FILE #######
load(file = "~/hmf/RData/geneCopyNumberDeletes.RData")
load(file = "~/hmf/RData/geneCopyNumberAmplifications.RData")

####### EXECUTE ALGORITHM #######
geneCopyNumberDeletes = geneCopyNumberDeletes %>% filter(germlineHetRegions == 0, germlineHomRegions == 0)
geneCopyNumberAmplifications = geneCopyNumberAmplifications %>% filter(germlineHetRegions == 0, germlineHomRegions == 0)


library(doParallel)
no_cores <- 7
cl<-makeCluster(no_cores, type="FORK"); date()
deletionsOutput = copy_number_drivers(genes, geneCopyNumberDeletes, cl = cl,  adjacent = "gene"); date()
amplificationOutput = copy_number_drivers(genes, geneCopyNumberAmplifications, cl = cl,  adjacent = "arm"); date()
stopCluster(cl)
date()

geneCopyNumberAmplificationSummary = amplificationOutput$summary
geneCopyNumberDeletesDriverSummary = deletionsOutput$summary

save(geneCopyNumberAmplificationSummary, file = "~/hmf/RData/geneCopyNumberAmplificationSummary.RData")
save(geneCopyNumberDeletesDriverSummary, file = "~/hmf/RData/geneCopyNumberDeletesDriverSummary.RData")
