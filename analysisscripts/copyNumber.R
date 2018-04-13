library(RMySQL)
detach("package:purple", unload=TRUE)
library(purple)
library(dplyr)

####### LOAD DATA FROM DATABASE #######
prodDB = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
cohort = purple::query_highest_purity_cohort(prodDB)
genes = purple::query_canonical_transcript(prodDB) %>% select(gene, chromosome, start = geneStart, end = geneEnd)
geneCopyNumberDeletes = purple::query_gene_copy_number_deletes(prodDB, cohort)
geneCopyNumberAmplifications = purple::query_gene_copy_number_amplifications(prodDB, cohort)
dbDisconnect(prodDB)
rm(prodDB)

geneCopyNumberDeletes$cancerType <- sapply(geneCopyNumberDeletes$sampleId, function(x) {cohort[match(x, cohort$sampleId), c("primaryTumorLocation")] })
geneCopyNumberAmplifications$cancerType <- sapply(geneCopyNumberAmplifications$sampleId, function(x) {cohort[match(x, cohort$sampleId), c("primaryTumorLocation")] })

save(cohort, genes, geneCopyNumberDeletes, file = "~/hmf/RData/geneCopyNumberDeletesData.RData")
save(cohort, genes, geneCopyNumberAmplifications, file = "~/hmf/RData/geneCopyNumberAmplificationsData.RData")


####### LOAD DATA FROM FILE #######
load(file = "~/hmf/RData/geneCopyNumberDeletesData.RData")
load(file = "~/hmf/RData/geneCopyNumberAmplificationsData.RData")


####### EXECUTE ALGORITHM #######
geneCopyNumberDeletes = geneCopyNumberDeletes %>% filter(germlineHetRegions > 0, germlineHomRegions > 0)
geneCopyNumberAmplifications = geneCopyNumberAmplifications %>% filter(germlineHetRegions > 0, germlineHomRegions > 0)


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
