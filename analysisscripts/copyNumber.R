library(RMySQL)
detach("package:purple", unload=TRUE)
library(purple)

####### LOAD DATA FROM DATABASE #######
prodDB = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
cohort = purple::query_highest_purity_cohort(prodDB)
geneCopyNumberDeletes = purple::query_gene_copy_number_deletes(prodDB, cohort)
geneCopyNumberAmplifications = purple::query_gene_copy_number_amplifications(prodDB, cohort)

allGenes = query_gene_copy_number(prodDB, geneCopyNumberAmplifications[1, c('sampleId')])
geneCopyNumberDeletes$cancerType <- sapply(geneCopyNumberDeletes$sampleId, function(x) {cohort[match(x, cohort$sampleId), c("cancerType")] })
geneCopyNumberAmplifications$cancerType <- sapply(geneCopyNumberAmplifications$sampleId, function(x) {cohort[match(x, cohort$sampleId), c("cancerType")] })

dbDisconnect(prodDB)
rm(prodDB)
save(cohort, allGenes, geneCopyNumberDeletes, file = "~/hmf/RData/geneCopyNumberDeletesData.RData")
save(cohort, allGenes, geneCopyNumberAmplifications, file = "~/hmf/RData/geneCopyNumberAmplificationsData.RData")


####### LOAD DATA FROM FILE #######
load(file = "~/hmf/RData/geneCopyNumberDeletesData.RData")
load(file = "~/hmf/RData/geneCopyNumberAmplificationsData.RData")


####### EXECUTE ALGORITHM #######
library(doParallel)
no_cores <- 7
cl<-makeCluster(no_cores, type="FORK")
date()
amplificationOutput = copy_number_drivers(allGenes, geneCopyNumberAmplifications, cl = cl)
date()
deletionsOutput = copy_number_drivers(allGenes, geneCopyNumberDeletes, cl = cl)
date()

stopCluster(cl)
date()


geneCopyNumberAmplificationSummary = amplificationOutput$summary
save(geneCopyNumberAmplificationSummary, file = "~/hmf/RData/geneCopyNumberAmplificationSummary.RData")

geneCopyNumberDeletesDriverSummary = deletionsOutput$summary
save(geneCopyNumberDeletesDriverSummary, file = "~/hmf/RData/geneCopyNumberDeletesDriverSummary.RData")
