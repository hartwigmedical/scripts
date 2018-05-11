detach("package:purple", unload=TRUE)
library(purple)
library(dplyr)


####### LOAD DATA FROM FILE #######
load(file = "~/hmf/RData/reference/geneCopyNumberDeletes.RData")
load(file = "~/hmf/RData/reference/geneCopyNumberAmplifications.RData")
load(file = "~/hmf/RData/reference/canonicalTranscripts.RData")
genes = canonicalTranscripts %>% select(gene, chromosome, start = geneStart, end = geneEnd)

####### EXECUTE ALGORITHM #######
geneCopyNumberDeletes = geneCopyNumberDeletes %>% filter(germlineHetRegions == 0, germlineHomRegions == 0)
geneCopyNumberAmplifications = geneCopyNumberAmplifications %>% filter(germlineHetRegions == 0, germlineHomRegions == 0)

deletionsOutput = copy_number_drivers(genes, geneCopyNumberDeletes, adjacent = "gene"); date()
amplificationOutput = copy_number_drivers(genes, geneCopyNumberAmplifications, adjacent = "arm"); date()

geneCopyNumberAmplificationSummary = amplificationOutput$summary
geneCopyNumberDeletionsSummary = deletionsOutput$summary

save(geneCopyNumberAmplificationSummary, file = "~/hmf/RData/processed/geneCopyNumberAmplificationSummary.RData")
save(geneCopyNumberDeletionsSummary, file = "~/hmf/RData/processed/geneCopyNumberDeletionsSummary.RData")



####### Can also be done in parallel #######
library(doParallel)
no_cores <- 7
cl<-makeCluster(no_cores, type="FORK"); date()
deletionsOutput = copy_number_drivers(genes, geneCopyNumberDeletes, cl = cl, adjacent = "gene"); date()
amplificationOutput = copy_number_drivers(genes, geneCopyNumberAmplifications, cl = cl, adjacent = "arm"); date()
stopCluster(cl)
date()

