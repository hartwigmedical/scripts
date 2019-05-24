library(vcfR)

mnvVcf <- read.vcfR("/Users/jon/hmf/analysis/mnv/mnv.output.vcf", verbose = F)
mnv  = vcfR2tidy(mnvVcf, single_frame = T)$dat %>%
  select(CHROM, POS, REF, ALT, mnvGT = gt_AD) %>%
  filter(nchar(REF) > 1, nchar(REF) == nchar(ALT))

sageVcf <- read.vcfR("/Users/jon/hmf/analysis/mnv/CPCT02030357T.sage.vcf", verbose = F)
sage  = vcfR2tidy(sageVcf, single_frame = F)$fix 


sage = vcfR2tidy(sageVcf, single_frame = T)$dat %>%
  filter(Indiv == "CPCT02030357T") %>%
  select(CHROM, POS, REF, ALT, sageGT = gt_AD)



combined = left_join(sage, mnv, by = c("CHROM", "POS", "REF", "ALT")) %>% mutate(mnvLength = nchar(REF))
combined$sageRefCount <- sapply(strsplit(combined$sageGT, split = ","), function (x) {as.numeric(x[1])})
combined$sageAltCount <- sapply(strsplit(combined$sageGT, split = ","), function (x) {as.numeric(x[2])})
combined$mnvRefCount <- sapply(strsplit(combined$mnvGT, split = ","), function (x) {as.numeric(x[1])})
combined$mnvAltCount <- sapply(strsplit(combined$mnvGT, split = ","), function (x) {as.numeric(x[2])})
combined$sageSupport <- pmin(5, combined$sageAltCount)

combinedSummary = combined %>% group_by(mnvLength, sageSupport) %>% count()
View(combinedSummary)
