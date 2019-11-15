# library(NMF)
library(dplyr)
library(MutationalPatterns)

args <- commandArgs(trailingOnly = TRUE)
print(args)

if(length(args) < 3)
{
  print("Requires arguments 1=sample counts, 2=signatures, 3=outputFile")
  stop()
}

countsFile = args[1]
sigsFile = args[2]
outputFile = args[3]

signatures = as.matrix(read.csv(sigsFile), stringsAsFactors=F)
sampleMatrixData = as.matrix(read.csv(countsFile), stringsAsFactors=F)
sampleId = colnames(sampleMatrixData)[1]


print(sprintf("Fitting sample(%s) SNV counts: input(%s) sigs(%s) outputFile(%s)", sampleId, countsFile, sigsFile, outputFile))

sigFitResult = fit_to_signatures(sampleMatrixData, signatures)
sigFitContributions = as.data.frame(sigFitResult$contribution)

# write out raw contribs and a user-friendly summary

sampleTotal = sum(sampleMatrixData[,1])
sigNames = rownames(sigFitContributions)
sigFitContributions = cbind(sigNames,sigFitContributions)
colnames(sigFitContributions) = c('Signature','SigContribution')
sigFitContributions = sigFitContributions %>% mutate(SigPercent=round(SigContribution/sampleTotal,4),
                                                     SigContribution = round(SigContribution,4))

print(sprintf("Signatures fit complete, writing to %s", outputFile))

write.csv(sigFitContributions %>% filter(SigContribution>0) %>% arrange(-SigPercent),outputFile,row.names=F,quote=F)

