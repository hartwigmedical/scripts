

# Standard SNV 96 buckets
create_empty_signature <- function() {
  DF <- data.frame(type = character(), context = character(), stringsAsFactors = FALSE)
  ref_bases = c("C", "T")
  bases = c("A", "C", "G", "T")
  for (ref in ref_bases) {
    for (alt in bases) {
      if (alt != ref) {
        type = paste(ref, alt, sep = ">")
        for (before in bases) {
          for (after in bases) {
            context = paste(before, after, sep = ref)
            DF = rbind(DF, data.frame(type, context, stringsAsFactors = FALSE))
          }
        }
      }
    }
  }
  return(DF)
}



###############
## PCAWG Sigs

pcawgSigs = as.matrix(read.csv("~/data/sigs/pcawg_sigs.csv"))
View(pcawgSigs)

pcawgSigCount = ncol(pcawgSigs)
print(pcawgSigCount)




###############
## PCAWG Sample Counts

pcawgSampleCounts = read.csv("~/data/sigs/pcawg_sample_counts_orig.csv")
View(pcawgSampleCounts[,1:10])
pcawgSampleCounts$Bucket = paste(pcawgSampleCounts$Type,pcawgSampleCounts$Context,sep='_')
colCount = ncol(pcawgSampleCounts)
pcawgSampleCounts2 = cbind(pcawgSampleCounts[,colCount],pcawgSampleCounts[,1:colCount-1])
View(pcawgSampleCounts2[,1:10])
pcawgSampleCounts2 = within(pcawgSampleCounts2, rm(Type))
pcawgSampleCounts2 = within(pcawgSampleCounts2, rm(Context))

print(ncol(pcawgSampleCounts2))

pcawgSampleNames = c('Bucket')
for(i in 2:ncol(pcawgSampleCounts2))
{
  pcawgSampleNames[i] = sprintf('PCAWG_%04d',i-1)
}

View(pcawgSampleNames)
colnames(pcawgSampleCounts2) = pcawgSampleNames
pcawgSampleCounts2 = pcawgSampleCounts2 %>% arrange(Bucket)
View(pcawgSampleCounts2[,1:10])

pcawgSampleCounts2 = within(pcawgSampleCounts2, rm(Bucket))
write.csv(pcawgSampleCounts2,"~/data/sigs/pcawg_matrix_data.csv", row.names = F, quote = F)
rm(pcawgSampleCounts2)



#############
## Load and report on PCAWG fit

get_pcawg_named_sigs<-function()
{
  # Platinum sigs - 48 in total
  pcawgSigsNamesStr = c("01", "02", "03", "04", "05", "06", "07a", "07b", "07c", "07d", "08", "09", "10a", "10b", # 14 in total
                         "11", "12", "13", "14", "15", "16", "17a", "17b", "18", "19", "20", # 11 in total
                         "21", "22", "24", "26", "28", "30", "33", "34", "35", "36", "37", "38", "39", "40", # 14 in total 
                         "5R1", "5R2", "5R3", "5R4", "5R5", "5R6", "5R7", "5R8", "5R9") # 9 in total

                           # "49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59", "60")
  
  return (pcawgSigsNamesStr)
}


## Re-run PCAWG sample counts on PCAWG signatures
pcawgMatrixData = read.csv('~/data/sigs/pcawg_matrix_data.csv')
View(pcawgMatrixData[,1:20])

pcawgSampleCounts = matrix_to_sample_counts(pcawgMatrixData,snvBuckets)
pcawgSampleCancerTypes = pcawgSampleCounts %>% group_by(SampleId) %>% count()
pcawgSampleCancerTypes$CancerType = "Unknown"
pcawgSampleCancerTypes = within(pcawgSampleCancerTypes, rm(n))
View(pcawgSampleCancerTypes)

snvPcawgDataContribs = as.matrix(read.csv('~/data/sigs/logs/snv_pcawg_data_fit_nmf_contribs.csv', stringsAsFactors = F))
snvPcawgDataSigs = as.matrix(read.csv('~/data/sigs/logs/snv_pcawg_data_fit_nmf_sigs.csv', stringsAsFactors = F))

snvPcawgSigsNamesStr = get_pcawg_named_sigs()
length(snvPcawgSigsNamesStr)

View(snvPcawgDataContribs[,1:20])

evaluate_nmf_data("SNV", "pcawg_data_fit", snvPcawgDataSigs, snvPcawgDataContribs, pcawgMatrixData, pcawgSampleCounts,
                  pcawgSampleCancerTypes, snvBuckets, snvPcawgSigsNamesStr, F, F, 0, F)



## Use PCAWG's own contributions
pcawgContribs = read.csv('~/data/sigs/pcawg_sample_contribs_orig.csv', stringsAsFactors = F)
View(pcawgContribs)

# set sampleId and signature back to an index form 0
for(i in 1:nrow(pcawgContribs))
{
  pcawgContribs[i,1] = sprintf('PCAWG_%04d',i)
}

newColumns = c('SampleId')
for(i in 2:ncol(pcawgContribs))
{
  newColumns[i] = i-1
}

View(newColumns)
colnames(pcawgContribs) = newColumns
View(pcawgContribs)

gatherIndex = ncol(pcawgContribs)
pcawgContribs2 = pcawgContribs %>% gather("Signature", "Count", 2:gatherIndex) %>% select(SampleId,Signature,Count)
View(pcawgContribs2)

pcawgContribs3 = pcawgContribs2 %>% spread(SampleId,Count)
pcawgContribs3[is.na(pcawgContribs3)] <- 0
pcawgContribs3$Signature = as.numeric(pcawgContribs3$Signature)
pcawgContribs3 = pcawgContribs3 %>% arrange(Signature)
View(pcawgContribs3[,1:10])
View(pcawgContribs3[,1560:1561])
ncol(pcawgContribs3)
pcawgContribs3 = within(pcawgContribs3, rm(Signature))
write.csv(pcawgContribs3, '~/data/sigs/pcawg_sample_contribs.csv', row.names = F, quote = F)

snvPcawgOwnContribs = as.matrix(read.csv('~/data/sigs/pcawg_sample_contribs.csv', stringsAsFactors = F))
View(snvPcawgOwnContribs[,1:10])
nrow(snvPcawgOwnContribs)

evaluate_nmf_data("SNV", "pcawg_data_own_contribs", snvPcawgDataSigs, snvPcawgOwnContribs, pcawgMatrixData, pcawgSampleCounts,
                  pcawgSampleCancerTypes, snvBuckets, snvPcawgSigsNamesStr, F, F, 0, F)




# PCAWG Residuals own data and own fit

pcawgSampleSigData = get_sig_summary(snvDpPcawgSigs, snvPcawgOwnContribs, pcawgMatrixData, snvPcawgSigsNamesStr, snvBuckets)
View(pcawgSampleSigData)
# pcawgSampleTotals = pcawgSampleSigData %>% group_by(SampleId) %>% 

pcawgSampleSigData$CancerType = "Unknown"

pcawgResSummary = (pcawgSampleSigData %>% group_by(SampleId,CancerType) 
                   %>% summarise(SampleTotal=sum(Count), 
                                 ResidualTotal=round(sum(ifelse(SigName=='Excess',-Count,ifelse(SigName=='Unalloc',Count,0))),0))
                   %>% mutate(ResidualPerc=round(ResidualTotal/SampleTotal,4)))

View(pcawgResSummary)

pcawgResPlot = (ggplot(pcawgResSummary %>% filter(SampleTotal>=1000), aes(CancerType, ResidualPerc))
                + theme(axis.text.x = element_text(angle = 90, hjust=1,size=7))
                + geom_boxplot(varwidth=T, fill="plum") +labs(title="Residuals % by Cancer Type", x="Cancer Type", y="Residuals %"))

print(pcawgResPlot)





###############
## Cosmis Sigs

getCOSMICSignatures <- function()
{
  sp_url = "http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"
  cancer_signatures = read.table(sp_url, sep = "\t", header = T)
  # reorder (to make the order of the trinucleotide changes the same)
  cancer_signatures = cancer_signatures[order(cancer_signatures[, 1]),]
  # only signatures in matrix
  cancer_signatures = as.matrix(cancer_signatures[, 4:33])
  
  return (cancer_signatures)
}

cosmicSignatures = getCOSMICSignatures()
View(cosmicSignatures)

write.csv(cosmicSignatures, "~/data/r_data/snv_cosmic_sigs.csv", quote=F, row.names=F)
rm(cosmicSignatures)
cosmicSigs = as.matrix(read.csv("~/data/sigs/snv_cosmic_sigs.csv"))
View(cosmicSigs)


# fit the samples counts to the cosmic signatures
snvSampleSigs = list()
for(s in snvSampleNames)
  # for(s in sampleCountResults)
{
  print(paste("fitting signatures for ", s, sep=''))
  
  if (!is.null(sampleCountResults[[s]]))
  {
    # we need to slice out only the mutation count columns (delete col 1 and 2)
    res = fit_to_signatures(sampleCountResults[[s]][, -c(1, 2)], cosmicSignatures)
    snvSampleSigs[[s]] <- res$contribution
  }
}
