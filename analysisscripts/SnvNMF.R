library(devtools) #; install_github("im3sanger/dndscv")
library(MutationalPatterns)
library(RMySQL)
library(data.table)
library("NMF")

# plotting
library(grid)
library(gridExtra)
library(ggplot2)

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

standard_mutation <- function(types) {
    types = gsub("G>T", "C>A", types)
    types = gsub("G>C", "C>G", types)
    types = gsub("G>A", "C>T", types)
    types = gsub("A>T", "T>A", types)
    types = gsub("A>G", "T>C", types)
    types = gsub("A>C", "T>G", types)
    return(types)
}

standard_context <- function(raw_type, standard_type, context) {
    x = which(raw_type != standard_type)
    context[x] = reverse(chartr("ATGC", "TACG", context[x]))
    return(context)
}

query_variants <- function(dbConnect, sampleIdStr) {
    query = paste(
        "SELECT sampleId, trinucleotideContext as context, concat(ref,'>', alt) as snv, adjustedVaf * adjustedCopyNumber as ploidy, clonality",
        "FROM somaticVariant",
        "WHERE filter = 'PASS' and length(alt) = length(ref) and length(alt) = 1 and trinucleotideContext not like '%N%'",
        sep = " ")

    if(sampleIdStr != "")
    {
      query = paste(query, " and sampleId in (", sampleIdStr, ")", sep='')
    }

    raw_data = dbGetQuery(dbConnect, query)
    raw_types = raw_data$snv
    standard_types = standard_mutation(raw_types)
    raw_context = raw_data$context
    context = standard_context(raw_types, standard_types, raw_context)

    DT = data.table(
        sample = raw_data$sampleId,
        type = standard_types,
        context = context,
        ploidy = raw_data$ploidy,
        clonality = raw_data$clonality)

    return(DT)
}

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

getSampleIdsStr<-function(samples)
{
  sampleIdsStr = ""

  for(i in 1:nrow(samples))
  {
    sample <- samples[i,]

    if(i > 1)
      sampleIdStr = paste(",'", sample$sampleId, "'", sep="")
    else
      sampleIdStr = paste("'", sample$sampleId, "'", sep="")

    # sampleIdStr = stringi::stri_replace_all_fixed(sampleIdStr, " ", "")
    sampleIdsStr = paste(sampleIdsStr, sampleIdStr)
  }

  return (sampleIdsStr)
}

## NMF functions
sigColours = c("#ff994b", "#463ec0", "#88c928", "#996ffb", "#68b1c0", "#e34bd9", "#106b00", "#d10073", "#98d76a",
               "#6b3a9d", "#d5c94e", "#0072e2", "#ff862c", "#31528d", "#d7003a", "#323233", "#ff4791", "#01837a",
               "#ff748a", "#777700", "#ff86be", "#4a5822", "#ffabe4", "#6a4e03", "#c6c0fb", "#ffb571", "#873659",
               "#dea185", "#a0729d", "#8a392f")

plot_top_snv_samples_by_sig<-function(sampleSigData, sigNames, topN = 50)
{
  sigSamplePlots = list()
  plotIndex = 1
  
  # merge cancer type with sampleId
  sampleSigData = unite(sampleSigData, "SampleId", SampleId, CancerType, sep='_')
  
  for(sigName in sigNames)
  {
    topNSamplesBySig = head(sampleSigData %>% filter(SigName==sigName) %>% arrange(-SnvCount),topN)
    
    # now grab all sig data for these top-N samples
    topNSampleSigData = sampleSigData %>% filter(SampleId %in% topNSamplesBySig$SampleId)
    
    title = paste("Top Samples for Signature ", sigName, sep="")
    
    sampleSigPlot <- (ggplot(topNSampleSigData, aes(x = reorder(SampleId, -SnvCount), y = SnvCount, fill = SigName))
                      + geom_bar(stat = "identity", colour = "black")
                      + labs(x = "", y = "SNV Count by Sample")
                      + scale_fill_manual(values = sigColours)
                      + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
                      + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
                      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=7))
                      + ggtitle(title))
    
    if(plotIndex > 1)
    {
      # remove legend after the first plot
      sampleSigPlot <- sampleSigPlot + theme(legend.position="none")
    }
    
    sigSamplePlots[[plotIndex]] <- sampleSigPlot
    
    if(plotIndex >=4)
    {
      multiplot(plotlist = sigSamplePlots, cols = 2)
      sigSamplePlots = list()
      plotIndex = 1
    }
    else
    {
      plotIndex = plotIndex + 1
    }
  }
  
  if(plotIndex > 1)
  {
    # now print all plots for this cancer type
    multiplot(plotlist = sigSamplePlots, cols = 2)
  }
}

get_snv_sig_stats<-function(sampleSigData) {
  
  # key stats per signature
  sigStats = (sampleSigData %>% group_by(SigName)
              %>% summarise(SampleCount=sum(SnvCount>0),
                            SamplePerc=round(sum(SnvCount>0)/n_distinct(SampleId),2),
                            SnvCount=sum(SnvCount),
                            MaxPercent=round(max(SigPercent),3),
                            AvgPercent=round(sum(ifelse(SigPercent>0.001,SigPercent,0))/sum(SigPercent>0),3),
                            PercGT75=sum(SigPercent>0.75),
                            Perc50_75=sum(SigPercent>0.5&SigPercent<=0.75),
                            Perc25_50=sum(SigPercent>0.25&SigPercent<=0.5),
                            Perc10_25=sum(SigPercent>0.1&SigPercent<=0.25),
                            Perc5_10=sum(SigPercent>0.05&SigPercent<=0.1),
                            PercLT5=sum(SigPercent>0.001&SigPercent<=0.05))
              %>% arrange(-SampleCount))
  
  return (sigStats)
}

plot_snv_sig_samples<-function(sampleSigData, cancerType)
{
  cancerSigData = sampleSigData
  
  if(cancerType != "")
  {
    cancerSigData = cancerSigData %>% filter(CancerType==cancerType)
  }
  
  cancerSampleSigData = cancerSigData %>% arrange(-SampleSnvCount, SampleId) %>% select('SampleId', 'SigName', 'SnvCount')
  
  sigCancerPlots = list()
  plotIndex = 1
  
  cancerSigStats = get_snv_sig_stats(cancerSigData)
  
  if(nrow(cancerSigStats) > 0) {
    
    if(cancerType == "")
    {
      title = "Signature SNV Counts"
    }
    else
    {
      title = paste("Signature SNV Counts for ", cancerType, sep="")
    }
    
    sigStatsPlot = (ggplot(data = cancerSigStats, aes(x = SigName, y = SnvCount, group = 1), fill = SigName)
                    + geom_bar(stat = "identity", colour = "black", size = 0.2)
                    + theme(axis.text.x = element_text(angle = 90, hjust = 1))
                    + ylab("SV Count") + xlab("Signature") + ggtitle(title)
                    + theme(legend.position="none"))
    
    axisRatio = max(cancerSigStats$SnvCount) / max(cancerSigStats$SampleCount) * 0.6
    
    sigStatsPlot <- (sigStatsPlot + geom_line(aes(y = SampleCount*axisRatio, color = "red"))
                     + scale_y_continuous(sec.axis = sec_axis(~.*(1/axisRatio), name = "Sample Count")))
    
    sigCancerPlots[[plotIndex]] <- sigStatsPlot
    plotIndex = plotIndex + 1
  }
  
  if(nrow(cancerSampleSigData) > 0) {
    
    # only plot 50 samples at a time
    numSamples = n_distinct(cancerSampleSigData$SampleId)
    numSigs = n_distinct(cancerSampleSigData$SigName)
    samplesPerPlot = 50
    rowsPerPlot = samplesPerPlot * numSigs
    plotCount = ceiling(numSamples/samplesPerPlot)
    
    for (n in 1:plotCount) {
      rowStart = ((n-1) * rowsPerPlot + 1)
      rowEnd = min((n * rowsPerPlot), numSamples * numSigs)
      
      if(cancerType == "")
      {
        title = "Sig SNV Counts by Sample"
      }
      else
      {
        title = paste("Sig SNV Counts by Sample for ", cancerType, sep="")
      }
      
      sampleSigPlot <- (ggplot(cancerSampleSigData[rowStart:rowEnd,], aes(x = reorder(SampleId, -SnvCount), y = SnvCount, fill = SigName))
                        + geom_bar(stat = "identity", colour = "black")
                        + labs(x = "", y = "SNV Count by Sample")
                        + scale_fill_manual(values = sigColours)
                        + theme_bw() + theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
                        + theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
                        + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=7)))
      
      if(n == 1)
      {
        sampleSigPlot <- sampleSigPlot + ggtitle(title)
      }
      if(n > 1)
      {
        # remove legend after the first plot
        sampleSigPlot <- sampleSigPlot + theme(legend.position="none")
      }
      
      sigCancerPlots[[plotIndex]] <- sampleSigPlot
      
      if(plotIndex >=4)
      {
        multiplot(plotlist = sigCancerPlots, cols = 2)
        sigCancerPlots = list()
        plotIndex = 1
      }
      else
      {
        plotIndex = plotIndex + 1
      }
    }
    
    if(plotIndex > 1)
    {
      # now print all plots for this cancer type
      multiplot(plotlist = sigCancerPlots, cols = 2)
    }
  }
}



### START -> DATA SETUP

load("~/data/highestPurityCohortSummary.RData")
View(highestPurityCohortSummary)

cosmicSignatures = getCOSMICSignatures()
View(cosmicSignatures)

dbProd = dbConnect(MySQL(), user='hmf', password='HMFhmf@1', dbname='hmfpatients', groups = "RAnalysis")
# dbDisconnect(dbProd)

cancerTypes = highestPurityCohortSummary %>% group_by(cancerType) %>% count()
View(cancerTypes)

# download SNVs by cancer type (to keep the queries smaller)
for(i in 1:nrow(cancerTypes))
{
  cancerTypeRow = cancerTypes[i,]
  cancerTypeStr = cancerTypeRow$cancerType

  print(paste(i, ": retrieving SNVs for cancer=", cancerTypeStr, sep=''))

  sampleIdsStr = getSampleIdsStr(highestPurityCohortSummary %>% filter(cancerType==cancerTypeStr))
  # print(sampleIdsStr)
  snvVariants = query_variants(dbProd, sampleIdsStr)
  print(paste(i, ":cancer=", cancerTypeStr, ", SNV count=", nrow(snvVariants),  sep=''))

  ctStr = stringi::stri_replace_all_fixed(cancerTypeStr, '/', '')
  snvFilename = paste("~/logs/r_output/snv_", ctStr, ".csv", sep='')

  write.csv(snvVariants, snvFilename, row.names=F, quote=F)
}

rm(snvVariants)

# now reload these same variants from file into a single massive set
allSnvVariants = data.frame(matrix(ncol = 5, nrow = 0))
allSnvVariants = setNames(allSnvVariants, c("sampleId", "context", "snv", "ploidy", "clonality"))

for(i in 1:nrow(cancerTypes))
{
  cancerTypeRow = cancerTypes[i,]
  cancerTypeStr = cancerTypeRow$cancerType
  print(paste(i, ": loading SNVs for cancer=", cancerTypeStr, sep=''))

  ctStr = stringi::stri_replace_all_fixed(cancerTypeStr, '/', '')
  snvFilename = paste("~/logs/r_output/snv_", ctStr, ".csv", sep='')

  allSnvVariants = rbind(allSnvVariants, read.csv(snvFilename))
}

nrow(allSnvVariants)

svnSampleStats = allSnvVariants %>% group_by(sample) %>% summarise(Count=n())
View(svnSampleStats)

samplesList = unique(highestPurityCohortSummary$sampleId)
View(samplesList)

# load("~/logs/r_output/snvSampleCounts.RData")
# View(result)

# sampleCountResults = result
index = 1
for(s in samplesList)
{
  print(paste(index, ": collating data for ", s, sep=''))
  index = index + 1

  emptySigs = create_empty_signature()

  # allSnvVariants
  sampleCounts = (allSnvVariants %>% filter(sample==s) %>% group_by(type,context)
                  %>% summarise(total=n(),
                                subclonal=sum(clonality=='SUBCLONAL'),
                                clonal=sum(clonality=='CLONAL')))

  sampleResult = merge(emptySigs, sampleCounts, all=T)
  sampleResult[is.na(sampleResult)] <- 0

  sampleCountResults[[s]] = sampleResult
}

save(sampleCountResults, file="~/logs/r_output/snvSampleCounts2.RData")

# release the mammoth
rm(allSnvVariants)

# fit the samples counts to the cosmic signatures
sampleSigs = list()
for(s in samplesList)
{
  print(paste("fitting signatures for ", s, sep=''))

  if (!is.null(result[[s]]))
  {
    # we need to slice out only the mutation count columns (delete col 1 and 2)
    res = fit_to_signatures(sampleCountResults[[s]][, -c(1, 2)], cosmicSignatures)
    sampleSigs[[s]] <- res$contribution
  }
}

View(sampleSigs)
View(cosmicSignatures)

# required format:
# contribution - sigs on the rows (no row numbers?), sampleIds on the column names, data being the contribution value, nothing else
# signatures - as-is I think, 30x sigs for the columns, 96 buckets for the rows
# sigNames - just 1 -> 30 is fine
# sampleIds - ensure naming is the same

contributionTotal = data.frame(matrix(ncol = 0, nrow = 30))
contributionClonal = data.frame(matrix(ncol = 0, nrow = 30))
contributionSubclonal = data.frame(matrix(ncol = 0, nrow = 30))

for(s in samplesList)
{
  # print(paste("extracting sig data for ", s, sep=''))

  if (!is.null(sampleSigs[[s]]))
  {
    sampleSigData = sampleSigs[[s]]
    totalValues = sampleSigData[,1]
    contributionTotal = cbind(contributionTotal, totalValues)

    subclonalValues = sampleSigData[,2]
    contributionSubclonal = cbind(contributionSubclonal, subclonalValues)

    clonalValues = sampleSigData[,3]
    contributionClonal = cbind(contributionClonal, clonalValues)
  }
}

contributionTotal = setNames(contributionTotal, samplesList)
contributionSubclonal = setNames(contributionSubclonal, samplesList)
contributionClonal = setNames(contributionClonal, samplesList)
# View(contributionTotal)


# should now be ready to run through NMF
sampleNames = colnames(contributionTotal)
sigNames = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30")
View(sampleNames)

sampleSigData = svnmf::get_sig_data(cosmicSignatures, contributionTotal, sigNames, sampleNames)
names(sampleSigData)[names(sampleSigData) == 'SvCount'] <- 'SnvCount'
View(sampleSigData)

sigStats = get_snv_sig_stats(sampleSigData)
View(sigStats)

# cancerTypes = highestPurityCohortSummary %>% group_by(cancerType) %>% count()
# View(cancerTypes)

sampleCancerTypes = highestPurityCohortSummary %>% select(sampleId, cancerType)
View(sampleCancerTypes)

sampleSigData = merge(sampleSigData, sampleCancerTypes, by.x="SampleId", by.y="sampleId", all.x=TRUE)
names(sampleSigData)[names(sampleSigData) == 'cancerType'] <- 'CancerType'

sampleSnvCounts = sampleSigData %>% group_by(SampleId) %>% summarise(SampleSnvCount=sum(SnvCount))
sampleSigData = merge(sampleSigData, sampleSnvCounts, by.x="SampleId",by.y="SampleId",all.x=TRUE)
View(sampleSigData)

# save for use in gene-correlations
save(sampleSigData, file="~/logs/r_output/snvSampleSigData.RData")

# bucket and estimate analysis
bucketNames = create_empty_signature()
bucketNames = unite(bucketData, "Bucket", type, context, sep='_')
View(bucketNames)
View(cosmicSignatures)

# get back to raw counts by sample and bucket
sampleBucketCounts = data.frame(matrix(ncol = 0, nrow = 96))

# index = 1
for(s in samplesList)
{
  # print(paste("extracting sig data for ", s, sep=''))
  sampleCounts = sampleCountResults[[s]]
  sampleBucketCounts = cbind(sampleBucketCounts, sampleCounts[,3])
  # index = index + 1
  # if(index > 10)
  #   break
}

sampleBucketCounts = setNames(sampleBucketCounts, samplesList)

View(sampleCounts)
View(sampleBucketCounts)

# nmfMatrixData = svnmf::convert_summary_counts_to_nmf(sampleCounts)
# bucketNames = svnmf::get_bucket_names(nmfMatrixData)
# nmfMatrixData = svnmf::remove_bucket_names(nmfMatrixData)
# View(nmfMatrixData)

nmfEstimate <- nmf(sampleBucketCounts, rank=20:35, method="brunet", nrun=4, seed=123456, .opt='vp6')
plot(nmfEstimate)
svnNmfEstimate = nmfEstimate
save(svnNmfEstimate, file="~/logs/r_output/snvNmfEstimate.RData")
rm(nmfEstimate)

# View(contributionTotal)
# sigBucketData = svnmf::get_bucket_data(cosmicSignatures, contributionTotal, bucketNames)
# View(sigBucketData)
# sigBucketStats = svnmf::get_sig_bucket_stats(sigBucketData)
# View(sigBucketStats)

# report top contributing buckets to aid with signature naming
sigNamesNamed = c("01_CL", "02_LINE", "03_LongDUP", "04_Stressed", "05_DBs", "06_ShortDUP", "07_BND_CN", "08_LongDEL_INV", "09_ShortDEL", "10_MidDUP", "11_MidDEL")

sigBucketTopN = svnmf::get_top_buckets(sigBucketData, sigNamesUnamed, sigNamesNamed)
View(sigBucketTopN)

# key bucket stats
bucketSummaryData = svnmf::get_bucket_stats(sigBucketData)
View(bucketSummaryData)

# least contributing 10 buckets
leastContribBuckets = svnmf::get_least_contrib_buckets(bucketSummaryData)
View(leastContribBuckets)


## Data Output to PDF

pdf(file=paste("~/logs/r_output/snvNmf_01.pdf", sep = ""), height = 14, width = 20)
par(mar=c(1,1,1,1))

# 5. Top 50 samples by signature, but include all other signatures as well
plot_top_snv_samples_by_sig(sampleSigData, sigNames)

# 6. Sigs with Samples by cancer type
plot_snv_sig_samples(sampleSigData, "") # all samples

for(cancerType in cancerTypes$cancerType)
{
  plot_snv_sig_samples(sampleSigData, cancerType)
}

dev.off()





tmp1 = read.csv("~/logs/r_output/snv_CNS.csv")
nrow(tmp1)
View(head(tmp1,10))

snvData = tmp1 %>% filter(clonality!='UNKNOWN')
nrow(snvData)
n_distinct(snvData$sample)


# for (s in samples) {
#
#   # slice for our variants
#   sample_variants = variants[sample == s]
#
#   # TODO: do we want to ignore unknown clonality??
#   total = sample_variants[clonality != 'UNKNOWN', .(total = .N), keyby = .(type, context)]
#   subclonal = sample_variants[clonality == 'SUBCLONAL', .(subclonal = .N), keyby = .(type, context)]
#   clonal = sample_variants[clonality == 'CLONAL', .(clonal = .N), keyby = .(type, context)]
#   clonalA = sample_variants[clonality == 'CLONAL' & ploidy < 1.5, .(clonalLowPloidy = .N), keyby = .(type, context)]
#   clonalB = sample_variants[clonality == 'CLONAL' & ploidy >= 1.5, .(clonalHighPloidy = .N), keyby = .(type, context)]
#
#   # cleanup
#   rm(sample_variants)
#
#   tmp = merge(empty, total, all=TRUE)
#   tmp = merge(tmp, subclonal, all=TRUE)
#   tmp = merge(tmp, clonal, all=TRUE)
#   tmp = merge(tmp, clonalA, all=TRUE)
#   tmp = merge(tmp, clonalB, all=TRUE)
#   tmp[is.na(tmp)] <- 0 # TODO check this works
#   stopifnot(nrow(tmp) == 96)
#
#   result[[s]] = tmp
# }


# list of patients -> data.table of mutation counts
# mutation_vectors = process_variants(snvVariants)

save(cancer_signatures,
     cohort,
     mutation_vectors,
     signatures,
     file = dataFile)
