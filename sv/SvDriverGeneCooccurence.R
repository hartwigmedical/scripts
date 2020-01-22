library(dplyr)
library(tidyr)
library(ggplot2)
library(stringi)
library(devtools)
library(data.table)
library(IRanges)
library(stringr)

###########################
# Driver Gene Co-occurrence
# with SV bucket lengths
# with whole genome duplication


#####
# Cohort sample list

load('~/data/hmf_cohort_may_2019.RData')
View(cohort)
View(highestPurityCohort) # 3524 has multiple biopsy samples removed

sampleCancerTypes = highestPurityCohort %>% select(SampleId=sampleId,CancerType=cancerType)
View(sampleCancerTypes)

#####
# Driver Gene preparation

# For each gene in the driver catalog or SV germline filter for genes with driverlikelihood > 0.5 OR not present (ie exclude 0 to 0.5 since we are unsure about whether there is a driver)

# select d.sampleId as SampleId, gene as Gene, category as Category, driver as Driver, driverLikelihood as DriverLikelihood 
# from driverCatalog d, purity p
# WHERE d.sampleId = p.sampleId and qcStatus = 'PASS' and status <> 'NO_TUMOR' and purity >= 0.195;
driverCatalog = read.csv('~/data/sv/drivers/driver_catalog_20190923.csv')
nrow(driverCatalog)
View(driverCatalog)

View(driverCatalog %>% group_by(cancerType) %>% count())

geneCount = n_distinct(driverCatalog$Gene)
print(geneCount)

driverLikelihoodCutoff=0.8
driverCatalog = driverCatalog %>% mutate(DriverStatus=ifelse(DriverLikelihood>driverLikelihoodCutoff,'TRUE','UNCLEAR'))
View(driverCatalog %>% group_by(Gene) %>% count())

# select g.sampleId as SampleId, gene as Gene, 1 as DriverLikelihood 
# from germlineVariant g, purity p
# WHERE g.sampleId = p.sampleId and qcStatus = 'PASS' and status <> 'NO_TUMOR' and purity >= 0.195;
germlineCatalog = read.csv('~/data/sv/drivers/germline_drivers_20190923.csv')
germlineCatalog$DriverStatus = 'TRUE'
View(germlineCatalog)

allDriverGenes = rbind(germlineCatalog %>% select(SampleId,Gene,DriverStatus),driverCatalog %>% select(SampleId,Gene,DriverStatus))
allDriverGenes = merge(allDriverGenes,sampleCancerTypes,by='SampleId')
nrow(allDriverGenes)
View(allDriverGenes)

# remove any duplicates
allDriverGenes = allDriverGenes %>% group_by(SampleId,Gene,CancerType) %>% summarise(DriverStatus=first(DriverStatus),Count=n())
nrow(allDriverGenes)
View(allDriverGenes)
# View(sampleCancerTypes)
allDriverGenes = allDriverGenes %>% filter(SampleId %in% sampleCancerTypes$SampleId)
nrow(allDriverGenes)
View(allDriverGenes %>% group_by(SampleId) %>% count())


#####
# Drive Gene and SV Type/Length co-occurrence

driverSignatureCoocurrence = read.csv('~/Dropbox/HMF Australia team folder/Structural Variant Analysis/coc_sv_results_extended.csv') 
View(driverSignatureCoocurrence %>% arrange(FDR) %>% filter(CancerType!='All',Category!='SIMPLE_SV',FDR<0.1))

svData = read.csv('~/data/sv/SVA_SVS_PROD.csv')
svData$Length = ifelse(as.character(svData$ChrStart)!=as.character(svData$ChrEnd)|svData$Type=='INS'|svData$ArmEnd!=svData$ArmStart, -1, svData$PosEnd-svData$PosStart)

View(svData %>% group_by(ResolvedType) %>% count())

driverSvData = svData %>% filter(ResolvedType=='DEL'|ResolvedType=='DUP'|ResolvedType=='LINE') %>% 
  mutate(Length=ifelse(Type=='DEL'|Type=='DUP',PosEnd-PosStart,0))

driverSvData = driverSvData %>% filter(SampleId %in% sampleCancerTypes$SampleId)
nrow(driverSvData)

View(driverSvData %>% group_by(ResolvedType) %>% count())

# Categories by type and length:
# Simple DUPs: <=1K bases, >1k - <=80k >80k - <=1M bases, >1M
# Simple DEL <= 500 bases, >500 - <= 5K, >5k <=1M, >1M bases
# LINE

shortDupThreshold=500
medDupThreshold=8e4
longDupThreshold=1e6
shortDelThreshold=500
medDelThreshold=1e4
longDelThreshold=5e5

dgSampleCounts = (driverSvData %>% group_by(SampleId)
                  %>% summarise(DUP_SHORT=sum(ResolvedType=='DUP'&Length<=shortDupThreshold),
                                DUP_MEDIUM=sum(ResolvedType=='DUP'&Length>shortDupThreshold&Length<=medDupThreshold),
                                DUP_LONG=sum(ResolvedType=='DUP'&Length>medDupThreshold&Length<=longDupThreshold),
                                DUP_V_LONG=sum(ResolvedType=='DUP'&Length>longDupThreshold),
                                DEL_SHORT=sum(ResolvedType=='DEL'&Length<=shortDelThreshold),
                                DEL_MEDIUM=sum(ResolvedType=='DEL'&Length>shortDelThreshold&Length<=medDelThreshold),
                                DEL_LONG=sum(ResolvedType=='DEL'&Length>medDelThreshold&Length<=longDelThreshold),
                                DEL_V_LONG=sum(ResolvedType=='DEL'&Length>longDelThreshold),
                                LINE=sum(ResolvedType=='LINE')))

# View(dgSampleCounts)

dgSampleCounts2 = gather(dgSampleCounts,"Category","Count", 2:ncol(dgSampleCounts))
View(dgSampleCounts2)

pdfPlot = (ggplot(data = dgSampleCounts2 %>% filter(Category=='DUP_LONG'&Count>=50), aes(x=reorder(SampleId, -Count), y=Count))
                 + geom_point()
                 + labs(title = "DEL Counts All Samples") + xlab("Samples")
                 + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=7)))

print(pdfPlot)

# work out median rates for each category
categoryData = dgSampleCounts2 %>% group_by(Category) %>% summarise(SampleCount=n(), Median=median(Count), Max=max(Count))
View(categoryData)

# find a set of samples which are enriched (eg. poisson of 99%+ compared to median rate) and not enriched (eg. < median rate).   
dgSampleCounts2 = merge(dgSampleCounts2, categoryData, by='Category', all.x=T)
dgSampleCounts2$PoisProb = ifelse(dgSampleCounts2$Count>dgSampleCounts2$Median,1-ppois(dgSampleCounts2$Count-1, dgSampleCounts2$Median), 1)
dgSampleCounts2$Enriched = ifelse(dgSampleCounts2$Count>=10&dgSampleCounts2$PoisProb<1e-10,'TRUE',ifelse(dgSampleCounts2$Count<=dgSampleCounts2$Median,'FALSE','UNCLEAR'))
View(dgSampleCounts2)
View(dgSampleCounts2 %>% group_by(Category,Enriched) %>% count())

# instead taken from HPC file
# SELECT s.sampleId as SampleId, primaryTumorLocation as CancerType
# FROM hmfpatients.clinical c, sample s, purity p
# WHERE s.sampleId = c.sampleId and s.sampleId = p.sampleId and qcStatus = 'PASS' and status <> 'NO_TUMOR' and purity >= 0.195;
# sampleCancerTypes = read.csv('~/data/sample_cancer_types.csv')

allCategoryData = merge(dgSampleCounts2,sampleCancerTypes,by='SampleId')
View(allCategoryData)
# allCategoryData = dgSampleCounts2


sampleSummary = allCategoryData %>% select(SampleId,CancerType,Category,Enriched)
View(sampleSummary)
View(sampleSummary %>% spread(Category,Enriched))

# Run a fisher exact test to see if the bucket counts are positively or negatively correlated within cancer types and pan-cancer.
View(allDriverGenes)
View(allDriverGenes %>% group_by(Gene) %>% summarise(DriverTrue=sum(DriverStatus=='TRUE'), DriverUnclear=sum(DriverStatus=='UNCLEAR')))

#catGeneResults = calc_category_gene_probs(sampleCancerTypes, allDriverGenes, dgSampleCounts2, T)
# View(catGeneResults)

write.csv(allDriverGenes %>% select(SampleId,CancerType,Gene,DriverStatus), '~/data/sv/drivers/coc_sv_driver_genes.csv', row.names = F, quote = F)
write.csv(allCategoryData %>% select(SampleId,CancerType,Category,Enriched,Count), '~/data/sv/drivers/coc_sv_sample_counts.csv', row.names = F, quote = F)

# load results produced by Java app
cocResults = read.csv('~/data/sv/drivers/STATS_GENE_CATEGORY.csv')

cocResults = add_fdr_by_cancer_type(cocResults)

write.csv(cocResults,'~/data/sv/drivers/sv_driver_gene_cooccurence_paper_hpc.csv',row.names = F, quote = F)

View(cocResults)
View(cocResults %>% filter(WithCatWithGene > ExpectedCount) %>% arrange(FETProb))
View(cocResults %>% filter(CancerType=='Prostate'&Gene=='CDK12'))

spopData = driverCatalog %>% filter(Gene=='SPOP')
spopData = merge(spopData,allCategoryData,by='SampleId',all.x=T)
View(spopData)


## 2-variable test - CancerType vs SV Length Category
# write out data for StatsCalc  to process using its 2-variable co-occurrence routine
twoVarData = allCategoryData %>% group_by(CancerType,Category) %>% summarise(Count=sum(Enriched=='TRUE'))
View(twoVarData)
write.csv(twoVarData,'~/data/sv/drivers/coc_sv_cat_vs_ct_2var_data.csv', row.names = F, quote = F)

# doesn't take the unclear category into consideration so still using the R method below

# confim same results with stat_calcs 2-variable test

# cocResults = read.csv('~/data/sv/drivers/STATS_2VAR.csv')
# cocResults = add_fdr_by_cancer_type(cocResults)
# View(cocResults %>% select(FDR,CountGtExp,everything()))

write.csv(cocResults,'~/data/sv/drivers/sv_driver_gene_cooccurence_paper_hpc.csv',row.names = F, quote = F)




# category vs cancer type co-occurrence
calc_category_cancer_probs<-function(sampleCancerTypes, allCategoryData, logCalcs = T)
{
  allCancerCatProbs = data.frame(matrix(ncol = 14, nrow = 0))
  categoryNames = unique(allCategoryData$Category)
  # categoryNames = c('DUP_LT_100')
  
  cancerTypesList = unique(sampleCancerTypes$CancerType)
  # cancerTypesList = c('Breast')
  
  totalSampleCount = nrow(sampleCancerTypes)
  
  for(cancerTypeStr in cancerTypesList)
  {
    catData = allCategoryData %>% filter(CancerType==cancerTypeStr)
    scCancerType = nrow(sampleCancerTypes %>% filter(CancerType==cancerTypeStr))
    
    print(paste("cancerType=", cancerTypeStr, ", sampleCount=", scCancerType, ", totalSampleCount=", totalSampleCount, sep=''))
    
    cancerCatProbs = data.frame(matrix(ncol = 10, nrow = 0))
    colnames(cancerCatProbs) = c("CatName", "CancerType", "CancerSC", "CategorySC", "WithCatWithCancerSC", "NoCatWithCancerSC", "WithCatNoCancerSC", "NoCatNoCancerSC", "CatAndCancerExp", "FisherET")
    
    for(catName in categoryNames)
    {
      # split up the sample category data into whether enriched, unclear or not enriched
      specificCataData = allCategoryData %>% filter(Category==catName)
      catEnrichedData = specificCataData %>% filter(Enriched=='TRUE')
      catUnclearEnrichedData = specificCataData %>% filter(Enriched=='UNCLEAR')
      
      scWithCat = nrow(catEnrichedData)
      scUnclearCat = nrow(catUnclearEnrichedData)
      scNoCat = totalSampleCount - scWithCat - scUnclearCat
      
      if(logCalcs)
      {
        print(paste("category=", catName, " enriched=", scWithCat, " not-enriched=", scNoCat, " unclear=", scUnclearCat, sep=''))
      }
      
      # now consider co-occurrence
      scWithCatWithCancer = nrow(catEnrichedData %>% filter(CancerType==cancerTypeStr))
      scUnclearCatWithCancer = nrow(catUnclearEnrichedData %>% filter(CancerType==cancerTypeStr))
      scNoCatWithCancer = scCancerType - scWithCatWithCancer - scUnclearCatWithCancer

      scWithCatNoCancer = nrow(catEnrichedData %>% filter(CancerType!=cancerTypeStr))
      scUnclearCatNoCancer = nrow(catUnclearEnrichedData %>% filter(CancerType!=cancerTypeStr))
      
      scNoCancer = (totalSampleCount - scCancerType)
      scNoCatNoCancer = scNoCancer - scWithCatNoCancer - scUnclearCatNoCancer

      cancerSamplesPerc = round(scCancerType/totalSampleCount,6) # % of samples with this cancer type

      if(scWithCatWithCancer < 0 | scNoCatWithCancer < 0 | scWithCatNoCancer < 0 | scNoCatNoCancer < 0)
      {
        print(paste("INVALID category=", catName, " Cancer=", cancerTypeStr, " CancerSC=", scCancerType, " catSC=", scWithCat, " expectedBothSC=", CancerCatExpected,
                    " fetProb=", round(fetProb,10), sep=''))
        
        print(paste("withCancer=", scCancerType, " noCatWithCancer=", scNoCatWithCancer, " withCatWithCancer=", scWithCatWithCancer, " unclearCatWithCancer=", scUnclearCatWithCancer, sep=''))
        print(paste("noCancer=", scNoCancer, " noCatNoCancer=", scNoCatNoCancer, " withCatNoCancer=", scWithCatNoCancer, " unclearCatNoCancer=", scUnclearCatNoCancer, sep=''))
        return (allcancerCatProbs)
      }
      
      fishMatrix = rbind(c(scWithCatWithCancer,scNoCatWithCancer), c(scWithCatNoCancer,scNoCatNoCancer))
      
      # expected count of this Cancer within enriched samples
      cancerCatExpected = round(scWithCat * cancerSamplesPerc, 2)
      
      if(scWithCatWithCancer < cancerCatExpected)
        fetProb = fisher.test(fishMatrix, alternative="less")$p.value
      else
        fetProb = fisher.test(fishMatrix, alternative="greater")$p.value
      
      if(fetProb < 1e-6 | logCalcs)
      {
        print(paste("catName=", catName, ", Cancer=", cancerTypeStr, " CancerSC=", scCancerType, " catSC=", scWithCat, " expectedBothSC=", cancerCatExpected,
                    " fetProb=", round(fetProb,10), sep=''))
        
        print(paste("withCancer=", scCancerType, " noCatWithCancer=", scNoCatWithCancer, " withCatWithCancer=", scWithCatWithCancer, " unclearCatWithCancer=", scUnclearCatWithCancer, sep=''))
        print(paste("noCancer=", scNoCancer, " noCatNoCancer=", scNoCatNoCancer, " withCatNoCancer=", scWithCatNoCancer, " unclearCatNoCancer=", scUnclearCatNoCancer, sep=''))
      }
      
      # colnames(cancerCatProbs) = c("CatName", "CancerType", "CancerSC", "CategorySC", "WithCatWithCancerSC", "NoCatWithCancerSC", "WithCatNoCancerSC", "NoCatNoCancerSC", "CatAndCancerExp", "FisherET")
      
      rowIndex = nrow(cancerCatProbs)+1
      cancerCatProbs[rowIndex,1] = catName
      cancerCatProbs[rowIndex,2] = as.character(cancerTypeStr)
      cancerCatProbs[rowIndex,3] = scCancerType
      cancerCatProbs[rowIndex,4] = scWithCat
      cancerCatProbs[rowIndex,5] = scWithCatWithCancer
      cancerCatProbs[rowIndex,6] = scNoCatWithCancer
      cancerCatProbs[rowIndex,7] = scWithCatNoCancer
      cancerCatProbs[rowIndex,8] = scNoCatNoCancer
      cancerCatProbs[rowIndex,9] = cancerCatExpected
      cancerCatProbs[rowIndex,10] = fetProb
    }

    if(nrow(cancerCatProbs) > 0)
    {
      cancerCatProbs$Count_GT_Exp = cancerCatProbs$WithCatWithCancerSC > cancerCatProbs$CatAndCancerExp
      cancerCatProbs = cancerCatProbs %>% arrange(FisherET)
      
      # set ranking values
      rowIndex = data.frame(as.numeric(as.character(rownames(cancerCatProbs))))
      colnames(rowIndex) <- c("Rank")
      cancerCatProbs = cbind(cancerCatProbs, rowIndex)
      
      cancerCatProbs$TestCount = length(categoryNames) * length(cancerTypesList)
      cancerCatProbs$FDR = cancerCatProbs$FisherET*cancerCatProbs$TestCount/cancerCatProbs$Rank

      allCancerCatProbs = rbind(allCancerCatProbs, cancerCatProbs)
    }
  }
  
  return (allCancerCatProbs)
}

cancerCategoryCOC = calc_category_cancer_probs(sampleCancerTypes, allCategoryData, F)
View(cancerCategoryCOC)
View(cancerCategoryCOC %>% select(FDR,Count_GT_Exp,everything()))

write.csv(cancerCategoryCOC, '~/data/sv/drivers/sv_length_vs_ct_cooccurence.csv', row.names = F, quote = F)



#####
# Whole Genome Duplication vs Driver Genes

sampleWGD = read.csv('~/data/sample_wgd_20190626.csv')
nrow(sampleWGD)
View(sampleWGD)
View(sampleWGD %>% group_by(WGD) %>% count())

sampleWGD = sampleWGD %>% filter(SampleId %in% sampleCancerTypes$SampleId)
nrow(sampleWGD)

nrow(highestPurityCohortSummary)
wgdSamples = sampleWGD %>% filter(WGD==1) %>% select(SampleId)
nrow(wgdSamples)
View(wgdSamples)

colnames(allDriverGenes)
View(head(allDriverGenes,100))

# Manual check

cancerType = 'Breast'
dgData = allDriverGenes %>% filter(CancerType==cancerType)
View(dgData)
View(dgData %>% group_by(SampleId) %>% count())
cancerSamples = sampleCancerTypes %>% filter(CancerType==cancerType)
scCancerType = nrow(cancerSamples)
print(scCancerType)

View(sampleCancerTypes %>% group_by(CancerType) %>% count())
cancerSamples = sampleCancerTypes %>% filter(CancerType==cancerType)
scCancerType = nrow(cancerSamples)
print(scCancerType)

View(dgData %>% group_by(Gene,DriverStatus) %>% count())


gene = 'BRCA2'
geneSamples = dgData %>% filter(Gene==gene)
print(nrow(geneSamples))
scGene = nrow(geneSamples %>% filter(DriverStatus=='TRUE'))
scUnclearGene = nrow(geneSamples %>% filter(DriverStatus=='UNCLEAR'))
noGeneSamples = cancerSamples %>% filter(!(SampleId %in% geneSamples$SampleId))
scNoGene = scCancerType - scGene -scUnclearGene
print(paste('gene=',gene,' scDriverGene=',scGene,' unclear=',scUnclearGene,' noDriver=', scNoGene,sep=''))
print(scUnclearGene)
print(scNoGene)

scWGD = nrow(cancerSamples %>% filter(SampleId %in% wgdSamples$SampleId))
scNoWGD = scCancerType - scWGD
print(paste('WGD=',scWGD,' scNoWGD=',scNoWGD,sep=''))

scWithWGDWithGene = nrow(geneSamples %>% filter(DriverStatus=='TRUE'&SampleId %in% wgdSamples$SampleId))
scWithWGDNoGene = nrow(noGeneSamples %>% filter(SampleId %in% wgdSamples$SampleId))
scNoWGDWithGene = nrow(geneSamples %>% filter(DriverStatus=='TRUE'&!(SampleId %in% wgdSamples$SampleId)))
scNoWGDNoGene = nrow(noGeneSamples %>% filter(!(SampleId %in% wgdSamples$SampleId)))
print(paste('withWGDwithGene=',scWithWGDWithGene,' withWGDnoGene=',scWithWGDNoGene,' noWGDwithGene=',scNoWGDWithGene,' noWGDnoGene=', scNoWGDNoGene,sep=''))

geneSamplesPerc = round(scGene/scCancerType,4) # % of samples with this gene driver
print(geneSamplesPerc)

# expected count of this gene within samples with WGD
expectedCount = round(scWGD * geneSamplesPerc,2)
print(expectedCount)

fishMatrix = rbind(c(scWithWGDWithGene,scNoWGDWithGene), c(scWithWGDNoGene,scNoWGDNoGene))

fetProb = fisher.test(fishMatrix, alternative="less")$p.value
fetProb = fisher.test(fishMatrix, alternative="greater")$p.value
print(fetProb)

# withWGDwithGene = dgData %>% filter(Gene==gene & SampleId %in% wgdSamples$SampleId))

calc_wdg_driver_probs<-function(sampleCancerTypes, wgdSamples, allDriverGenes, logCalcs = T)
{
  resultsCols = c('CancerType','CancerSC','WGDSC','Gene','WithGenesSC','NoGeneSC','WithWGDWithGene','WithWGDNoGene','NoWGDWithGene','NoWGDNoGene','Expected','FisherET')
  colCount = length(resultsCols)
  
  allResults = data.frame(matrix(ncol = colCount+4, nrow = 0))

  cancerTypesList = unique(sampleCancerTypes$CancerType)

  totalSampleCount = nrow(sampleCancerTypes)

  for(cancerTypeStr in cancerTypesList)
  {
    dgData = allDriverGenes %>% filter(CancerType==cancerTypeStr)
    cancerSamples = sampleCancerTypes %>% filter(CancerType==cancerTypeStr)
    
    scCancerType = nrow(cancerSamples)
    
    genesList = unique(dgData$Gene)
    
    print(paste("cancerType=", cancerTypeStr, " sampleCount=", scCancerType, " totalSampleCount=", totalSampleCount, sep=''))
    
    results = data.frame(matrix(ncol = colCount, nrow = 0))
    colnames(results) = resultsCols

    geneCount = length(genesList)

    scWGD = nrow(cancerSamples %>% filter(SampleId %in% wgdSamples$SampleId))
    scNoWGD = scCancerType - scWGD
    
    print(paste("cancer=", cancerTypeStr, " sampleCount=", scCancerType, " wgdCount=", scWGD, " geneCount=", geneCount, sep=''))

    for(gene in genesList)
    {
      geneSamples = dgData %>% filter(Gene==gene)
      scGene = nrow(geneSamples %>% filter(DriverStatus=='TRUE'))
      scUnclearGene = nrow(geneSamples %>% filter(DriverStatus=='UNCLEAR'))
      noGeneSamples = cancerSamples %>% filter(!(SampleId %in% geneSamples$SampleId))
      scNoGene = scCancerType - scGene -scUnclearGene
      
      if(logCalcs)
      {
        print(paste('gene=',gene,' scDriverGene=',scGene,' unclear=',scUnclearGene,' noDriver=', scNoGene,sep=''))
      }

      scWithWGDWithGene = nrow(geneSamples %>% filter(DriverStatus=='TRUE'&SampleId %in% wgdSamples$SampleId))
      scWithWGDNoGene = nrow(noGeneSamples %>% filter(SampleId %in% wgdSamples$SampleId))
      scNoWGDWithGene = nrow(geneSamples %>% filter(DriverStatus=='TRUE'&!(SampleId %in% wgdSamples$SampleId)))
      scNoWGDNoGene = nrow(noGeneSamples %>% filter(!(SampleId %in% wgdSamples$SampleId)))
      
      if(logCalcs)
      {
        print(paste('withWGDwithGene=',scWithWGDWithGene,' withWGDnoGene=',scWithWGDNoGene,' noWGDwithGene=',scNoWGDWithGene,' noWGDnoGene=', scNoWGDNoGene,sep=''))
      }
      
      geneSamplesPerc = round(scGene/scCancerType,4) # % of samples with this gene driver
      # print(geneSamplesPerc)
      
      # expected count of this gene within samples with WGD
      expectedCount = round(scWGD * geneSamplesPerc,4)

      fishMatrix = rbind(c(scWithWGDWithGene,scNoWGDWithGene), c(scWithWGDNoGene,scNoWGDNoGene))
      
      if(scWithWGDWithGene < expectedCount)
        fetProb = fisher.test(fishMatrix, alternative="less")$p.value
      else
        fetProb = fisher.test(fishMatrix, alternative="greater")$p.value

      if(logCalcs)
      {
        print(paste('expected=',expectedCount,' actual=',scWithWGDWithGene,' fetProb=',round(fetProb,4),sep=''))
      }
      
      # colnames(results) = c("CancerType","CancerSC","WGDSC","WithGenesSC","NoGeneSC","WithWGDWithGene","WithWGDNoGene","NoWGDWithGene","NoWGDNoGene","FisherET")
      
      rowIndex = nrow(results)+1
      results[rowIndex,1] = as.character(cancerTypeStr)
      results[rowIndex,2] = gene
      results[rowIndex,3] = scCancerType
      results[rowIndex,4] = scWGD
      results[rowIndex,5] = scGene
      results[rowIndex,6] = scNoGene
      results[rowIndex,7] = scWithWGDWithGene
      results[rowIndex,8] = scWithWGDNoGene
      results[rowIndex,9] = scNoWGDWithGene
      results[rowIndex,10] = scNoWGDNoGene
      results[rowIndex,11] = expectedCount
      results[rowIndex,12] = fetProb
    }
    
    if(nrow(results) > 0)
    {
      results = results %>% mutate(Count_GT_Exp=WithWGDWithGene>Expected) %>% arrange(FisherET)
      
      # set ranking values
      rowIndex = data.frame(as.numeric(as.character(rownames(results))))
      colnames(rowIndex) <- c("Rank")
      results = cbind(results, rowIndex)
      
      results = results %>% mutate(TestCount=geneCount*length(cancerTypesList),
                                   FDR=FisherET*TestCount/Rank)
      
      allResults = rbind(allResults, results)
    }
  }
  
  return (allResults)
}

cancerGeneCounts = allDriverGenes %>% group_by(CancerType,Gene) %>% count() %>% group_by(CancerType) %>% summarise(GeneCount=n())
View(cancerGeneCounts)

wgdDriverGeneCOC = calc_wdg_driver_probs(sampleCancerTypes, wgdSamples, allDriverGenes, F)

wgdDriverGeneCOC = merge(wgdDriverGeneCOC,cancerGeneCounts,by='CancerType',all.x=T)

wgdDriverGeneCOC = wgdDriverGeneCOC %>% mutate(Count_GT_Exp=WithWGDWithGene>Expected)

wgdDriverGeneCOC = wgdDriverGeneCOC %>% mutate(TestCount=GeneCount,
                                               FDR=FisherET*TestCount/Rank)

wgdDriverResults = wgdDriverGeneCOC %>% 
  select(CancerType,Samples=CancerSC,WGD=Gene,Gene=WGDSC,PositiveCorrelation=Count_GT_Exp,Expected,PValue=FisherET,QValue=FDR,
         WithWGDWithGene,WithWGDNoGene,NoWGDWithGene,NoWGDNoGene) %>% 
  mutate(Expected=round(Expected,1)) %>%
  arrange(QValue) %>% filter(QValue<0.1)

View(wgdDriverResults)

colnames(wgdDriverGeneCOC)

write.csv(wgdDriverGeneCOC %>% select(CancerType,Gene=,)
          , '~/data/drivers_wgd_cooccurrence.csv', row.names = F, quote = F)


View(wgdDriverGeneCOC)

wgdDriverGeneCOC = calc_wdg_driver_probs(
  sampleCancerTypes %>% filter(CancerType=='Breast'), 
  wgdSamples,
  allDriverGenes %>% filter(Gene=='BRCA2'),T)

wgdDriverGeneCOC = calc_wdg_driver_probs(
  sampleCancerTypes %>% filter(CancerType=='Breast'), 
  wgdSamples,
  allDriverGenes,T)



#####
# LOH by event
View(tsgDrivers)

tsgChrArmData = tsgDrivers %>% filter(Chromosome!='X'&Chromosome!='Y'&LohType!='DEL_MIN') %>% group_by(Chromosome,Arm) %>% summarise(Count=n())
tsgChrArmData$ChrNumber = as.numeric(as.character(tsgChrArmData$Chromosome))
View(tsgChrArmData %>% arrange(ChrNumber))

print(ggplot(data = tsgChrArmData, aes(x=reorder(Chromosome,ChrNumber),y=Count))
        + geom_col(aes(x=Chromosome, y=Count, fill=Arm)))

print(ggplot() + geom_col(data=tsgChrArmData, aes(x=reorder(Chromosome,ChrNumber), y=Count, fill=Arm), position = "dodge"))


tsgChrArmData2 = tsgDrivers %>% filter(Chromosome!='X'&Chromosome!='Y'&LohType!='DEL_MIN') %>% group_by(Chromosome,Arm,LohType) %>% summarise(Count=n())
tsgChrArmData2$ChrNumber = as.numeric(as.character(tsgChrArmData2$Chromosome))

print(ggplot() 
      + geom_col(data=tsgChrArmData2, aes(x=reorder(Chromosome,ChrNumber), y=Count, fill=Arm), position = "dodge")
      + facet_wrap(~LohType))

# all LOH events
lohEvents = read.csv('~/data/sv/CN_LOH_EVENTS.csv')
rm(lohEvents)
View(head(lohEvents,1000))
lohData = lohEvents %>% filter(IsValid=='true'&Skipped=='false'&Chromosome!='X'&Chromosome!='Y') %>% select(Chromosome,SegStart,SegEnd)
lohData$ChrNumber = as.numeric(as.character(lohData$Chromosome))
lohSummaryData = lohData %>% group_by(Chromosome,ChrNumber) %>% summarise(Count=n())
# lohPlotData$LohType = ifelse(lohPlotData$SegStart=='TELOMERE'&lohPlotData$SegEnd=='TELOMERE','TELO_TELO')
rm(lohData)
rm(lohPlotData)

print(ggplot() + geom_col(data=lohSummaryData, aes(x=reorder(Chromosome,ChrNumber), y=Count)))


# DRIVERS with LOH events

svaDrivers = read.csv('~/data/sv/drivers/SVA_DRIVERS.csv')
View(svaDrivers)

tsgDrivers = svaDrivers %>% filter(Category=='TSG'&(FullyMatched=='true'|DriverType=='MUTATION'))
View(tsgDrivers)
nrow(tsgDrivers)
View(tsgDrivers %>% group_by(SampleId,Gene) %>% count())
View(tsgDrivers %>% group_by(EventType) %>% count())

tsgDrivers = tsgDrivers %>% filter(SampleId %in% sampleCancerTypes$SampleId)
tsgDrivers = merge(tsgDrivers,sampleCancerTypes,by='SampleId')
View(tsgDrivers %>% group_by(CancerType) %>% count())

tsgSampleData = tsgDrivers %>% group_by(CancerType,SampleId,Gene) %>% summarise(LohChr=sum(EventType=='LOH_CHR'),
                                                                                LohArm=sum(EventType=='LOH_ARM'))

tsgSampleData = tsgSampleData %>% mutate(LohType=ifelse(LohChr>0,'LOH_CHR',ifelse(LohArm>0,'LOH_ARM','LOH_FOCAL')))
tsgSampleData = tsgSampleData %>% select(-LohChr,-LohArm)
View(tsgSampleData %>% group_by(LohType) %>% count())

View(tsgSampleData)



# look at co-occurrence for:
# LOH type vs cancer 
# LOH vs gene per cancer type and pan-cancer

# manual validation

cancerType = 'Ovary'
catName = 'LOH'
geneName = 'TP53'

totalSampleCount = n_distinct(tsgSampleData$SampleId)
print(totalSampleCount)

catData = tsgSampleData %>% filter(Gene==geneName)

cancerData = tsgSampleData %>% filter(CancerType==cancerType&Gene==geneName)
scWithCancer = nrow(cancerData %>% group_by(SampleId) %>% count())
print(scWithCancer)



# "gene=TP53catName=TELO_TELO cancer=Ovary cancerSC=96 catSC=473 expected=26.293 fetProb=0"
# "withCancer=96 noCatWithCancer=16 withCatWithCancer=80"
# "noCancer=2868 noCatNoCancer=1238 withCatNoCancer=393"

geneList = unique(tsgSampleData$Gene)
View(geneList)

# write out data for LINX stats to process using its 3-variable co-occurrence routine
View(tsgSampleData %>% select(SampleId,Gene,LohType,CancerType) %>% arrange(Gene,SampleId))

View(tsgSampleData)
write.csv(tsgSampleData %>% select(SampleId,Gene,CancerType,LohType) %>% arrange(Gene,SampleId), 
          '~/data/sv/coc_3var_gene_loh_data.csv', row.names = F, quote = F)

cocGeneLohCancerResults = read.csv('~/data/sv/SVA_STATS_DATA.csv')
View(cocGeneLohCancerResults)

# supplement with FDR info

cocGeneLohCancerResultsExtended = data.frame()
for(geneName in geneList)
{
  cocData = cocGeneLohCancerResults %>% filter(Gene==geneName) %>% arrange(FETProb)
  
  if(nrow(cocData))
  {
    # set ranking values
    rowIndex = data.frame(as.numeric(as.character(rownames(cocData))))
    colnames(rowIndex) <- c("Rank")
    cocData = cbind(cocData, rowIndex)
    cocData$FDR = cocData$FETProb*cocData$TestCount/cocData$Rank
    cocGeneLohCancerResultsExtended = rbind(cocGeneLohCancerResultsExtended,cocData)
  }
}

View(cocGeneLohCancerResultsExtended)

View(cocGeneLohCancerResultsExtended  %>% filter(CountGtExp=='true') %>% group_by(Gene,CancerType) %>% mutate(mostSig=min(FDR)==FDR) %>% ungroup() 
     %>% select(FDR,everything()) %>% filter(mostSig==T,FDR<0.1) %>% arrange(FDR))

write.csv(cocGeneLohCancerResultsExtended, '~/data/sv/coc_gene_loh_results_extended.csv', row.names = F, quote = F)


# previous results
cocGeneLohCancerResultsPrev = read.csv('~/data/sv/coc_gene_loh_results_extended.csv')
View(cocGeneLohCancerResultsPrev)

View(cocGeneLohCancerResultsPrev %>% group_by(Gene,CancerType) %>% mutate(mostSig=min(FDR)==FDR) %>% ungroup() 
     %>% select(FDR,positive=CountGtExp,everything()) %>% filter(mostSig==T,FDR<0.1) %>% arrange(FDR))


# swap gene and cancer type around
cancerList = unique(tsgSampleData$CancerType)
View(cancerList)

write.csv(tsgSampleData %>% select(SampleId,CancerType,Gene,LohType) %>% arrange(CancerType,SampleId), 
          '~/data/sv/coc_3var_gene_loh_by_cancer_data.csv', row.names = F, quote = F)

cocCancerLohGeneResults = read.csv('~/data/sv/SVA_STATS_3VAR.csv')
View(cocCancerLohGeneResults)

cocCancerLohGeneResultsExtended = data.frame()
for(cancerType in cancerList)
{
  cocData = cocCancerLohGeneResults %>% filter(CancerType==cancerType) %>% arrange(FETProb)
  
  if(nrow(cocData))
  {
    # set ranking values
    rowIndex = data.frame(as.numeric(as.character(rownames(cocData))))
    colnames(rowIndex) <- c("Rank")
    cocData = cbind(cocData, rowIndex)
    cocData = cocData %>% mutate(FDR=FETProb*TestCount/Rank)
    cocCancerLohGeneResultsExtended = rbind(cocCancerLohGeneResultsExtended,cocData)
  }
}

View(cocCancerLohGeneResultsExtended)

View(cocCancerLohGeneResultsExtended  %>% filter(CountGtExp=='true') %>% group_by(CancerType,Gene) %>% mutate(mostSig=min(FDR)==FDR) %>% ungroup() 
     %>% select(FDR,everything()) %>% filter(mostSig==T,FDR<0.1) %>% arrange(FDR))




#####
# 2-variable test - CancerType vs LohType
# write out data for StatsCalc to process using its 3-variable co-occurrence routine
View(tsgSampleData %>% group_by(SampleId,LohType,Gene) %>% count())

twoVarData = tsgSampleData %>% group_by(LohType,Gene) %>% summarise(Count=n())
sum(twoVarData$Count)

View(twoVarData)
write.csv(twoVarData, '~/data/sv/coc_2var_gene_loh_data.csv', row.names = F, quote = F)

# manual validation
cancerType = 'Biliary'
lohType = 'LOH_ARM'

scCancer = nrow(tsgSampleData %>% filter(CancerType==cancerType)) # 269
scLohType = nrow(tsgSampleData %>% filter(LohType==lohType)) # 2295
scWithCat1WithCat2 = nrow(tsgSampleData %>% filter(LohType==lohType&CancerType==cancerType)) # 162

sampleWithCat1 = tsgSampleData %>% filter(CancerType==cancerType)
sampleWithCat2 = tsgSampleData %>% filter(LohType==lohType)
scWithCat1NoCat2 = nrow(tsgSampleData %>% filter(CancerType==cancerType&LohType!=lohType)) # 107
scNoCat1WithCat2 = nrow(tsgSampleData %>% filter(LohType==lohType&CancerType!=cancerType)) # 2133
scNoCat1NoCat2 = nrow(tsgSampleData %>% filter(CancerType!=cancerType&LohType!=lohType)) # 5397
totalRecords = nrow(tsgSampleData)
print(totalRecords)
print(totalRecords - scWithCat1WithCat2 - scWithCat1NoCat2 - scNoCat1WithCat2)


View(tsgSampleData %>% filter(CancerType==cancerType & !(SampleId %in% sampleWithCat2$SampleId)))
View(tsgSampleData %>% filter(LohType==lohType & !(SampleId %in% sampleWithCat1$SampleId)))
View(tsgSampleData %>% filter(!(SampleId %in% sampleWithCat2$SampleId) & !(SampleId %in% sampleWithCat1$SampleId)))

cocGeneLohResults = read.csv('~/data/sv/SVA_STATS_2VAR.csv')
View(cocGeneLohResults)

# supplement with FDR info
cocGeneLohResults = cocGeneLohResults %>% arrange(FETProb)

# set ranking values
rowIndex = data.frame(as.numeric(as.character(rownames(cocGeneLohResults))))
colnames(rowIndex) <- c("Rank")
cocGeneLohResults = cbind(cocGeneLohResults, rowIndex)
cocGeneLohResults = cocGeneLohResults %>% mutate(FDR = FETProb*TestCount/Rank)

View(cocGeneLohResults)

# Aug 2019 cohort paper results

View(cocGeneLohResults %>% filter(CountGtExp=='true') %>% group_by(Gene) %>% mutate(mostSig=min(FDR)==FDR) %>% ungroup() 
     %>% select(FDR,everything()) %>% filter(mostSig==T,FDR<0.1) %>% arrange(FDR))

write.csv(cocGeneLohCancerResultsExtended, '~/data/sv/coc_gene_loh_results_extended.csv', row.names = F, quote = F)





calc_loh_gene_coc_v2<-function(cancerTypesList, tsgSampleData, logCalcs = T)
{
  allResults = data.frame()
  
  # cancerTypesList = c('Ovary')

  geneList = unique(tsgSampleData$Gene)

  categoryList = unique(tsgData$LohType)
  categoryCount = length(categoryList)
  
  totalSampleCount = unique(tsgSampleData$SampleId)
  
  for(geneName in geneList)
  {
    subResults = data.frame(matrix(ncol = 9, nrow = 0))
    colnames(subResults) = c("CatType", "CancerType", "SampleCount", "CatSC", "CancerSC", "WithCatWithCancerSC", "NoCatNoCancerSC", "ExpectedSC", "FisherET")
    
    geneData = tsgSampleData %>% filter(Gene==geneName)
    sampleCount = n_distinct(geneData$SampleId)
    
    print(paste("gene=", geneName, ", sampleCount=", sampleCount, sep=''))

    for(catName in categoryList)
    {
      # pan-cancer rates for each category (within this gene)
      scWithCat = nrow(geneData %>% filter(LohType==catName) %>% group_by(SampleId) %>% count())

      for(cancerType in cancerTypesList)
      {
        cancerData = tsgSampleData %>% filter(CancerType==cancerType&Gene==geneName)
        
        # scCancerType = nrow(tsgData %>% group_by(SampleId) %>% count())
    
        catCancerSummary = geneData %>% group_by(SampleId) %>% summarise(WithCat=sum(LohType==catName), WithCancer=sum(CancerType==cancerType))
        
        scNoCatNoCancer = nrow(catCancerSummary %>% filter(WithCat==0&WithCancer==0))
        scWithCatWithCancer = nrow(catCancerSummary %>% filter(WithCat>0&WithCancer>0))
        scWithCatNoCancer = nrow(catCancerSummary %>% filter(WithCat>0&WithCancer==0))
        scNoCatWithCancer = nrow(catCancerSummary %>% filter(WithCat==0&WithCancer>0))
        
        scWithCancer = nrow(catCancerSummary %>% filter(WithCancer>0))
        # scWithCat = nrow(catCancerSummary %>% filter(WithCat>0))
        
        expectedCount = round(scWithCat/sampleCount*scWithCancer,4) 
        
        if(scWithCatWithCancer < 0 | scNoCatWithCancer < 0 | scWithCatNoCancer < 0 | scNoCatNoCancer < 0)
        {
          print(paste("INVALID catName=", catName, " cancer=", cancerType, " cancerSC=", scWithCancer, " catSC=", scWithCat, sep=''))
          
          print(paste("withCancer=", scWithCancer, " noCatWithCancer=", scNoCatWithCancer, " withCatWithCancer=", scWithCatWithCancer, sep=''))
          print(paste("noCancer=", scNoCancer, " noCatNoCancer=", scNoCatNoCancer, " withCatNoCancer=", scWithCatNoCancer, sep=''))
          return (allResults)
        }
        
        fishMatrix = rbind(c(scWithCatWithCancer,scNoCatWithCancer), c(scWithCatNoCancer,scNoCatNoCancer))
        
        if(scWithCatWithCancer < expectedCount)
          fetProb = fisher.test(fishMatrix, alternative="less")$p.value
        else
          fetProb = fisher.test(fishMatrix, alternative="greater")$p.value
        
        if(fetProb < 1e-6 | logCalcs)
        {
          print(paste("gene=", geneName, " catName=", catName, " cancer=", cancerType, " cancerSC=", scWithCancer, " catSC=", scWithCat, " expected=", expectedCount,
                      " fetProb=", round(fetProb,10), sep=''))
          
          print(paste("withCancer=", scWithCancer, " noCatWithCancer=", scNoCatWithCancer, " withCatWithCancer=", scWithCatWithCancer, sep=''))
          print(paste("noCancer=", scNoCancer, " noCatNoCancer=", scNoCatNoCancer, " withCatNoCancer=", scWithCatNoCancer, sep=''))
        }
        
        #     colnames(subResults) = c("CatType", "CancerType", "SampleCount", "CatSC", "CancerSC", "WithCatWithCancerSC", "NoCatNoCancerSC", "ExpectedSC", "FisherET")

        rowIndex = nrow(subResults)+1
        subResults[rowIndex,1] = catName
        subResults[rowIndex,2] = as.character(cancerType)
        subResults[rowIndex,3] = sampleCount
        subResults[rowIndex,4] = scWithCat
        subResults[rowIndex,5] = scWithCancer
        subResults[rowIndex,6] = scWithCatWithCancer
        subResults[rowIndex,7] = scNoCatNoCancer
        subResults[rowIndex,8] = expectedCount
        subResults[rowIndex,9] = fetProb
      }
    }
    
    if(nrow(subResults) > 0)
    {
      subResults$Count_GT_Exp = subResults$WithCatWithCancerSC > subResults$ExpectedSC
      subResults = subResults %>% arrange(FisherET)
      
      # set ranking values
      rowIndex = data.frame(as.numeric(as.character(rownames(subResults))))
      colnames(rowIndex) <- c("Rank")
      subResults = cbind(subResults, rowIndex)
      
      subResults$TestCount = categoryCount * geneCount
      subResults$FDR = subResults$FisherET*subResults$TestCount/subResults$Rank
      
      subResults$Gene = geneName
      
      allResults = rbind(allResults, subResults)
    }
  }
  
  return (allResults)
}

geneLohTypevsCancerCocResults = calc_loh_gene_coc_v2(cancerTypesList, tsgSampleData, F)
View(geneLohTypevsCancerCocResults)

View(allResults)








View(tsgDataSummary)
sum(tsgDataSummary$WithGeneNoCat)

tsgNoGeneSamples = tsgDataSummary %>% filter(WithGene==0)
View(tsgNoGeneSamples)
nrow(tsgNoGeneSamples)
tsgNoCatSamples = tsgDataSummary %>% filter(WithCat==0)
nrow(tsgNoCatSamples)
scNoCatNoGene = nrow(tsgDataSummary %>% filter(WithCat==0&WithGene==0))
print(scNoCatNoGene)
scWithCatWithGene = nrow(tsgDataSummary %>% filter(WithCat>0&WithGene>0))
print(scWithCatWithGene)
nrow(tsgDataSummary %>% filter(WithGeneWithCat>0))

scWithCatNoGene = nrow(tsgDataSummary %>% filter(WithCat>0&WithGene==0))
print(scWithCatNoGene)
View(tsgDataSummary %>% filter(WithCat>0&WithGene==0))

scNoCatWithGene = nrow(tsgDataSummary %>% filter(WithCat==0&WithGene>0))
print(scNoCatWithGene)

scWithGene = nrow(tsgDataSummary %>% filter(WithGene>0))
print(scWithGene)
scWithCat = nrow(tsgDataSummary %>% filter(WithCat>0))
print(scWithCat)


View(cancerTypesPlusAllList)

lohTypeGeneCocResults = calc_loh_gene_coc(cancerTypesPlusAllList, tsgSampleData, F)
View(lohTypeGeneCocResults)

calc_loh_gene_coc<-function(cancerTypesList, tsgSampleData, logCalcs = T)
{
  allResults = data.frame()

  cancerTypesList = cancerTypesPlusAllList

  for(cancerTypeStr in cancerTypesList)
  {
    if(cancerTypeStr!='All')
      tsgData = tsgSampleData %>% filter(CancerType==cancerTypeStr)
    
    scCancerType = nrow(tsgData %>% group_by(SampleId) %>% count())
    
    geneList = unique(tsgData$Gene)
    geneCount = length(geneList)
    
    categoryList = unique(tsgData$LohType)
    categoryCount = length(categoryList)
    
    print(paste("cancerType=", cancerTypeStr, ", sampleCount=", scCancerType, ' lohTypes=', categoryCount, ' geneCount=', geneCount, sep=''))
    
    geneCatProbs = data.frame(matrix(ncol = 9, nrow = 0))
    colnames(geneCatProbs) = c("CatType", "Gene", "SampleCount", "CatSC", "GeneSC", "WithCatWithGeneSC", "NoCatNoGeneSC", "ExpectedSC", "FisherET")
    
    for(catName in categoryList)
    {
      for(geneName in geneList)
      {
        tsgDataSummary = tsgData %>% group_by(SampleId) %>% summarise(WithCat=sum(LohType==catName),
                                                                      WithGene=sum(Gene==geneName))
        
        scNoCatNoGene = nrow(tsgDataSummary %>% filter(WithCat==0&WithGene==0))
        scWithCatWithGene = nrow(tsgDataSummary %>% filter(WithCat>0&WithGene>0))
        scWithCatNoGene = nrow(tsgDataSummary %>% filter(WithCat>0&WithGene==0))
        scNoCatWithGene = nrow(tsgDataSummary %>% filter(WithCat==0&WithGene>0))

        scWithGene = nrow(tsgDataSummary %>% filter(WithGene>0))
        scWithCat = nrow(tsgDataSummary %>% filter(WithCat>0))
        
        expectedCount = round(scWithGene/scCancerType*scWithCat,4) 

        if(scWithCatWithGene < 0 | scNoCatWithGene < 0 | scWithCatNoGene < 0 | scNoCatNoGene < 0)
        {
          print(paste("INVALID catName=", catName, " gene=", geneName, " geneSC=", scWithGene, " catSC=", scWithCat, sep=''))
          
          print(paste("withGene=", scWithGene, " noCatWithGene=", scNoCatWithGene, " withCatWithGene=", scWithCatWithGene, sep=''))
          print(paste("noGene=", scNoGene, " noCatNoGene=", scNoCatNoGene, " withCatNoGene=", scWithCatNoGene, sep=''))
          return (allResults)
        }
        
        fishMatrix = rbind(c(scWithCatWithGene,scNoCatWithGene), c(scWithCatNoGene,scNoCatNoGene))
        
        if(scWithCatWithGene < expectedCount)
          fetProb = fisher.test(fishMatrix, alternative="less")$p.value
        else
          fetProb = fisher.test(fishMatrix, alternative="greater")$p.value
        
        if(fetProb < 1e-6 | logCalcs)
        {
          print(paste("catName=", catName, ", gene=", geneName, " geneSC=", scWithGene, " catSC=", scWithCat, " expectedBothSC=", expectedCount,
                      " fetProb=", round(fetProb,10), sep=''))
          
          print(paste("withGene=", scWithGene, " noCatWithGene=", scNoCatWithGene, " withCatWithGene=", scWithCatWithGene, sep=''))
          print(paste("noGene=", scNoGene, " noCatNoGene=", scNoCatNoGene, " withCatNoGene=", scWithCatNoGene, sep=''))
        }
        
        rowIndex = nrow(geneCatProbs)+1
        geneCatProbs[rowIndex,1] = catName
        geneCatProbs[rowIndex,2] = as.character(geneName)
        geneCatProbs[rowIndex,3] = scCancerType
        geneCatProbs[rowIndex,4] = scWithCat
        geneCatProbs[rowIndex,5] = scWithGene
        geneCatProbs[rowIndex,6] = scWithCatWithGene
        geneCatProbs[rowIndex,7] = scNoCatNoGene
        geneCatProbs[rowIndex,8] = expectedCount
        geneCatProbs[rowIndex,9] = fetProb
      }
    }
    
    if(nrow(geneCatProbs) > 0)
    {
      geneCatProbs$Count_GT_Exp = geneCatProbs$WithCatWithGeneSC > geneCatProbs$ExpectedSC
      geneCatProbs = geneCatProbs %>% arrange(FisherET)
      
      # set ranking values
      rowIndex = data.frame(as.numeric(as.character(rownames(geneCatProbs))))
      colnames(rowIndex) <- c("Rank")
      geneCatProbs = cbind(geneCatProbs, rowIndex)
      
      geneCatProbs$TestCount = categoryCount * geneCount
      geneCatProbs$FDR = geneCatProbs$FisherET*geneCatProbs$TestCount/geneCatProbs$Rank
      
      geneCatProbs$CancerType = cancerTypeStr
      
      allResults = rbind(allResults, geneCatProbs)
    }
  }
  
  return (allResults)
}

lohTypeGeneCocResults = calc_loh_gene_coc(cancerTypesPlusAllList, tsgSampleData, F)
View(lohTypeGeneCocResults)

write.csv(lohTypeGeneCocResults, '~/data/sv/coc_sv_loh_gene.csv', row.names = F, quote = F)


# LOH vs cancer type
calc_loh_cancer_coc<-function(cancerTypesList, tsgSampleData, logCalcs = T)
{
  allResults = data.frame()
  
  totalSamples = n_distinct(tsgSampleData$SampleId)

  categoryList = unique(tsgData$LohType)
  categoryCount = length(categoryList)
  
  # cancerTypesList = c('Prostate')
  
  for(cancerTypeStr in cancerTypesList)
  {
    withCancerTsgData = tsgSampleData %>% filter(CancerType==cancerTypeStr)
    noCancerTsgData = tsgSampleData %>% filter(CancerType!=cancerTypeStr)
    
    scWithCancer = nrow(withCancerTsgData %>% group_by(SampleId) %>% count())
    scNoCancer = nrow(noCancerTsgData %>% group_by(SampleId) %>% count())

    print(paste("cancerType=", cancerTypeStr, ", withCancer=", scWithCancer, ", noCancer=", scNoCancer, ' lohTypes=', categoryCount, sep=''))
    
    cancerCatProbs = data.frame(matrix(ncol = 8, nrow = 0))
    colnames(cancerCatProbs) = c("CatType", "CancerType", "CatSC", "CancerSC", "WithCatWithCancerSC", "NoCatNoCancerSC", "ExpectedSC", "FisherET")
    
    for(catName in categoryList)
    {
      catSummaryData = tsgSampleData %>% group_by(SampleId) %>% summarise(WithCat=sum(LohType==catName))
      scWithCat = nrow(catSummaryData %>% filter(WithCat>0))
      scNoCat = nrow(catSummaryData %>% filter(WithCat==0))
      
      withCancerData = withCancerTsgData %>% group_by(SampleId) %>% summarise(WithCat=sum(LohType==catName))
      noCancerData = noCancerTsgData %>% group_by(SampleId) %>% summarise(WithCat=sum(LohType==catName))

      scWithCatWithCancer = nrow(withCancerData %>% filter(WithCat>0))
      scWithCatNoCancer = nrow(noCancerData %>% filter(WithCat>0))
      scNoCatWithCancer = nrow(withCancerData %>% filter(WithCat==0))
      scNoCatNoCancer = nrow(noCancerData %>% filter(WithCat==0))
      
      expectedCount = round(scWithCancer/totalSamples*scWithCat,4) 
      
      if(scWithCatWithCancer < 0 | scNoCatWithCancer < 0 | scWithCatNoCancer < 0 | scNoCatNoCancer < 0)
      {
        print(paste("INVALID catName=", catName, " cancerSC=", scWithCancer, " catSC=", scWithCat, sep=''))
        
        print(paste("withCancer=", scWithCancer, " noCatWithCancer=", scNoCatWithCancer, " withCatWithCancer=", scWithCatWithCancer, sep=''))
        print(paste("noCancer=", scNoCancer, " noCatNoCancer=", scNoCatNoCancer, " withCatNoCancer=", scWithCatNoCancer, sep=''))
        return (allResults)
      }
      
      fishMatrix = rbind(c(scWithCatWithCancer,scNoCatWithCancer), c(scWithCatNoCancer,scNoCatNoCancer))
      
      if(scWithCatWithCancer < expectedCount)
        fetProb = fisher.test(fishMatrix, alternative="less")$p.value
      else
        fetProb = fisher.test(fishMatrix, alternative="greater")$p.value
      
      if(fetProb < 1e-6 | logCalcs)
      {
        print(paste("catName=", catName, ", cancerSC=", scWithCancer, " catSC=", scWithCat, " expectedSC=", expectedCount,
                    " fetProb=", round(fetProb,10), sep=''))
        
        print(paste("withCancer=", scWithCancer, " noCatWithCancer=", scNoCatWithCancer, " withCatWithCancer=", scWithCatWithCancer, sep=''))
        print(paste("noCancer=", scNoCancer, " noCatNoCancer=", scNoCatNoCancer, " withCatNoCancer=", scWithCatNoCancer, sep=''))
      }
      
      # colnames(cancerCatProbs) = c("CatType", "CancerType", "CatSC", "CancerSC", "WithCatWithCancerSC", "NoCatNoCancerSC", "ExpectedSC", "FisherET")

      rowIndex = nrow(cancerCatProbs)+1
      cancerCatProbs[rowIndex,1] = catName
      cancerCatProbs[rowIndex,2] = as.character(cancerTypeStr)
      cancerCatProbs[rowIndex,3] = scWithCat
      cancerCatProbs[rowIndex,4] = scWithCancer
      cancerCatProbs[rowIndex,5] = scWithCatWithCancer
      cancerCatProbs[rowIndex,6] = scNoCatNoCancer
      cancerCatProbs[rowIndex,7] = expectedCount
      cancerCatProbs[rowIndex,8] = fetProb
    }
    
    if(nrow(cancerCatProbs) > 0)
    {
      cancerCatProbs$Count_GT_Exp = cancerCatProbs$WithCatWithCancerSC > cancerCatProbs$ExpectedSC
      cancerCatProbs = cancerCatProbs %>% arrange(FisherET)
      
      # set ranking values
      rowIndex = data.frame(as.numeric(as.character(rownames(cancerCatProbs))))
      colnames(rowIndex) <- c("Rank")
      cancerCatProbs = cbind(cancerCatProbs, rowIndex)
      
      cancerCatProbs$TestCount = categoryCount * geneCount
      cancerCatProbs$FDR = cancerCatProbs$FisherET*cancerCatProbs$TestCount/cancerCatProbs$Rank
      
      cancerCatProbs$CancerType = cancerTypeStr
      
      allResults = rbind(allResults, cancerCatProbs)
    }
  }
  
  return (allResults)
}

lohTypeCancerCocResults = calc_loh_cancer_coc(cancerTypesList, tsgSampleData, F)
View(lohTypeCancerCocResults)

write.csv(lohTypeCancerCocResults, '~/data/sv/coc_sv_loh_cancer.csv', row.names = F, quote = F)

# manual validation
totalSamples = n_distinct(tsgSampleData$SampleId)
print(totalSamples)

categoryList = unique(tsgData$LohType)
categoryCount = length(categoryList)
print(categoryCount)

cancerTypeStr = 'Prostate'
withCancerTsgData = tsgSampleData %>% filter(CancerType==cancerTypeStr)
noCancerTsgData = tsgSampleData %>% filter(CancerType!=cancerTypeStr)
View(withCancerTsgData)

scWithCancer = nrow(withCancerTsgData %>% group_by(SampleId) %>% count())
scNoCancer = nrow(noCancerTsgData %>% group_by(SampleId) %>% count())

print(paste("cancerType=", cancerTypeStr, ", withCancer=", scWithCancer, ", noCancer=", scNoCancer, ' lohTypes=', categoryCount, sep=''))

catName = 'TELO_TELO'

withCatData = tsgSampleData %>% group_by(SampleId) %>% summarise(WithCat=sum(LohType==catName))
scWithCat = nrow(withCatData %>% filter(WithCat>0))
print(scWithCat)
scNoCat = nrow(withCatData %>% filter(WithCat==0))
print(scNoCat)


withCancerData = withCancerTsgData %>% group_by(SampleId) %>% summarise(WithCat=sum(LohType==catName))
#View(withCancerData)
noCancerData = noCancerTsgData %>% group_by(SampleId) %>% summarise(WithCat=sum(LohType==catName))
#View(noCancerData)

scWithCatWithCancer = nrow(withCancerData %>% filter(WithCat>0))
scWithCatNoCancer = nrow(noCancerData %>% filter(WithCat>0))
scNoCatWithCancer = nrow(withCancerData %>% filter(WithCat==0))
scNoCatNoCancer = nrow(noCancerData %>% filter(WithCat==0))

expectedCount = round(scWithCancer/totalSamples*scWithCat,4) 
    
print(paste("catName=", catName, " cancerSC=", scWithCancer, " catSC=", scWithCat, " expected=", expectedCount, sep=''))
print(paste("withCancer=", scWithCancer, " noCatWithCancer=", scNoCatWithCancer, " withCatWithCancer=", scWithCatWithCancer, sep=''))
print(paste("noCancer=", scNoCancer, " noCatNoCancer=", scNoCatNoCancer, " withCatNoCancer=", scWithCatNoCancer, sep=''))




