library(data.table)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(stringi)
library(devtools)
library(hexbin)



transCatCounts = read.csv(formFilename('~/logs/',sampleId,'category_counts.csv'))
transCombo = read.csv('~/logs/RNA_EXP_TRANS_COMBO_DATA.csv')
View(transCatCounts)
View(transCatCounts %>% filter(is.na(EmFitCount)))
View(transCombo)
View(transCombo %>% filter(is.na(EmFitCount)))



## Transcript Abundance & Expression

transCohort = read.csv('~/data/rna/logs/isfox_ut_138.transcript_distribution.csv')
transCohort = read.csv('~/logs/isfox_test.transcript_distribution.csv')
View(transCohort)

# View(transCohort %>% group_by(Tpm100Bucket=2**round(log(Pct_1.00,2))) %>% count)

transCohort2 = transCohort %>% gather('TpmPctBucket','TPM',4:ncol(transCohort))
View(head(transCohort2,100))
transCohort2 = transCohort2 %>% mutate(Percent=stri_replace_all_fixed(TpmPctBucket,'Pct_',''))
transCohort2 = transCohort2 %>% mutate(Count=as.numeric(as.character(TPM)))

geneCohort = transCohort2 %>% group_by(GeneId,GeneName,Percent) %>% summarise(TPM=sum(TPM)) 
View(geneCohort)
View(head(geneCohort,1000))

geneCohort = geneCohort %>% mutate(LogTPM=ifelse(TPM>0,log(TPM),0))
nrow(geneCohort)

write.csv(geneCohort %>% select(GeneId,GeneName,Percent,TPM),'~/data/rna/logs/isfox_ut_138.gene_distribution.csv',row.names = F,quote = F)
View(geneCohort %>% filter(GeneName=='AHR'))

print(ggplot(geneCohort %>% filter(GeneName=='AHR'), aes(x=Percent, y=TPM))
      + geom_bar(stat = "identity", colour = "black")
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=9))
      + facet_wrap(~GeneName))

print(ggplot(geneCohort %>% filter(GeneName %in% c('ALB','TP53','AHR','KRAS','NRAS')), aes(x=Percent, y=LogTPM))
      + geom_bar(stat = "identity", colour = "black")
      + theme(axis.text.x = element_text(angle = 90, hjust=1,size=7))
      + facet_wrap(~GeneName))


print(ggplot(gcRatioData2 %>% filter(SampleId=='CPCT02210029T') %>% 
               filter(GeneName=='ALL_PERC'|GeneName=='TRANS_FIT_EXPECTED_PERC'), aes(x=RatioDec, y=Count))
      + geom_bar(stat = "identity", colour = "black")
      + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=9))
      + facet_wrap(~GeneName))


genesCohort = read.csv('~/data/rna/logs/isfox_all.gene_distribution.csv')
View(genesCohort %>% filter(GeneName=='HRAS'))
View(head(genesCohort %>% arrange(-Pct_1.00),500))

geneCohort2 = genesCohort %>% gather('FpmPctBucket','FPM',3:ncol(genesCohort))
View(head(geneCohort2,100))
geneCohort2 = geneCohort2 %>% mutate(Percent=stri_replace_all_fixed(FpmPctBucket,'Pct_',''))
geneCohort2 = geneCohort2 %>% mutate(Count=as.numeric(as.character(FPM)))


print(ggplot(geneCohort2 %>% filter(GeneName %in% c('MYC','CCND1','CCNE1','ERBB2','AR','FGFR1','MYC','EGFR','TERT')), aes(x=Percent, y=log(FPM,10)))
      + geom_bar(stat = "identity", colour = "black")
      + theme(axis.text.x = element_text(angle = 90, hjust=1,size=7))
      + facet_wrap(~GeneName))

transCohort = read.csv('~/data/rna/logs/isfox_test.transcript_cohort.csv')
transCohort = read.csv('~/data/rna/logs/isfox_ut_138.transcript_distribution.csv')

View(transCohort)
View(transCohort %>% group_by(TpmBucket=2**round(log(TPM,2))) %>% count)

transCohortUt130 = read.csv('~/logs/isfox_ut_130.transcript_cohort.csv')
View(transCohortUt130 %>% group_by(TpmBucket=2**round(log(TPM,2))) %>% count)
View(transCohortUt130 %>% filter(GeneName=='AHR'))
View(transCohortUt130)
print(ggplot(transCohortUt130 %>% filter(GeneName %in% c('KRAS','NRAS','HRAS','APC')), aes(x=GeneName,y=log(TPM)))
      + geom_violin(scale="count",fill="#6baed6"))






geneTpmSummary = geneDataCombined %>% group_by(GeneName,SampleId) %>% summarise(AdjTPM=sum(AdjTPM)+0.01) %>% spread(SampleId,AdjTPM)
View(geneTpmSummary)

print(ggplot(data=geneTpmSummary %>% filter(`CPCT02010419T`+`CPCT02120066T`>100,`CPCT02010419T`+`CPCT02120066T`<1e6),aes(`CPCT02010419T`,`CPCT02120066T`))
      + geom_hex(bins=120) 
      + geom_abline(slope=1, intercept=0)+scale_x_log10() +scale_y_log10())

ggplot(data=allSamples %>% filter(`CPCT02210029T`+`CPCT02010419T`>100,`CPCT02210029T`+`CPCT02010419T`<1e6),aes(`CPCT02210029T`,`CPCT02010419T`))+geom_hex(bins=120) +
  geom_abline(slope=1, intercept=0)+scale_x_log10() +scale_y_log10()
ggplot(data=allSamples %>% filter(`CPCT02210029T`+`CPCT02010419T`>100,`CPCT02210029T`+`CPCT02010419T`<10000),aes(`CPCT02210029T`,`CPCT02010419T`))+geom_point() + 
  geom_text(aes(label = GeneName), size = 2, nudge_x = 0.25, color = "red") +scale_x_log10() +scale_y_log10()

geneData = read.csv(formFilename('~/data/rna/gcp_data/','CPCT02210029T','gene_data.csv'))
geneData = geneData %>% mutate(EnrichedGene=(GeneName %in% excludedGenes))
View(geneData)
View(geneData %>% filter(EnrichedGene))
View(geneData %>% filter(!EnrichedGene))
nonEnrichedTPMTotal = sum(ifelse(geneData$GeneName %in% excludedGenes,0,geneData$TPM))

geneData = geneData %>% mutate(AdjTPM=round(1e6*TPM/sum(pmin(0.01*nonEnrichedTPMTotal,TPM))),
                               PercOfNonEnrichedTotal=round(AdjTPM/nonEnrichedTPMTotal,6))
View(geneData %>% select(GeneId,GeneName,TPM,AdjTPM,IncludedTPM,SplicedFragments,UnsplicedFragments,FitResiduals,everything()))


# dfExpData = read.csv('~/data/rna/cohort/dana_faber_Expression_MutVsWt_050720.csv')
dfExpData = read.csv('~/data/rna/cohort/isofox_sample_mut_vs_wt_data.csv')
dfExpData = dfExpData %>% filter(SampleId!='')
View(dfExpData)
View(dfExpData %>% filter(SampleId==''))

# correct gene names

# PLPP1 = ENSG00000067113, PPAP2A
# RETREG1 = ENSG00000154153, FAM134B
# SELENOP = ENSG00000250722, SEPP1

dfExpData = dfExpData %>% mutate(GeneName=ifelse(GeneName=='PLPP1','PPAP2A',
                                                 ifelse(GeneName=='RETREG1','FAM134B',
                                                        ifelse(GeneName=='SELENOP','SEPP1',as.character(GeneName)))))

dfGenes = dfExpData %>% group_by(GeneName) %>% count 
dfGenes = merge(dfGenes,ensemblGeneData %>% select(GeneName,GeneId),by='GeneName',all.x=T)
View(dfGenes)
View(dfGenes %>% filter(is.na(GeneId)))
write.csv(dfGenes %>% select(GeneId,GeneName),'~/data/rna/cohort/dana_faber_gene_ids.csv', row.names = F, quote = F)
# View(dfExpData %>% group_by(GeneName) %>% count)

dfSamples = dfExpData %>% group_by(SampleId,DfCancerType=CancerType) %>% count 
View(dfSamples)
View(dfSamples %>% filter(SampleId %in% rnaSampleData$SampleId))
dfRnaSamples = dfSamples %>% filter(SampleId %in% rnaSampleData$SampleId) %>% select(SampleId,DfCancerType)
View(dfRnaSamples)

dfRnaSamples = merge(dfRnaSamples,rnaSampleData %>% select(SampleId,CancerType),by='SampleId',all.x=T)
dfRnaSamples = dfRnaSamples %>% mutate(CancerType=ifelse())

View(dfRnaSamples %>% group_by(CancerType) %>% count)

write.csv(dfRnaSamples %>% select(SampleId,CancerType),'~/data/rna/cohort/dana_faber_rna_samples.csv',row.names = F,quote = F)

dfGeneExpPerc = read.csv('~/data/rna/cohort/isofox_dana_faber.sample_gene_perc_data.csv')
View(dfGeneExpPerc)


allGeneExp = read.csv('~/data/rna/cohort/isfox_all.gene_distribution.csv')
View(allGeneExp)
rm(allGeneExp)

dfCombined = merge(dfExpData %>% filter(SampleId %in% dfRnaSamples$SampleId),dfGeneExpPerc,by=c('SampleId','GeneName'),all.x=T)
View(dfCombined)
View(dfCombined %>% filter(is.na(TPM)|is.na(Cohort)))

# replace HMF sampleId with generic Id
rowIndex = data.frame(as.numeric(as.character(rownames(dfRnaSamples))))
colnames(rowIndex) <- c("RowIndex")
dfRnaSamples = cbind(rowIndex,dfRnaSamples)
dfRnaSamples = dfRnaSamples %>% mutate(HmfSampleId=sprintf('SampleId_%04d',RowIndex))
View(dfRnaSamples)

dfCombinedOut = merge(dfCombined,dfRnaSamples %>% select(SampleId,HmfSampleId),by='SampleId',all.x=T)
View(dfCombinedOut)

write.csv(dfCombinedOut %>% select(HmfSampleId,GeneId,GeneName,CancerType,Cohort,TPM),
          '~/data/rna/cohort/df_mwu_gene_exp_data.csv',row.names = F, quote = F)



##########
## Homozygous Disruptions and Gene Expression
homDis = read.csv('~/data/rna/rna_hom_disruptions.csv')
homDisGenes = homDis %>% group_by(GeneName) %>% count 
View(homDisGenes)
View(homDis %>% group_by(SampleId) %>% count)
homDisGenes = merge(homDisGenes ,ensemblGeneData %>% select(GeneId,GeneName),by='GeneName',all.x=T)
write.csv(homDisGenes %>% select(GeneId,GeneName),'~/data/rna/hom_disruption_gene_ids.csv',row.names = F, quote=F)

homDisExp = read.csv('~/data/rna/cohort/isofox_hom_dis.sample_gene_perc_data.csv')
View(homDisExp)
homDisExp = merge(homDisExp,homDis %>% select(SampleId,GeneName,Driver),by=c('SampleId','GeneName'),all.x=T)
homDisExp = homDisExp %>% mutate(DriverStatus=ifelse(is.na(Driver),'NO_HOM_DIS','HOM_DIS'))
View(homDisExp %>% group_by(GeneName) %>% summarise(Sample=n(),
                                                    WithDriver=sum(DriverStatus=='HOM_DIS')))

write.csv(homDisExp %>% select(-Driver),'~/data/rna/cohort/hom_disruption_gene_expression.csv',row.names = F, quote = F)

geneNames = unique(homDisGenes$GeneName)
View(geneNames)

homDisResults = data.frame(matrix(ncol = 4, nrow = 0))
colnames(homDisResults) = c('GeneName','NoHomDisSampleCount','HomDisSampleCount','PValue') 

for(geneName in geneNames)
{
  geneData = homDisExp %>% filter(GeneName==geneName)
  
  if(nrow(geneData %>% group_by(DriverStatus) %>% count) == 2)
  {
    mwwGt = wilcox.test(TPM ~ DriverStatus, data=geneData, alternative="greater") # refers to Mutation being great instead of WildType
    mwwLt = wilcox.test(TPM ~ DriverStatus, data=geneData, alternative="less")
    pValMin = pmin(mwwGt$p.value, mwwLt$p.value)
    
    rowIndex = nrow(homDisResults) + 1
    homDisResults[rowIndex,1] = geneName
    homDisResults[rowIndex,2] = nrow(geneData %>% filter(DriverStatus=='NO_HOM_DIS'))
    homDisResults[rowIndex,3] = nrow(geneData %>% filter(DriverStatus=='HOM_DIS'))
    homDisResults[rowIndex,4] = pValMin # mwwGt$p.value # pValMin
  }
  else
  {
    print(sprintf("gene(%s) has no hom-dis drivers with RNA", geneName))
  }
}

View(homDisResults)

write.csv(homDisExp %>% select(-Driver),'~/data/rna/cohort/hom_dis_mwu_results.csv',row.names = F, quote = F)


# manual
geneData = homDisExp %>% filter(GeneName=='TP53')
mwwGt = wilcox.test(TPM ~ DriverStatus, data=geneData, alternative="greater") # refers to Mutation being great instead of WildType



# manual calc using Mann-Whitney
geneName = 'ZNF184'
cancerType = 'Bladder'
sampleGroup = dfCombined %>% filter(GeneName==geneName&CancerType==cancerType)
View(sampleGroup)

group1Mut = 'Group1_Mut'
group1Wt = 'Group1_Wt'
group2Mut = 'Group2_Mut'
group2Wt = 'Group2_Wt'

groupVals = sampleGroup %>% filter(Cohort==group1Mut|Cohort==group1Wt) %>% arrange(Cohort)
groupVals = sampleGroup %>% filter(Cohort==group2Mut|Cohort==group2Wt) %>% arrange(Cohort)
View(groupVals)

mwwResult = wilcox.test(TPM ~ Cohort, data=groupVals) 
mwwResult = wilcox.test(TPM ~ Cohort, data=groupVals, alternative="greater")
print(mwwResult$p.value)
mwwResult = wilcox.test(TPM ~ Cohort, data=groupVals %>% arrange(ifelse(grepl('Wt',Cohort),1,2)), alternative="greater")
print(mwwResult$p.value)
mwwResult = wilcox.test(TPM ~ Cohort, data=groupVals %>% arrange(ifelse(grepl('Mut',Cohort),1,2)), alternative="less")
print(mwwResult$p.value)

tmp = data.frame(matrix(ncol = 2, nrow = 0))
colnames(tmp) = c('Cohort','Value')
tmp[1,1] = 'A'
tmp[1,2] = 1
tmp[2,1] = 'A'
tmp[2,2] = 3
tmp[3,1] = 'A'
tmp[3,2] = 5
tmp[4,1] = 'B'
tmp[4,2] = 2
tmp[5,1] = 'B'
tmp[5,2] = 4
tmp[6,1] = 'B'
tmp[6,2] = 6

mwwResult = wilcox.test(Value ~ Cohort, data=tmp, alternative="greater")
print(mwwResult$p.value)
mwwResult = wilcox.test(Value ~ Cohort, data=tmp, alternative="less")
print(mwwResult$p.value)



View(mwwResult)
print(mwwResult$p.value)
print(mwwResult$statistic)

cancerTypes = unique(dfCombined$CancerType)
geneNames = unique(dfCombined$GeneName)
View(geneNames)

dfResults = data.frame(matrix(ncol = 6, nrow = 0))
colnames(dfResults) = c('CancerType','GeneName','Cohort','WtSampleCount','MutSampleCount','PValue') # ,'Result'

for(cancerType in cancerTypes)
{
  for(geneName in geneNames)
  {
    sampleGroup = dfCombined %>% filter(GeneName==geneName&CancerType==cancerType)
    
    if(nrow(sampleGroup) > 0)
    {
      for(i in 1:2)
      {
        wtGroup = ifelse(i==1,group1Wt,group2Wt)
        mutGroup = ifelse(i==1,group1Mut,group2Mut)
        
        groupValues = sampleGroup %>% filter(Cohort==mutGroup|Cohort==wtGroup) %>% arrange
        
        if(nrow(groupValues %>% group_by(Cohort) %>% count) == 2)
        {
          mwwGt = wilcox.test(TPM ~ Cohort, data=groupValues, alternative="greater") # refers to Mutation being great instead of WildType
          # mwwLt = wilcox.test(TPM ~ Cohort, data=groupValues, alternative="less")
          # pValMin = pmin(mwwGt$p.value, mwwLt$p.value)
          
          rowIndex = nrow(dfResults) + 1
          dfResults[rowIndex,1] = cancerType
          dfResults[rowIndex,2] = geneName
          dfResults[rowIndex,3] = sprintf('Group_%d',i)
          dfResults[rowIndex,4] = nrow(groupValues %>% filter(Cohort==wtGroup))
          dfResults[rowIndex,5] = nrow(groupValues %>% filter(Cohort==mutGroup))
          dfResults[rowIndex,6] = mwwGt$p.value # pValMin
          # dfResults[rowIndex,7] = ifelse(mwwGt$p.value < mwwLt$p.value,'Elevated','None')
        }
        else
        {
          print(sprintf("ct=%s g=%s grp=1 has too few cohorts", cancerType, geneName))
        }
      }
    }
  }
}

View(dfResults)

testSampleCounts = dfExpData %>% group_by(CancerType,GeneName) %>% summarise(Group1SampleCount=sum(grepl('Group1',Cohort)),
                                                                             Group2SampleCount=sum(grepl('Group2',Cohort)))
View(testSampleCounts)
dfResultsWithSc = merge(dfResults,testSampleCounts %>% select(CancerType,GeneName,OrigSampleCount=Group1SampleCount),by=c('CancerType','GeneName'),all.x=T)
View(dfResultsWithSc)

write.csv(dfResultsWithSc,'~/data/rna/cohort/df_mwu_results.csv',row.names = F, quote = F)
write.csv(dfCombined,'~/data/rna/cohort/df_mwu_input_data.csv',row.names = F, quote = F)


wilcox.test(Likert ~ Speaker, 
            data=Data)
