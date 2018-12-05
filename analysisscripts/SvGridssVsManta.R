library(tidyr)
library(dplyr)
library(GenomicRanges)
library(MutationalPatterns)
detach("package:purple", unload=TRUE)
library(purple)
library(ggplot2)

combined_signature <- function(GT1, GT2, MT1, MT2) {
  sampleFactors = c("GT1","GT2","MT1","MT2",
                    "GT1,GT2","GT1,MT1","GT1,MT2","GT2,MT1","GT2,MT2","MT1,MT2",
                    "GT1,GT2,MT1","GT1,GT2,MT2","GT1,MT1,MT2","GT2,MT1,MT2",
                    "GT1,GT2,MT1,MT2")

  GT1Contrib = GT1 %>%
    unite(scope, GT1, GT2, MT1, MT2, sep = ",") %>%
    mutate(scope = gsub(',NA',"", scope)) %>%
    dplyr::select(sampleId = scope, startPosition, endPosition, type)

  GT2Contrib = GT2 %>%
    filter(is.na(GT1)) %>%
    unite(scope, GT2, MT1, MT2, sep = ",") %>%
    mutate(scope = gsub(',NA',"", scope)) %>%
    dplyr::select(sampleId = scope, startPosition, endPosition, type)

  MT1Contrib = MT1 %>%
    filter(is.na(GT1), is.na(GT2)) %>%
    unite(scope, MT1, MT2, sep = ",") %>%
    mutate(scope = gsub(',NA',"", scope)) %>%
    dplyr::select(sampleId = scope, startPosition, endPosition, type)

  MT2Contrib = MT2 %>%
    filter(is.na(GT1), is.na(GT2), is.na(MT1)) %>%
    mutate(scope = "MT2") %>%
    dplyr::select(sampleId = scope, startPosition, endPosition, type)

  total = bind_rows(GT1Contrib, GT2Contrib) %>% bind_rows(MT1Contrib) %>% bind_rows(MT2Contrib) %>%  mutate(sampleId = factor(sampleId, rev(sampleFactors)))
  return (sv_signature(total))
}


combined_scope <- function(GT1, GT2, MT1, MT2) {
  result = list()

  GT1 = GT1 %>% mutate(GT1 = "GT1", GT2 = NA, MT1 = NA, MT2 = NA, GT1ClusterCount = ClusterCount, GT2ClusterCount = NA, MT1ClusterCount = NA, MT2ClusterCount = NA)
  GT2 = GT2 %>% mutate(GT1 = NA, GT2 = "GT2", MT1 = NA, MT2 = NA, GT1ClusterCount = NA, GT2ClusterCount = ClusterCount, MT1ClusterCount = NA, MT2ClusterCount = NA)
  MT1 = MT1 %>% mutate(GT1 = NA, GT2 = NA, MT1 = "MT1", MT2 = NA, GT1ClusterCount = NA, GT2ClusterCount = NA, MT1ClusterCount = ClusterCount, MT2ClusterCount = NA)
  MT2 = MT2 %>% mutate(GT1 = NA, GT2 = NA, MT1 = NA, MT2 = "MT2", GT1ClusterCount = NA, GT2ClusterCount = NA, MT1ClusterCount = NA, MT2ClusterCount = ClusterCount)

  scope = append_scope("GT1", "GT2", GT1, GT2); GT1 = scope[["S1"]]; GT2 = scope[["S2"]]
  scope = append_scope("GT1", "MT1", GT1, MT1); GT1 = scope[["S1"]]; MT1 = scope[["S2"]]
  scope = append_scope("GT1", "MT2", GT1, MT2); GT1 = scope[["S1"]]; MT2 = scope[["S2"]]
  scope = append_scope("GT2", "MT1", GT2, MT1); GT2 = scope[["S1"]]; MT1 = scope[["S2"]]
  scope = append_scope("GT2", "MT2", GT2, MT2); GT2 = scope[["S1"]]; MT2 = scope[["S2"]]
  scope = append_scope("MT1", "MT2", MT1, MT2); MT1 = scope[["S1"]]; MT2 = scope[["S2"]]

  result[["GT1"]] <- GT1
  result[["GT2"]] <- GT2
  result[["MT1"]] <- MT1
  result[["MT2"]] <- MT2

  return (result)
}

append_scope <- function(sample1Name, sample2Name, S1, S2) {
  S1_S2 = sv_overlaps(S1, S2, maxgap = 100)

  S1[S1_S2$queryHits, sample2Name] <- sample2Name
  S1[S1_S2$queryHits, paste0(sample2Name, "ClusterCount")] <- S2[S1_S2$subjectHits, "ClusterCount"]

  S2[S1_S2$subjectHits, sample1Name] <- sample1Name
  S2[S1_S2$subjectHits, paste0(sample1Name, "ClusterCount")] <- S1[S1_S2$queryHits, "ClusterCount"]

  result = list()
  result[["S1"]] <- S1
  result[["S2"]] <- S2

  return (result)
}


old_column_names <- function(variants) {
  variants %>%
    mutate(
      sampleId = SampleId,
      startChromosome = ChrStart,
      endChromosome = ChrEnd,
      startPosition = PosStart,
      endPosition = PosEnd,
      startOrientation = OrientStart,
      endOrientation = OrientEnd,
      type = Type) %>%
    dplyr::select(-SampleId, -ChrStart, - ChrEnd, -PosStart, -PosEnd, -OrientStart, -OrientEnd) %>%
    dplyr::select(sampleId, startChromosome, endChromosome, startPosition,endPosition,startOrientation,endOrientation,type, everything())
}


######################## LOADING STARTS HERE ###############################
gridssCohortVariants = old_column_names(read.csv('/Users/peterpriestley/hmf/analyses/cluster/CLUSTER_GRIDSS_LATEST.csv', header = T, stringsAsFactors = F))
gridssCohort = gridssCohortVariants %>%
  dplyr::select(sampleId) %>%
  distinct() %>%
  mutate(patientId = substring(sampleId, 1, 12)) %>%
  group_by(patientId) %>%
  filter(n() > 1)

prodCohortVariants = old_column_names(read.csv('/Users/peterpriestley/hmf/analyses/gridssAssessment/CLUSTER_V23.csv', header = T, stringsAsFactors = F)) %>%
  filter(sampleId %in% gridssCohort$sampleId)
#save(prodCohortVariants, file = "/Users/jon/hmf/gridss/RData/mantaVariants.RData")
#load(file = "/Users/jon/hmf/gridss/RData/mantaVariants.RData")
#prodCohortVariants = old_column_names(read.csv('/Users/jon/hmf/gridss/clustering/CLUSTER_V24_COMPACT.csv', header = T, stringsAsFactors = F)) %>%
#  filter(sampleId %in% gridssCohort$sampleId)
prodCohort = prodCohortVariants %>%
  dplyr::select(sampleId) %>%
  distinct() %>%
  mutate(patientId = substring(sampleId, 1, 12)) %>%
  group_by(patientId) %>%
  filter(n() > 1)

cohort = inner_join(gridssCohort, prodCohort, by = c("patientId", "sampleId")) %>%
  group_by(patientId) %>%
  filter(n() ==2) %>%
  arrange(sampleId) %>%
  mutate(scope = paste0("Sample", row_number())) %>%
  spread(scope, sampleId)

rm(gridssCohort, prodCohort)
#View(gridssCohortVariants %>% filter(!(type=='DEL'&(InexactHOEnd-InexactHOStart>=6|nchar(Homology)<6)&endPosition-startPosition<1000)) %>% group_by(sampleId) %>% count())


#View(gridssCohortVariants %>% filter(sampleId=='CPCT02010255T'))

#### TEMP FILTER for SHORT DELS
#View(gridssCohortVariants %>% filter(grepl('asm',AsmbStart)|grepl('asm',AsmbEnd)|!(type=='DEL'&(InexactHOEnd-InexactHOStart>=6|nchar(Homology)<6)&endPosition-startPosition<1000)))
gridssCohortVariants = gridssCohortVariants #%>% 
  #filter(type!='DEL'|(InexactHOEnd-InexactHOStart<6&nchar(Homology)<6)|endPosition-startPosition>800|endPosition-startPosition<100) %>% 
  #filter(type!='INV'|endPosition-startPosition>40|nchar(Homology)<6)
View(gridssCohortVariants %>% filter(type=='INV',endPosition-startPosition<=40,nchar(Homology)>=6))
  #filter(type!='DEL'|nchar(Homology)<6|endPosition-startPosition>1000,type=='BND'|type=='INV'|endPosition-startPosition>100) %>%
  #filter((AdjCNChgEnd>0.15|AdjCNChgStart>0.15))

prodCohortVariants = prodCohortVariants %>% 
  filter(type!='DEL'|nchar(Homology)<6|endPosition-startPosition>1000) %>%
  filter((AdjCNChgEnd>0.15|AdjCNChgStart>0.15),DupBEStart!='true',DupBEEnd!='true',TransType!='SPAN')

######## START PROCESSING SINGLE PATIENT
annotatedGT1 = data.frame(stringsAsFactors = F)
annotatedGT2 = data.frame(stringsAsFactors = F)
annotatedMT1 = data.frame(stringsAsFactors = F)
annotatedMT2 = data.frame(stringsAsFactors = F)

library(cowplot)
#pdf('~/hmf/analyses/gridssAssessment/gridsVsManta.pdf',width=6,height=4,paper='special') 
for (i in 1:nrow(cohort)) {
  patient = cohort[i, ]
  GT1 = gridssCohortVariants %>% filter(sampleId == patient$Sample1)
  GT2 = gridssCohortVariants %>% filter(sampleId == patient$Sample2)
  MT1 = prodCohortVariants %>% filter(sampleId == patient$Sample1)
  MT2 = prodCohortVariants %>% filter(sampleId == patient$Sample2)
  combinedScope = combined_scope(GT1, GT2, MT1, MT2)
  GT1 = combinedScope[["GT1"]]
  GT2 = combinedScope[["GT2"]]
  MT1 = combinedScope[["MT1"]]
  MT2 = combinedScope[["MT2"]]

  signature = combined_signature(GT1, GT2, MT1, MT2)

  #sigPlot = plot_sv_signature(signature) + ggtitle(patient$patientId)
  ####pdf(file=paste('~/hmf/analyses/gridssAssessment/',patient$patientId,".pdf",sep=""),width=10)
  ####print(sigPlot)

  #save_plot(paste0("/Users/peterpriestley/hmf/analyses/gridssAssessment/", patient$patientId, ".png"), sigPlot, base_width = 18, base_height = 6)

  annotatedGT1 = bind_rows(GT1, annotatedGT1)
  annotatedGT2 = bind_rows(GT2, annotatedGT2)
  annotatedMT1 = bind_rows(MT1, annotatedMT1)
  annotatedMT2 = bind_rows(MT2, annotatedMT2)

}
#dev.off


createBuckets <- function(cluster) {
  cluster$PloidyBucket=2**(pmin(5,pmax(-3,round(log(cluster$Ploidy,2),0))))
  cluster$CnChEndBucket=2**(pmin(5,pmax(-3,round(log(pmax(cluster$AdjCNChgEnd,0.01),2),0))))
  cluster$CnChStartBucket=2**(pmin(5,pmax(-3,round(log(pmax(cluster$AdjCNChgStart,0.01),2),0))))
  cluster$ClusterCountBucket=2**(pmin(5,pmax(-3,round(log(cluster$ClusterCount,2),0))))
  cluster$LnkLenStartBucket=ifelse(!cluster$LnkLenStart>0,0,2**(pmin(20,pmax(0,round(log(cluster$LnkLenStart,2),0)))))
  cluster$FoldBackLenStartBucket=ifelse(!cluster$FoldbackLenStart>0,0,2**(pmin(20,pmax(0,round(log(cluster$FoldbackLenStart,2),0)))))
  cluster$LnkLenEndBucket=ifelse(!cluster$LnkLenEnd>0,0,2**(pmin(20,pmax(0,round(log(cluster$LnkLenEnd,2),0)))))
  cluster$FoldBackLenEndBucket=ifelse(!cluster$FoldbackLenEnd>0,0,2**(pmin(20,pmax(0,round(log(cluster$FoldbackLenEnd,2),0)))))
  cluster$LengthBucket=ifelse(cluster$Type=='BND'|cluster$Type=='INS'|cluster$endPosition-cluster$startPosition==0|cluster$ArmEnd!=cluster$ArmStart,
                              0,2**(pmin(25,pmax(0,round(log(cluster$endPosition-cluster$startPosition,2),0)))))
  cluster$HomLenBucket=2**(round(log(nchar(cluster$Homology),2),0))
  cluster$stressedArm=(cluster$ArmCountStart>1.3*cluster$ArmExpStart+6|cluster$ArmCountEnd>1.3*cluster$ArmExpEnd+6)
  cluster
}

annotatedMT1=createBuckets(annotatedMT1)
annotatedGT1=createBuckets(annotatedGT1)
annotatedMT2=createBuckets(annotatedMT2)
annotatedGT2=createBuckets(annotatedGT2)

#View(gridssCohortVariants %>% arrange(-Ploidy))
#####TESTING
View(annotatedMT2 %>% filter(sampleId=='CPCT02020351TII'))# %>% filter(DupBEStart=='false',DupBEEnd=='false',MantaPrecise=='true',is.na(GT1),type!='INS',Ploidy>0.3,AdjCNChgStart>0.3|AdjCNChgEnd>0.3))

View(annotatedGT2 %>% filter(type!='DEL'|nchar(Homology)<6|endPosition-startPosition>1000,type=='BND'|type=='INV'|endPosition-startPosition>100) %>% group_by(LengthBucket,type,GT2,MT1,MT2) %>% count() %>% spread(LengthBucket,n))

View(annotatedGT2 %>% 
       filter(type!='DEL'|nchar(Homology)<6|endPosition-startPosition>1000,type=='BND'|type=='INV'|endPosition-startPosition>100) %>% 
       group_by(ClusterCountBucket,sampleId) %>% count() %>% spread(ClusterCountBucket,n))
View(annotatedGT2%>% 
       filter(type!='DEL'|nchar(Homology)<6|endPosition-startPosition>1000,type=='BND'|type=='INV'|endPosition-startPosition>100))

View(annotatedGT2 %>% 
       filter(type!='DEL'|nchar(Homology)<6|endPosition-startPosition>1000,type=='BND'|type=='INV'|endPosition-startPosition>100) %>% 
       group_by(ClusterCount,ClusterDesc,stressedArm) %>% count() %>% spread(stressedArm,n)
       )

View(annotatedGT2 %>% filter(ClusterCount==1,stressedArm==F,type=='BND'|type=='INV'|type=='CRS'))

View(annotatedGT2 %>% filter(sampleId=='CPCT02080070TII',!is.na(GT1)|!is.na(MT1)|!is.na(MT2)) %>% mutate(len=endPosition-startPosition) %>% 
       select(Type,len,Ploidy,Homology,InsertSeq,everything()))

View(annotatedGT1 %>% filter(sampleId=='CPCT02010276T',is.na(GT2)) %>% mutate(len=endPosition-startPosition) %>% filter(Type=='INV') %>%
       select(Type,len,Ploidy,Homology,InsertSeq,MT1,MT2,everything()))

View(annotatedGT2 %>% filter(sampleId=='CPCT02140041TII') %>% mutate(len=endPosition-startPosition) %>% filter(is.na(GT1)) %>%
       select(Type,len,Ploidy,Homology,InsertSeq,GT1,MT1,MT2,everything()))

View(annotatedGT1 %>% mutate(len=endPosition-startPosition) %>% filter(startChromosome==11|endChromosome==11,sampleId=='CPCT02010276T') %>%
       select(Type,len,Ploidy,Homology,InsertSeq,GT2,MT1,MT2,everything()))


View(annotatedGT1 %>% mutate(len=endPosition-startPosition) %>% filter(type=='INV') %>%
       select(Type,len,Ploidy,Homology,InsertSeq,GT2,MT1,MT2,everything()))


View(annotatedGT2 %>% filter(sampleId=='CPCT02010641TII',is.na(MT2)) %>% mutate(len=endPosition-startPosition) %>%
       select(Type,len,Ploidy,Homology,InsertSeq,GT1,MT1,everything()))

View(annotatedGT2 %>% group_by(GT1,DupBEEnd,DupBEStart) %>% count() %>% spread(GT1,n))

########################################################
######### FULL COHORT ANALYSIS ############################
########################################################



#gridssCohortVariants = old_column_names(read.csv('~/Dropbox/HMF Australia team folder/Structural Variant Analysis/CLUSTER_GRIDSS.csv', header = T, stringsAsFactors = F))
gridssCohortVariants = old_column_names(read.csv('~/Dropbox/HMF Australia team folder/Structural Variant Analysis/CLUSTER_GRIDSS_484.csv', header = T, stringsAsFactors = F))
gridssCohortVariants=createBuckets(gridssCohortVariants)

View(gridssCohortVariants %>% filter(ResolvedType!='LowQual',AdjCNChgEnd<0.2|AdjCNChgStart<0.2))
View(gridssCohortVariants %>% filter(sampleId=='CPCT02220023T',ClusterId==586)%>% arrange(startChromosome,startPosition) 
     %>% select(ResolvedType,LengthBucket,Ploidy,type,ClusterCount,AsmbStart,AsmbEnd,everything()))# %>% group_by(type,startChromosome,ArmStart) %>% count() %>% spread(type,n))#select(LengthBucket,everything()))

############# THIS ONE #################
View(gridssCohortVariants %>% filter( grepl('CPCT02220023T',sampleId),ClusterCount>0,startChromosome=='3'|endChromosome=='3')%>% arrange(startChromosome,startPosition) 
     %>% select(ResolvedType,LengthBucket,Ploidy,type,ClusterCount,AsmbStart,AsmbEnd,everything()))# %>% group_by(ClusterCount,startChromosome,ArmStart,type) %>% count() %>% spread(type,n))#select(LengthBucket,everything()))

View(gridssCohortVariants %>% filter(ResolvedType!='LowQual',sampleId=='CPCT02010382TII'))#%>% group_by(ResolvedType,type) %>% count() %>% spread(type,n))

View(gridssCohortVariants %>% filter(ResolvedType!='LowQual',grepl('SglPair_DUP',ResolvedType))%>% select(ResolvedType,ClusterReason,Ploidy,AdjCNChgStart,AdjCNChgEnd,type,ClusterCount,everything())) #

View(gridssCohortVariants %>% filter(ResolvedType!='LowQual') %>% group_by(sampleId,type) %>% count() %>% spread(type,n))

View(gridssCohortVariants %>% filter(ResolvedType!='LowQual',ClusterCountBucket==1,type=='INV',ArmEnd==ArmStart))#  %>% group_by(sampleId,startChromosome) %>% count() %>% arrange(-n))


#####################
#### FOLDBACK RULES   #############
# adjacent BE with same orientation
# BE can form a continuous chain in same cluster OR simple INV
# copy number within 0.25 or 10%
# if multiple possible in a chain, use actual INV else shortest
View(gridssCohortVariants %>% filter(sampleId=='CPCT02010503T',!is.na(FoldbackLnkStart)|!is.na(FoldbackLnkEnd),AdjCNChgStart>0.5,AdjCNChgEnd>0.5) %>% select(ChainId,ChainIndex,ChainCount,FoldbackLnkStart,FoldbackLenStart,FoldbackLnkEnd,FoldbackLenEnd,LengthBucket,Ploidy,AdjCNChgStart,AdjCNChgEnd,startChromosome,endChromosome,startPosition,endPosition,startOrientation,endOrientation,everything()))
View(gridssCohortVariants %>% filter( grepl('DRUP01010024T',sampleId),ClusterCount>0,startChromosome=='8'|endChromosome=='8')%>% arrange(startChromosome,startPosition) 
     %>% select(ChainId,ChainIndex,ChainCount,FoldbackLnkStart,FoldbackLenStart,FoldbackLnkEnd,FoldbackLenEnd,LengthBucket,Ploidy,AdjCNChgStart,AdjCNChgEnd,startChromosome,endChromosome,startPosition,endPosition,startOrientation,endOrientation,everything()))
   #  %>% select(ResolvedType,LengthBucket,Ploidy,type,ClusterCount,AsmbStart,AsmbEnd,everything()))
     
View(gridssCohortVariants %>% filter(!is.na(FoldbackLnkStart)|!is.na(FoldbackLnkStart),AdjCNChgStart>0.8,AdjCNChgEnd>0.8) %>% group_by(sampleId) %>% count()) 

View(gridssCohortVariants %>% filter(!is.na(FoldbackLnkStart)|!is.na(FoldbackLnkStart),AdjCNChgStart>0.5,AdjCNChgEnd>0.5,sampleId =='CPCT02050120T') %>% group_by(sampleId,ClusterId,startChromosome,endChromosome) %>% count() %>% spread(startChromosome,n)) 


View(gridssCohortVariants %>% filter(type=="DUP"|type=='DEL',!is.na(FoldbackLnkStart)|!is.na(FoldbackLnkEnd)) %>% select(ChainId,ChainIndex,ChainCount,FoldbackLnkStart,FoldbackLenStart,FoldbackLnkEnd,FoldbackLenEnd,LengthBucket,Ploidy,AdjCNChgStart,AdjCNChgEnd,startChromosome,endChromosome,startPosition,endPosition,startOrientation,endOrientation,everything()))

View(gridssCohortVariants %>% filter(!is.na(FoldbackLnkStart)|!is.na(FoldbackLnkStart),AdjCNChgStart>0.5,AdjCNChgEnd>0.5) %>% group_by(FoldBackLenStartBucket,type) %>% count() %>% spread(type,n)) 

#########################
########## CHARTING ####################

chartSample = 'CPCT02040113T'; chartChr = '19';# chartChr2='22';
print(ggplot(data = gridssCohortVariants %>% filter(sampleId==chartSample,startChromosome==chartChr|endChromosome==chartChr) %>% mutate(
          modStart=ifelse(startChromosome==chartChr,startPosition,as.numeric(startChromosome)*-1e6),modEnd=ifelse(endChromosome==chartChr,endPosition,as.numeric(endChromosome)*-1e6)),aes(modStart,modEnd)) + 
        geom_point(aes(size = Ploidy,colour=ClusterDesc)) + theme_bw())
#View(gridssCohortVariants %>% filter( grepl(chartSample,sampleId),ClusterCount>31)%>% arrange(startChromosome,startPosition) 
#     %>% select(LengthBucket,AdjCNChgStart,AdjCNChgEnd,Ploidy,type,AsmbStart,AsmbEnd,LEStart,LEEnd,startChromosome,endChromosome,startPosition,endPosition,startOrientation,endOrientation,ArmStart,ArmEnd,ClusterId,SubClusterId))#,everything()))

View(gridssCohortVariants %>% filter(sampleId == chartSample,ClusterCount>8) %>% unite(chr_arm_start,startChromosome,ArmStart) %>% unite(chr_arm_end,endChromosome,ArmEnd) %>%
       group_by(ClusterId,sampleId,chr_arm_start,chr_arm_end) %>% count() %>% spread(chr_arm_end,n,fill=""))

View(gridssCohortVariants %>% filter(sampleId==chartSample,startChromosome==chartChr|endChromosome==chartChr))#

  View(gridssCohortVariants %>% filter(ClusterCount>0,Ploidy>0.5) %>% group_by(ResolvedType,type) %>% count() %>% spread(type,n))
View(gridssCohortVariants %>% filter(sampleId==chartSample,ClusterCount<=8,ClusterCount>0))#%>% group_by(ClusterCountBucket,CnChStartBucket) %>% count() %>% spread(CnChStartBucket,n))
View(temp)

View(gridssCohortVariants %>% filter(grepl('CCCCCCCCCCCC',InsertSeq)|grepl('GGGGGGGGGGGG',InsertSeq),ResolvedType!='LowQual'))
View(gridssCohortVariants %>% filter(grepl('GGGGGGGGGGGG',InsertSeq)|grepl('CCCCCCCCCCCC',InsertSeq),ResolvedType!='LowQual') %>% group_by(round(AdjCNStart,-1),Type) %>% count() %>% spread(Type,n))
#########################
#############################
View(rbind(gridssCohortVariants %>% unite(chr_arm,startChromosome,ArmStart) %>% filter(ResolvedType!='LowQual',type!='NONE') %>%
             group_by(sampleId,chr_arm,CC=pmin(ClusterCount,4),id=ClusterId) %>% count(),gridssCohortVariants %>%  unite(chr_arm,endChromosome,ArmEnd) %>% filter(ResolvedType!='LowQual',type!='NONE') %>%
             group_by(sampleId,chr_arm,CC=pmin(ClusterCount,4),id=ClusterId) %>% count()) %>% 
       group_by(sampleId,chr_arm,CC,id) %>% count() %>%
       group_by(sampleId,CC,chr_arm) %>% count() %>% spread(CC,n) %>%
       filter(chr_arm!='0_P'))

View(rbind(gridssCohortVariants %>% unite(chr_arm,startChromosome,ArmStart) %>%
             group_by(sampleId,ResolvedType,id=ClusterId) %>% count(),gridssCohortVariants %>%  unite(chr_arm,endChromosome,ArmEnd) %>%
             group_by(sampleId,ResolvedType,id=ClusterId) %>% count()) %>% 
       group_by(sampleId,ResolvedType,id) %>% count() %>%
       group_by(sampleId,ResolvedType) %>% count() %>% spread(ResolvedType,n))

View(rbind(gridssCohortVariants %>% unite(chr_arm,startChromosome,ArmStart) %>%
             group_by(sampleId,ResolvedType,ClusterCount,chr_arm,id=ClusterId) %>% count(),gridssCohortVariants %>%  unite(chr_arm,endChromosome,ArmEnd) %>%
             group_by(sampleId,ResolvedType,ClusterCount,chr_arm,id=ClusterId) %>% count()) %>% 
       filter(ResolvedType=='Line',ClusterCount>8) %>%
       group_by(sampleId,id,chr_arm) %>% summarise(n=sum(n)/2) %>%
       group_by(sampleId,id) %>% summarise(arms=n(),svCount=sum(n)) )

View(rbind(gridssCohortVariants %>% unite(chr_arm,startChromosome,ArmStart) %>%
             group_by(sampleId,chr_arm,id=ClusterId) %>% count(),gridssCohortVariants %>%  unite(chr_arm,endChromosome,ArmEnd) %>%
             group_by(sampleId,chr_arm,id=ClusterId) %>% count()) %>% 
       group_by(sampleId,chr_arm,id) %>% count() %>%
       group_by(sampleId,chr_arm) %>% count()  %>%
       filter(chr_arm!='0_P') %>% spread(chr_arm,n))

View(gridssCohortVariants %>% filter(ClusterCount>5) %>% group_by(sampleId,ClusterId) %>% count() %>% group_by(sampleId) %>% summarise(max(n),n()))

View(gridssCohortVariants %>% filter(ClusterCount>0) %>% 
       group_by(sampleId,startChromosome,type) %>% count() %>% spread(type,n) )


######################

View(gridssCohortVariants %>% filter(grepl('AAAAAAAAA',InsertSeq)|grepl('TTTTTTTTT',InsertSeq),type=='DUP'))
View(gridssCohortVariants %>% group_by(type,grepl('bebe',AsmbStart)|grepl('bebe',AsmbEnd)) %>% count() %>% spread(type,n))
View(gridssCohortVariants %>% filter(type=='BND') %>% group_by(sampleId,startChromosome,endChromosome,ArmStart,ArmEnd,isLine=LEStart!='false'|LEEnd!='false') %>% count() %>% spread(isLine,n))

View(gridssCohortVariants %>% filter(type=='BND') %>% group_by(sampleId,startChromosome,endChromosome,ArmStart,ArmEnd,isLine=LEStart!='false'|LEEnd!='false') %>% count() %>% group_by(sampleId) %>% count())

View(gridssCohortVariants %>% filter(ClusterCount>2) %>% unite(chr_arm,startChromosome,ArmStart) %>% group_by(ClusterId,sampleId,chr_arm) %>% count() %>% group_by(ClusterId,sampleId) %>% count())
View(gridssCohortVariants %>% filter(startChromosome!=endChromosome|ArmStart!=ArmEnd,endChromosome!=0,ResolvedType!='LowQual') %>% group_by(sampleId,startChromosome,endChromosome,ArmStart,ArmEnd) %>% count())# %>% group_by(sampleId) %>% count())

View(gridssCohortVariants %>% filter(startChromosome!=endChromosome|ArmStart!=ArmEnd,endChromosome!=0,ResolvedType!='LowQual') %>% group_by(sampleId,startChromosome,endChromosome,ArmStart,ArmEnd,isLine=LEStart!='false'|LEEnd!='false') %>% count() %>% spread(isLine,n))

temp = (gridssCohortVariants %>% filter(startChromosome!=endChromosome|ArmStart!=ArmEnd) %>% group_by(sampleId,CNSupport=(AdjCNChgEnd>0.15|AdjCNChgStart>0.15)) %>% count() %>% spread(CNSupport,n))

temp = (gridssCohortVariants %>% filter(startChromosome!=endChromosome|ArmStart!=ArmEnd,ResolvedType!='LowQual') %>% 
          group_by(sampleId,startChromosome,endChromosome,isLine=LEStart!='None'|LEEnd!='None') %>% count() %>% group_by(sampleId,isLine) %>% count() %>% spread(isLine,nn))
View(temp)
