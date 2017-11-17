 ### ssh -L 3307:localhost:3306 peter@ext-hmf-datastore
 library(RMySQL)
 library(data.table)
 library(ggplot2)
 library(dplyr)
 library(grid)
 library(devtools)
 library("NMF")
 library(tidyverse)

 ############## GENE PANEL & MAPPINGS ###############
 
 PATH="hmf/repos/scripts/analysisscripts/"
 genePanel<-read.table(paste(PATH,"genePanel.txt",sep=""), sep="\t",col.names="gene")
 tumorMapping<-read.table(paste(PATH,"tumorMappings.txt",sep=""), sep="\t",col.names=c("cancerType","category"))
 geneString<-paste(shQuote(genePanel$gene),collapse=",")
 

 
 get_HMF_panel<-function(dbConnect)
 {
   query = "SELECT distinct gene from geneCopyNumber"
   print(query)
   raw_data = dbGetQuery(dbConnect, query)
 }
 
 get_panel_somatics<-function(dbConnect,geneString,sampleString="")
 {
   query = paste("SELECT gene,sampleId,ref,alt,cosmicId,dbsnpId,effect,alleleReadCount,totalReadCount	",
            ",adjustedVaf,adjustedCopyNumber , highConfidence,trinucleotideContext,	clonality, loh ",
            "FROM somaticVariant ",
            "WHERE gene in (",geneString,") ",
            "AND effect NOT IN ('UTR variant','intron variant','sequence feature','synonymous variant') ",
            "AND filter = 'PASS' ",
            sep="")
   
   if (sampleString != ""){
      query=paste(query,"AND sampleId in (",sampleString,")" )
   }

   raw_data = dbGetQuery(dbConnect, query)
 }

 get_gene_copy_number<-function(dbConnect,sampleString="")
 {
   query = "SELECT gene,minCopyNumber,regions FROM gene_copy_number"
   
   if (sampleString != ""){
     query=paste(query,"AND sampleId in (",sampleString,")" )
   }
   raw_data = dbGetQuery(dbConnect, query)
 }
 
 
 get_panel_disruptions<-function(dbConnect,geneString,sampleString="")
 {
   query = paste("SELECT * FROM gene_disruption_view v ",
                 "WHERE gene in (",geneString,") ",
                 sep="")
   
   if (sampleString != ""){
     query=paste(query,"AND sampleId in (",sampleString,")" )
   }
   raw_data = dbGetQuery(dbConnect, query)
 }
 
  
 get_QCpass_samples<-function(dbConnect)
 {
   query = "SELECT sampleId from purity where qcstatus = 'PASS' and status <> 'NO_TUMOR'"
   raw_data = dbGetQuery(dbConnect, query)
 }
 
 get_sample_purity<-function(dbConnect,sampleId)
 {
   query = paste("SELECT purity FROM purity where sampleId = '",sampleId,"'",sep="")
   raw_data = dbGetQuery(dbConnect, query)
 }
 
 get_sample_data<-function(dbConnect,sampleString="")
 {
   query = paste("SELECT purity.*,hospital,clinical.gender as clinicalGender,",
          "cancerType,biopsyLocation,treatment,biopsyDate,registrationDate,sampleArrivalDate ",
          "FROM purity left join clinical on purity.sampleId = clinical.sampleId ")
   
   if (sampleString != ""){
      query=paste(query,"WHERE purity.sampleId in (",sampleString,")" )
   }
   raw_data = dbGetQuery(dbConnect, query)
 }
 
 
 
 get_multiple_biopsy_samples<-function(dbConnect)
 {
   query = paste(
     "select left(p.sampleId,12) as patient,cancerType,count(*), min(p.sampleId) as sampleId1,max(p.sampleId) as sampleId2,",
     "round(min(purity),2) as minPurity,round(max(purity),2) as maxPurity,",
     "round(min(ploidy),2) as minPloidy,round(max(ploidy),2) as maxPloidy ",
     "from purity p left join clinical c on p.sampleId = c.sampleId ",
     "group by 1,2 having count(*) > 1",
     sep = "")
   #print(query)
   raw_data = dbGetQuery(dbConnect, query)
 }

 get_somatic_variants<-function(dbConnect, sampleId)
 {
   query = paste("select chromosome,position,ref,alt,length(ref)=length(alt) as isSNV,round(greatest(0,adjustedVaf*adjustedCopyNumber),2) as ploidy ",
     "from hmfpatients.somaticVariant where filter = 'PASS' and sampleId = '",
     sampleId,
     "'",
     sep = "")
   #print(query)
   raw_data = dbGetQuery(dbConnect, query)
 }
  
 get_multiple_biopsy_somatic_variants<-function(dbConnect, multipleBiopsy)
 {
   query = paste(
     "select chromosome,position,ref,alt,length(ref)=length(alt) as isSNV,count(*) as sampleCount,",
     "round(sum(if(sampleId='",multipleBiopsy$sampleId1,"',(greatest(0,adjustedVaf*adjustedCopyNumber)),0)),2) as Sample1Ploidy, ",
     "round(sum(if(sampleId='",multipleBiopsy$sampleId2,"',(greatest(0,adjustedVaf*adjustedCopyNumber)),0)),2) as Sample2Ploidy ",
     "from hmfpatients.somaticVariant where filter = 'PASS' and sampleId in (",
     paste(c(shQuote(multipleBiopsy$sampleId1),shQuote(multipleBiopsy$sampleId2)),collapse = ","),
     ") group by 1,2,3,4,5",
     sep = "")
   #print(query)
   raw_data = dbGetQuery(dbConnect, query)
 }
 
 get_copy_number_data<-function(dbConnect, sampleId)
 {
   query = paste(
     "select c.*, sum(r.observedTumorRatioCount) as observedTumorRatioCount ",
     "from copyNumber c, copyNumberRegion r ",
     "where c.chromosome = r.chromosome and c.start <= r.start and c.end >= r.end ",
     "and c.sampleId = r.sampleId and c.sampleId = '",sampleId,"' group by c.id;",
      sep = "")
   raw_data = dbGetQuery(dbConnect, query)
   chrOrder <-c((1:22),"X","Y")
   raw_data$chromosome <- factor(raw_data$chromosome, chrOrder, ordered=TRUE)
   raw_data
 }
 
 head(CNV)
 
 get_SV_ends<-function(dbConnect, sampleId)
 {
   query = paste("select startChromosome as chromosome,startPosition as position,startOrientation as orientation,startAF as AF, type ",
                 "from structuralVariant where sampleId = '",sampleId,"' ", 
                 "union ",
                 "select endChromosome as chromosome,endPosition as position,endOrientation as orientation,endAF as AF, type ",
                 " from structuralVariant where sampleId = '",sampleId,"' " ,
                 sep = "")
   raw_data = dbGetQuery(dbConnect, query)
 }

 ############ PLOTTING ###############
 multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
   require(grid)
   
   plots <- c(list(...), plotlist)
   numPlots = length(plots)
   
   if (is.null(layout)) {
     layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),ncol = cols, nrow = ceiling(numPlots/cols))
   }
   
   if (numPlots==1) {
     print(plots[[1]])
   } else {
     grid.newpage()
     pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
     for (i in 1:numPlots) {
       matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
       print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                       layout.pos.col = matchidx$col))
     }
   }
 }
 
  chart_multiple_biopsies<-function(variants,multipleBiopsy)
 {
   binsize<-0.05
   p1<-ggplot() + xlim(0,3.5)+ labs(title = multipleBiopsy$sampleId1,x='') + 
     geom_histogram(aes(x=Sample1Ploidy),data=subset(variants,sampleCount == 2),fill = "red", alpha = 0.6,binwidth = binsize) + 
     geom_histogram(aes(x=Sample1Ploidy),data=variants[variants$sampleCount==1 & variants$Sample1Ploidy>0,],fill = "red", alpha = 0.2,binwidth = binsize)
   p2<-ggplot() + xlim(0,3.5)+ labs(title = multipleBiopsy$sampleId2,x='ploidy')  + 
     geom_histogram(aes(x=Sample2Ploidy),data=subset(variants,sampleCount == 2),fill = "blue", alpha = 0.6,binwidth = binsize) + 
     geom_histogram(aes(x=Sample2Ploidy),data=variants[variants$sampleCount==1 & variants$Sample2Ploidy>0,],fill = "blue", alpha = 0.2,binwidth = binsize)
   p3<-ggplot(variants,aes(Sample1Ploidy,Sample2Ploidy))+geom_point() +xlim(0,3)+ylim(0,3)+labs(x=multipleBiopsy$SampleId1,y=multipleBiopsy$SampleId2)
   multiplot(p1,p2,p3,cols=2)
 }

 chart_CNV_by_chromosome<-function(CNV)
 {
   p <- ggplot(CNV, aes(start, copyNumber, group=chromosome)) + geom_step(aes(),shape=20, size=0.2, alpha=0.5) + ylim(0,6) +  facet_wrap(~chromosome, nrow=3, scales="free_x")
   p <- p + theme(axis.line.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
                  panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),panel.grid.major.y=element_line(colour="grey80", linetype=2, size=0.2),
                  panel.background=element_blank(),
                  panel.border=element_blank(),
                  plot.background=element_blank())
   print(p)
 }
 
 chart_purity_and_ploidy<-function(sampleData)
 {
   sampleData$purityBucket = round(sampleData$purity,1)
   sampleData$ploidyBucket = round(sampleData$ploidy,1)
   p1<-ggplot(sampleData,aes(purityBucket)) + geom_bar(aes(fill=category))
   p2<-ggplot(sampleData,aes(ploidyBucket)) + geom_bar(aes(fill=category))
   multiplot(p1,p2)
 }
 ######### SINGLE BIOPSY LOGIC  ########## 
 dbConnect = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis")
 sampleId = "DRUP01090005T"
 variants<-get_somatic_variants(dbConnect,sampleId)
  ggplot() + xlim(0,3.5)+ labs(title = sampleId,x='') + 
   geom_histogram(aes(x=ploidy),data=variants,fill = "red", alpha = 0.6,binwidth = 0.1) 

 ######### SINGLE BIOPSY SV PLOIDIES ########## 
 dbConnect = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis")
 sampleId = "DRUP01090005T"
 PURITY<-as.double(get_sample_purity(dbConnect,sampleId))
 CNV<-as.data.table(get_copy_number_data(dbConnect,sampleId))
 CNV[, prevCopyNumber:=c(NA, copyNumber[-.N]), by=chromosome]
 CNV[, prevInferred:=c(NA, inferred[-.N]), by=chromosome]
 
 SV<-get_SV_ends(dbConnect,sampleId)
 SV<-merge(x = SV, y = CNV, by.x = c("chromosome","position"),by.y=c("chromosome","start"), all.x = TRUE)
 SV$ploidy<-ifelse(SV$orientation ==1, -(SV$prevCopyNumber*PURITY+(1-PURITY)*2)*SV$orientation*SV$AF/PURITY ,
             -(SV$copyNumber*PURITY+(1-PURITY)*2)*SV$orientation*SV$AF/PURITY)#*PURITY#+(1-PURITY)*2)*SV$orientation*SV$AF/PURITY
 SV$chCopyNumber<-SV$copyNumber-SV$prevCopyNumber
 ggplot() + xlim(0,3.5)+ labs(title = sampleId,x='') + 
   geom_histogram(aes(x=ploidy),data=SV,fill = "red", alpha = 0.6,binwidth = 0.1) 
 ggplot(SV[SV$inferred==0 & SV$prevInferred==0,],aes(ploidy,chCopyNumber))+geom_point()+labs(title="SV Ploidy vs ChCopyNumber Scatter")+ xlim(-4,4)+ylim(-4,4)
 dbDisconnect(dbConnect)
 
 ######### MULTIPLE BIOPSIES LOGIC  ########## 
 
 dbConnect = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
 multipleBiopsies<-get_multiple_biopsy_samples(dbConnect)
 for (i in (1:nrow(multipleBiopsies))) {
   variants<-get_multiple_biopsy_somatic_variants(dbConnect,multipleBiopsies[i,])
   chart_multiple_biopsies(variants,multipleBiopsies[i,])
 }
 dbDisconnect(dbConnect)
 
 ########### SAMPLE DATA + PURITY & PLOIDY ###############

 dbConnect = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
 sampleString<-paste(shQuote(get_QCpass_samples(dbConnect)$sampleId),collapse=",")
 sampleData<-get_sample_data(dbConnect,sampleString)
 sampleData<-merge(x = sampleData, y = tumorMapping, by = "cancerType", all.x = TRUE)
 chart_purity_and_ploidy(sampleData)
 ggplot(sampleData,aes(purity,ploidy))+geom_point()+labs(title="Ploidy vs Purity Scatter")
 
 #sampleData %>% group_by(qcStatus,status) %>% summarise(count=n())
 #sampleData[sampleData$gender!=toupper(sampleData$clinicalGender) & !is.na(sampleData$clinicalGender), ]
 
 ######### CNV LOGIC  ##########
 
 dbConnect = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis")
 CNV<-get_copy_number_data(dbConnect,"CPCT02010288T")
 dbDisconnect(dbConnect)
 chart_CNV_by_chromosome(CNV)
 
 ########## DRIVERS ############
 dbConnect = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis")
 
 sampleString<-paste(shQuote(get_QCpass_samples(dbConnect)$sampleId),collapse=",")
 geneString<-paste(shQuote(get_HMF_panel(dbConnect)$gene),collapse=",")
 
 #SOMATICS
 somatics<-get_panel_somatics(dbConnect,geneString,sampleString)
 dbDisconnect(dbConnect)
 
 #disruptions
 dbConnect = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
 disruptions<-get_panel_disruptions(dbConnect,geneString,sampleString)
 dbDisconnect(dbConnect)
 
 #################################
 # SCRATCH
 CNV %>% group_by(segmentStartSupport,segmentEndSupport) %>% summarise(sumCount=sum(observedTumorRatioCount),medCount=median(observedTumorRatioCount),countChromosome=n())
 CNV %>% group_by(chromosome) %>% summarise(sumCount=sum(observedTumorRatioCount),count=n())
 
 df = data.frame(x = c(1, 14, 3, 21, 11), y = c(102, 500, 40, 101, 189)) 
 apply(df, 2, function(x) x - c(PURITY)) 
