 ### ssh -L 3307:localhost:3306 peter@ext-hmf-datastore
 library(RMySQL)
 library(data.table)
 library(ggplot2)
 library(dplyr)
 library(grid)
 library(gridExtra)
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

 get_gene_copy_number<-function(dbConnect,sampleString="",minCopyNumber=-1e6,maxCopyNumber=1e6)
 {
   query = paste("SELECT sampleId,gene,minCopyNumber,regions FROM gene_copy_number ",
                 "WHERE minCopyNumber >",minCopyNumber," ",
                 "AND minCopyNumber <",minCopyNumber," ",
                 sep="")

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

 
 get_multiple_biopsy_structural_variants<-function(dbConnect, multipleBiopsy)
 {
   query = paste(
     "select startChromosome,startPosition,endChromosome,endPosition,type,count(*) as sampleCount,",
     "round(sum(if(sampleId='",multipleBiopsy$sampleId1,"',ploidy,0)),2) as Sample1Ploidy, ",
     "round(sum(if(sampleId='",multipleBiopsy$sampleId2,"',ploidy,0)),2) as Sample2Ploidy ",
     "from hmfpatients.structuralVariant where sampleId in (",
     paste(c(shQuote(multipleBiopsy$sampleId1),shQuote(multipleBiopsy$sampleId2)),collapse = ","),
     ") group by 1,2,3,4,5",
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
     "select c.*, sum(r.observedTumorRatioCount) as observedTumorRatioCount, ",
     "sum(gcContent*r.observedTumorRatioCount)/sum(r.observedTumorRatioCount) as GCContent ",
     "from copyNumber c, copyNumberRegion r ",
     "where c.chromosome = r.chromosome and c.start <= r.start and c.end >= r.end ",
     "and c.sampleId = r.sampleId and c.sampleId = '",sampleId,"' group by c.id;",
      sep = "")
   raw_data = dbGetQuery(dbConnect, query)
   chrOrder <-c((1:22),"X","Y")
   raw_data$chromosome <- factor(raw_data$chromosome, chrOrder, ordered=TRUE)
   raw_data
 }

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
   p3<-ggplot()+labs(title=c(multipleBiopsy$SampleId1,multipleBiopsy$SampleId2))+geom_point(data=variants,aes(Sample1Ploidy,Sample2Ploidy))+xlim(0,3)+ylim(0,3)
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
   p1<-ggplot(sampleData,aes(purityBucket)) + geom_bar(aes(fill=cancerType))
   p2<-ggplot(sampleData,aes(ploidyBucket)) + geom_bar(aes(fill=cancerType))
   multiplot(p1,p2)
 }
 ######### SINGLE BIOPSY SV PLOIDIES ##########
 sampleId = "CPCT02380019T"
 dbConnect = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
 #somatics
 variants<-get_somatic_variants(dbConnect,sampleId)
 ggplot() + xlim(0,3.5)+ labs(title = sampleId,x='') +
   geom_histogram(aes(x=ploidy),data=variants,fill = "red", alpha = 0.6,binwidth = 0.02)
 ggplot() + xlim(0,3.5)+ labs(title = sampleId,x='') +
   geom_histogram(aes(x=ploidy),data=variants,fill = "red", alpha = 0.6,binwidth = 0.02)+
   facet_wrap( ~chromosome )

 #CNV
 PURITY<-as.double(get_sample_purity(dbConnect,sampleId))
 CNV<-as.data.table(get_copy_number_data(dbConnect,sampleId))
 CNV[, prevCopyNumber:=c(NA, copyNumber[-.N]), by=chromosome]
 CNV[, prevGCContent:=c(NA, GCContent[-.N]), by=chromosome]
 CNV$inferred<-(CNV$copyNumberMethod=="STRUCTURAL_VARIANT")
 CNV[, prevInferred:=c(NA, inferred[-.N]), by=chromosome]
 CNV$chCopyNumber<-CNV$copyNumber-CNV$prevCopyNumber
 CNV$Length<-CNV$end-CNV$start
 CNV$chGCContent<-CNV$GCContent-CNV$prevGCContent
 CNV<-transform(CNV, percentChCopyNumber = chCopyNumber/pmax(copyNumber, prevCopyNumber))
 CNV[, prevLength:=c(NA, Length[-.N]), by=chromosome]
 CNV<-transform(CNV, minLength = pmin(Length, prevLength))

 # SV PLOIDY vs CopyNumber Change
 SV<-get_SV_ends(dbConnect,sampleId)
 SV<-merge(x = SV, y = CNV, by.x = c("chromosome","position"),by.y=c("chromosome","start"), all.x = TRUE)
 SV$ploidy<-ifelse(SV$orientation ==1, -(SV$prevCopyNumber*PURITY+(1-PURITY)*2)*SV$orientation*SV$AF/PURITY ,
                   -(SV$copyNumber*PURITY+(1-PURITY)*2)*SV$orientation*SV$AF/PURITY)
 ggplot() + xlim(0,3.5)+ labs(title = sampleId,x='') +
   geom_histogram(aes(x=ploidy),data=SV[(SV$inferred==FALSE & SV$prevInferred==FALSE),],fill = "red", alpha = 0.8,binwidth = 0.1) +
   geom_histogram(aes(x=ploidy),data=SV[(SV$inferred==TRUE | SV$prevInferred==TRUE),],fill = "blue", alpha = 0.2,binwidth = 0.1)
 ggplot(SV[(SV$inferred==FALSE & SV$prevInferred==FALSE),],aes(ploidy,chCopyNumber))+geom_point()+
   labs(title="SV Ploidy vs ChCopyNumber Scatter NORMAL")+ xlim(-4,4)+ylim(-4,4)+
   facet_wrap( ~segmentStartSupport )
 ggplot(SV[(SV$inferred==TRUE | SV$prevInferred==TRUE)&SV$minLength>300,],aes(ploidy,chCopyNumber))+geom_point()+
   labs(title="SV Ploidy vs ChCopyNumber Scatter INFERRED")+ xlim(-4,4)+ylim(-4,4)+
   facet_wrap( ~segmentStartSupport )
 dbDisconnect(dbConnect)


 # LENGTH PLOTTING

 ggplot(aes(chCopyNumber),data=CNV) + stat_ecdf(geom = "step", pad = FALSE) +
   facet_wrap( ~segmentStartSupport ) + xlim(-4,4) + labs(title=sampleId)
 ggplot(aes(percentChCopyNumber),data=CNV) + stat_ecdf(geom = "step", pad = FALSE) +
   facet_wrap( ~segmentStartSupport ) + xlim(-1,1) + labs(title=sampleId)
 ggplot(aes(minLength),data=CNV) + stat_ecdf(geom = "step", pad = FALSE) +
   scale_x_log10() + facet_wrap( ~segmentStartSupport )+ labs(title=sampleId)
 CNV[, .(count=.N), by = segmentStartSupport]
 PURITY
 ggplot(CNV[segmentStartSupport=="NONE",],aes(minLength,chGCContent))+geom_point()+labs(title="Length vs GC Content Scatter") + scale_x_log10() +
   facet_wrap( ~segmentStartSupport )
 ggplot(CNV,aes(minLength,percentChCopyNumber))+geom_point()+labs(title="Length vs GC Content Scatter") + scale_x_log10() +
   facet_wrap( ~segmentStartSupport )




 ######### MULTIPLE BIOPSIES LOGIC  ##########

 dbConnect = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
 multipleBiopsies<-get_multiple_biopsy_samples(dbConnect)
 multipleBiopsies<-data.frame(sampleId1='CPCT02020506T',sampleId2='CPCT02020506TII')
 for (i in (1:nrow(multipleBiopsies))) {
   #variants<-get_multiple_biopsy_somatic_variants(dbConnect,multipleBiopsies[i,])
   #chart_multiple_biopsies(variants,multipleBiopsies[i,])
   variants<-get_multiple_biopsy_structural_variants(dbConnect,multipleBiopsies[i,])
   chart_multiple_biopsies(variants,multipleBiopsies[i,])
 }

 dbDisconnect(dbConnect)
 multipleBiopsy=multipleBiopsies[1,]
 head(variants[variants$sampleCount==2 & variants$Sample1Ploidy < 0.5 & variants$Sample2Ploidy > 0.7,],100)
 ggplot()+labs(title=c(multipleBiopsy$SampleId1,multipleBiopsy$SampleId2))+geom_point(data=variants,aes(Sample1Ploidy,Sample2Ploidy))+ facet_wrap( ~chromosome)+xlim(0,3)+ylim(0,3)

 ########### SAMPLE DATA + PURITY & PLOIDY ###############

 dbConnect = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
 sampleString<-paste(shQuote(get_QCpass_samples(dbConnect)$sampleId),collapse=",")
 sampleData<-get_sample_data(dbConnect,sampleString)

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
 dbConnect = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
 geneString<-paste(shQuote(get_HMF_panel(dbConnect)$gene),collapse=",")
 dbDisconnect(dbConnect)
 dbConnect = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis")
 sampleString<-paste(shQuote(get_QCpass_samples(dbConnect)$sampleId),collapse=",")
 dbDisconnect(dbConnect)

 #SOMATICS
 dbConnect = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis")
 somatics<-get_panel_somatics(dbConnect,geneString,sampleString)
 dbDisconnect(dbConnect)

 #CNV DELETIONS


 #disruptions - CURRENTLY NOT WORKING
 dbConnect = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
 disruptions<-get_panel_disruptions(dbConnect,geneString,sampleString)
 dbDisconnect(dbConnect)

 #################################
 # SCRATCH
 CNV %>% group_by(segmentStartSupport,segmentEndSupport) %>% summarise(sumCount=sum(observedTumorRatioCount),medCount=median(observedTumorRatioCount),countChromosome=n())
 CNV %>% group_by(chromosome) %>% summarise(sumCount=sum(observedTumorRatioCount),count=n())

 df = data.frame(x = c(1, 14, 3, 21, 11), y = c(102, 500, 40, 101, 189))
 apply(df, 2, function(x) x - c(PURITY))

 library(plyr)
 d.f <- data.frame(
   grp = as.factor( rep( c("A","B"), each=120 ) ) ,
   grp2 = as.factor( rep( c("cat","dog","elephant"), 40 ) ) ,
   val = c( sample(c(2:4,6:8,12),120,replace=TRUE), sample(1:4,120,replace=TRUE) )
 )
 d.f <- arrange(d.f,grp,grp2,val)
 d.f.ecdf <- ddply(d.f, .(grp,grp2), transform, ecdf=ecdf(val)(val) )

 ggplot( d.f.ecdf, aes(val, ecdf, colour = grp) ) + geom_step() + facet_wrap( ~grp2 )

 d.f.ecdf



 cluster_s
