 ### ssh -L 3307:localhost:3306 peter@ext-hmf-datastore
 library(RMySQL)
 library(data.table)
 library(ggplot2)
 library(dplyr)
 library(grid)
 library(devtools)
 library("NMF")

 ############## QUERIES
 
 get_sample_data<-function(dbConnect)
 {
   query = paste("SELECT purity.*,hospital,clinical.gender as clinicalGender,",
          "primaryTumorLocation,biopsyLocation,treatment,biopsyDate,registrationDate,sampleArrivalDate   ",
          "from purity left join clinical on purity.sampleId = clinical.sampleId;")
   raw_data = dbGetQuery(dbConnect, query)
 }
 
 get_multiple_biopsy_samples<-function(dbConnect)
 {
   query = paste(
     "select left(p.sampleId,12) as patient,primaryTumorLocation,count(*), min(p.sampleId) as sampleId1,max(p.sampleId) as sampleId2,",
     "round(min(purity),2) as minPurity,round(max(purity),2) as maxPurity,",
     "round(min(ploidy),2) as minPloidy,round(max(ploidy),2) as maxPloidy ",
     "from purity p left join clinical c on p.sampleId = c.sampleId ",
     "group by 1,2 having count(*) > 1",
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
   print(query)
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
 
 ######### MULTIPLE BIOPSIES LOGIC  ########## 
 dbConnect = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis")
 multipleBiopsies<-get_multiple_biopsy_samples(dbConnect)
 for (i in (1:nrow(multipleBiopsies))) {
   variants<-get_multiple_biopsy_somatic_variants(dbConnect,multipleBiopsies[i,])
   chart_multiple_biopsies(variants,multipleBiopsies[i,])
 }
 dbDisconnect(dbConnect)
 
 ######## PURITY & PLOIDY ###############
 dbConnect = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
 sampleData<-get_sample_data(dbConnect)
 dbDisconnect(dbConnect) 
 sampleData %>% group_by(qcStatus,status) %>% summarise(count=n())
 out<-sampleData %>% filter(qcStatus=='PASS',status!='NO_TUMOR') %>% group_by(primaryTumorLocation) %>% summarise(sampleCount=n())
 out<-out[order(-out$sampleCount),]
 out$primaryTumorLocation <-factor(out$primaryTumorLocation, levels =  out$primaryTumorLocation[order(-out$sampleCount)])
 sampleData[sampleData$gender!=toupper(sampleData$clinicalGender) & !is.na(sampleData$clinicalGender), ]
 ggplot(out,aes(x=primaryTumorLocation,y=sampleCount)) + geom_bar(stat="identity") + theme_classic()
 
 ######### CNV LOGIC  ##########
 dbConnect = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis")
 CNV<-get_copy_number_data(dbConnect,"CPCT02010288T")
 dbDisconnect(dbConnect)
 chrOrder <-c((1:22),"X","Y")
 CNV$chromosome <- factor(CNV$chromosome, chrOrder, ordered=TRUE)
 chart_CNV_by_chromosome(CNV)

 # SCRATCH
 CNV %>% group_by(segmentStartSupport,segmentEndSupport) %>% summarise(sumCount=sum(observedTumorRatioCount),medCount=median(observedTumorRatioCount),countChromosome=n())
 CNV %>% group_by(chromosome) %>% summarise(sumCount=sum(observedTumorRatioCount),count=n())
 