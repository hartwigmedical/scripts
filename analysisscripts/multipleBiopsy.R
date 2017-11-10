 ### ssh -L 3307:localhost:3306 peter@ext-hmf-datastore
 library(RMySQL)
 library(data.table)
 library(ggplot2)
 library(grid)

 multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
   require(grid)
   
   # Make a list from the ... arguments and plotlist
   plots <- c(list(...), plotlist)
   
   numPlots = length(plots)
   
   # If layout is NULL, then use 'cols' to determine layout
   if (is.null(layout)) {
     # Make the panel
     # ncol: Number of columns of plots
     # nrow: Number of rows needed, calculated from # of cols
     layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                      ncol = cols, nrow = ceiling(numPlots/cols))
   }
   
   if (numPlots==1) {
     print(plots[[1]])
     
   } else {
     # Set up the page
     grid.newpage()
     pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
     
     # Make each plot, in the correct location
     for (i in 1:numPlots) {
       # Get the i,j matrix positions of the regions that contain this subplot
       matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
       
       print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                       layout.pos.col = matchidx$col))
     }
   }
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
 
 dbConnect = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
 multipleBiopsies<-get_multiple_biopsy_samples(dbConnect)
 for (i in (1:nrow(multipleBiopsies))) {
  variants<-get_multiple_biopsy_somatic_variants(dbConnect,multipleBiopsies[i,])
  chart_multiple_biopsies(variants,multipleBiopsies[i,])
 }

 
 

         