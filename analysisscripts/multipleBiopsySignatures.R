detach("package:purple", unload=TRUE)
library(purple)
library(RMySQL)
library(data.table)
library(IRanges)

dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")

## Get cohort
cohort = purple::query_purity(dbProd)
patientIdLookups = query_patient_id_lookup(dbProd)
patientIds = purple::apply_to_cohort(cohort, function(x) {purple::sample_to_patient_id(x$sampleId, patientIdLookups)})
cohort$patientId <- patientIds$V1
multipleBiopsyCohort = purple::multiple_biopsy(cohort)
multipleBiopsyPatientsId = unique(multipleBiopsyCohort$patientId)
save(cohort, multipleBiopsyCohort, multipleBiopsyPatientsId, file="~/hmf/multipleBiopysyCohort.RData")

## Get somatic variants
multipleBiopsySomaticVariants = purple::query_somatic_variants(dbProd, multipleBiopsyCohort)
save(multipleBiopsySomaticVariants, file="~/hmf/multipleBiopsySomaticVariants.RData")

## Get structural variants
multipleBiopsyStructuralVariants = purple::query_structural_variants(dbProd, multipleBiopsyCohort)
save(multipleBiopsyStructuralVariants, file="~/hmf/multipleBiopsyStructuralVariants.RData")

### Clean up
dbDisconnect(dbProd)
rm(dbProd)

#### LOAD DATA
load(file="~/hmf/multipleBiopysyCohort.RData")
load(file="~/hmf/multipleBiopsySomaticVariants.RData")
load(file="~/hmf/multipleBiopsyStructuralVariants.RData")

### Mutational Signatures
multipleBiopsyMutationalSignature = list()
multipleBiopsyIndelSignature = list()
multipleBiopsySVSignature = list()

for (patientId in multipleBiopsyPatientsId) {
  cat("Processing", patientId, "\n")
  
  # Extract patient data
  patientSampleIds = unique(multipleBiopsyCohort[multipleBiopsyCohort$patientId == patientId, c("sampleId")])
  patientSomaticVariants = scope(multipleBiopsySomaticVariants[multipleBiopsySomaticVariants$sample %in% patientSampleIds,])
  patientStructuralVariants = scope(multipleBiopsyStructuralVariants[multipleBiopsyStructuralVariants$sample %in% patientSampleIds,])
  
  # Signatures
  patientMutationalSignature = purple::mutational_signature_by_scope(patientSomaticVariants)
  patientIndelSignature = purple::indel_signature_by_scope(patientSomaticVariants)
  patientSVSignature = purple::sv_signature_by_scope(patientStructuralVariants)
  
  # Save Signatures
  multipleBiopsyMutationalSignature[[patientId]] <- patientMutationalSignature
  multipleBiopsyIndelSignature[[patientId]] <- patientIndelSignature
  multipleBiopsySVSignature[[patientId]] <- patientSVSignature  

  # Plots
  colnames(patientSomaticVariants)
  
}

save(multipleBiopsyMutationalSignature, multipleBiopsyIndelSignature, multipleBiopsySVSignature, file="~/hmf/multipleBiopsySignatures.RData")


patientId = multipleBiopsyPatientsId[1]
#PLOTTING
signatures = purple::getCOSMICSignatures()
myCOLORS = c("#ff994b","#463ec0","#88c928","#996ffb","#68b1c0","#e34bd9","#106b00","#d10073","#98d76a",
             "#6b3a9d","#d5c94e","#0072e2","#ff862c","#31528d","#d7003a","#323233","#ff4791","#01837a",
             "#ff748a","#777700","#ff86be","#4a5822","#ffabe4","#6a4e03","#c6c0fb","#ffb571","#873659",
             "#dea185","#a0729d","#8a392f")
library(devtools) #; install_github("im3sanger/dndscv")
library(MutationalPatterns)
library(RMySQL)
library(data.table)
library("NMF")
library(grid)
library(gridExtra)
library(ggplot2)

binsize<-0.05
patientSampleIds[1]

patientVariants = add_scope_to_somatic_variants(patientSomaticVariants)
patient1Variants = patientVariants[patientVariants$scope == patientSampleIds[1], ]
patient2Variants = patientVariants[patientVariants$scope == patientSampleIds[2], ]
sharedVariants = patientVariants[patientVariants$scope == "Shared", ]
sharedPloidy = patientVariants[scope == "Shared", .(minPloidy = min(ploidy), maxPloidy = max(ploidy) ), by=.(chromosome, position, alt, ref)]

p1 <- ggplot() + xlim(0,3.5)+ labs(title = patientSampleIds[1],x='Somatic Ploidy') +
  geom_histogram(aes(x=ploidy),data=patient1Variants,fill = "red", alpha = 0.6,binwidth = binsize) +
  geom_histogram(aes(x=ploidy),data=sharedVariants,fill = "blue", alpha = 0.2,binwidth = binsize)

p2 <- ggplot() + xlim(0,3.5)+ labs(title = patientSampleIds[2],x='Somatic Ploidy')+
  geom_histogram(aes(x=ploidy),data=patient2Variants,fill = "red", alpha = 0.6,binwidth = binsize)+
  geom_histogram(aes(x=ploidy),data=sharedVariants,fill = "blue", alpha = 0.2,binwidth = binsize)

p3 <- ggplot()+labs(title="Shared Somatic Variant Ploidy")+geom_point(data=sharedPloidy,aes(minPloidy,maxPloidy))+xlim(0,3)+ylim(0,3)

p4 <-plot_contribution(patientIndelSignature) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10),legend.text=element_text(size=5),axis.title.y = element_text(size=10))+
  ggtitle("Indel Signatures")

p5 <-plot_contribution(patientSVSignature) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10),legend.text=element_text(size=5),axis.title.y = element_text(size=10))+
  ggtitle("Structural Variant Signatures")


patientSignature = fit_to_signatures(patientMutationalSignature[, -c(1, 2)], signatures)
p6 <- plot_contribution(patientSignature$contribution, signatures)+theme(axis.text.x = element_text(angle = 90, hjust = 1,size=10),legend.text=element_text(size=5),axis.title.y = element_text(size=10))+
  scale_fill_manual( values= myCOLORS)+labs(fill="")+ggtitle("Mutational Signatures")

multiplot(p1,p2,p3,p4,p5,p6, cols = 2)

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

patientSVSignature$type


DT = data.table(patientStructuralVariants)

## Length
DT$length <- pmax(DT$startPosition, DT$endPosition) - pmin(DT$startPosition, DT$endPosition)
DT$length <- ifelse(DT$type == "BND", 0, DT$length)
DT$length <- ifelse(DT$type == "INS", 0, DT$length)
DT$length <- cut(DT$length, right = FALSE, breaks = c(-Inf, 1, 1000, 10000, 1e+05, 1e+06, Inf), labels = sv_length_buckets())
DT$type <- factor(paste(DT$type, DT$length, sep="_"), levels=sv_type_length_buckets(), ordered = TRUE)

result = dcast(DT, type ~ scope, value.var = "sampleId", fun.aggregate = length)
result = merge(create_empty_sv_signature(), result, all.x = TRUE)
result[is.na(result)] <- 0

sv_length_buckets<-function() {
  return (c("0", "<1k", "1k-10k", "10k-100k", "100k-1M", ">1M"))
}

sv_type_length_buckets<-function() {
  buckets = length_buckets()[-1]
  labels = c("INS_0","BND_0", paste("DEL", buckets, sep = "_"), paste("DUP", buckets, sep = "_"), paste("INV", buckets, sep = "_"))
  return (labels)
}

create_empty_sv_signature<-function() {
  buckets = sv_type_length_buckets()
  type = factor(buckets,levels=buckets,ordered=TRUE)
  return (data.frame(type = labels, stringsAsFactors = FALSE))
}

