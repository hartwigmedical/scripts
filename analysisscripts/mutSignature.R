library(devtools)#; install_github("im3sanger/dndscv")
library(MutationalPatterns)
library(RMySQL)
library(data.table)
library(grid)
library(gridExtra)
library("NMF")
library(ggplot2)

install.packages("Bio")
library(BiocInstaller)

myCOLORS = c("#ff994b","#463ec0","#88c928","#996ffb","#68b1c0","#e34bd9","#106b00","#d10073","#98d76a",
             "#6b3a9d","#d5c94e","#0072e2","#ff862c","#31528d","#d7003a","#323233","#ff4791","#01837a",
             "#ff748a","#777700","#ff86be","#4a5822","#ffabe4","#6a4e03","#c6c0fb","#ffb571","#873659",
             "#dea185","#a0729d","#8a392f")
OLDmyCOLORS = c("#ff994b","#463ec0","#88c928","#996ffb","#68a100","#e34bc9","#106b00","#d10073","#98d76a",
            "#6b3a9d","#d5c94e","#0072e2","#ff862c","#31528d","#d7003a","#00825b","#ff4791","#01837a",
            "#ff748a","#777700","#ff86be","#4a5822","#ffabe4","#6a4e03","#c6c0fb","#ffb571","#873659",
            "#dea185","#a0729d","#8a392f")

########### FUNCTION DEFINITIONS ######################

standard_mutation<-function(types){
  types = gsub("G>T", "C>A", types)
  types = gsub("G>C", "C>G", types)
  types = gsub("G>A", "C>T", types)
  types = gsub("A>T", "T>A", types)
  types = gsub("A>G", "T>C", types)
  types = gsub("A>C", "T>G", types)
  return(types)
}

standard_context<-function(raw_type, standard_type, context) {
  x = which(raw_type != standard_type)
  y = context[x]
  y = reverse(chartr("ATGC", "TACG", y))
  context[x] = y
  return(context)
}

create_empty_signature2<-function() {
  DF <- data.frame(type=rep("", 96), context=rep("   ", 96), count=rep(0L, 96), stringsAsFactors=FALSE)
  ref_bases = c("C", "T")
  alt_bases = list(c("A", "G", "T"), c("A", "C", "G"))
  bases = c("A", "C", "G", "T")
  for (i in  0:1) {
    ref = ref_bases[i+1]
    for (j in 0:2) {
      alt = alt_bases[[i+1]][j+1]
      type = paste(ref, alt, sep = ">")
      for (x in 0:3) {
        for (y in 0:3) {
          before = bases[x+1]
          after = bases[y+1]
          index = i * 48 + j * 16 + x * 4 + y + 1
          context = paste(before, after, sep=ref)
          DF[index, ] <- list(type, context, 0)
        }
      }
    }
  }
  return (DF)
}

create_empty_signature<-function() {
  DF <- data.frame(type=character(), context=character(), count=double(), stringsAsFactors=FALSE)
  ref_bases = c("C", "T")
  bases = c("A", "C", "G", "T")
  for (ref in ref_bases) {
    for (alt in bases) {
      if (alt != ref) {
        type = paste(ref, alt, sep = ">")
        for (before in bases) {
          for (after in bases) {
            context = paste(before, after, sep=ref)
            DF = rbind(DF, data.frame(type, context, count=0, stringsAsFactors=FALSE))
          }
        }
      }
    }
  }
  return(DF)
}

select_cancer_types<-function(dbConnect){
  query = "select distinct cancerType from clinical c where cancerType is not null;"
  return (dbGetQuery(dbConnect, query))
}

select_DRUP<-function(dbConnect, type){
  query = paste(
    "select p.sampleId,'NONE' as biopsySite from purity p where ",
    "qcStatus = 'PASS' and status <> 'NO_TUMOR' and sampleId like 'DRUP%'",
    sep = "")
  return (dbGetQuery(dbConnect, query))
}

select_cohort_with_subclones<-function(dbConnect, type,minMutationCount=0,subclonalProportion=0){
  query = paste(
    "select s.sampleId,biopsySite from hmfpatients.clinical c, sample s, purity p,somaticVariant sv where s.sampleId = c.sampleId and s.sampleId = p.sampleId ",
    "and sv.sampleId = s.sampleId ",
    "and qcStatus = 'PASS' and status <> 'NO_TUMOR' and cancerType like '%",
    type,"%' AND filter = 'PASS' ",
    "group by s.sampleId,biopsySite having sum(if(clonality='SUBCLONAL',1,0))/count(*) >",subclonalProportion,
    " and count(*) >",minMutationCount,
    sep = "")
  return (dbGetQuery(dbConnect, query))
}

sample_signature<-function(dbConnect, sample, empty_signature,clonality=""){
  query = paste(
    "select sampleId, trinucleotideContext as context, concat(ref,'>', alt) as snv, count(*) as count from somaticVariant where sampleId = '",
    sample,
    "' and filter = 'PASS' and length(alt) = length(ref) and length(alt) = 1 and trinucleotideContext not like '%N%'",
    sep = "")
  if (clonality!=""){
    query=paste(query," AND clonality = '",clonality,"'",sep="")
  }
  query=paste(query," group by 1, 2, 3",sep="")
  raw_data = dbGetQuery(dbConnect, query)
  raw_types = raw_data$snv
  standard_types = standard_mutation(raw_types)
  raw_context = raw_data$context
  standard_context = standard_context(raw_types, standard_types, raw_context)

  DT=data.table(type=standard_types, context=standard_context, count=raw_data$count)
  DT=rbind(DT, empty_signature)
  grouped=DT[, .(count=sum(count)), keyby=.(type,context)]
  return(grouped)
}

cohort_signature<-function(dbConnect, cohort,clonality=""){
  empty_signature = create_empty_signature()
  result = matrix(, nrow = 96, ncol = 0)
  for (i in 1:nrow(cohort)) {
    sample = cohort[i,]$sampleId
    signature = sample_signature(dbConnect, sample, empty_signature,clonality)
    total_count = sum(signature$count)
    if (total_count > 0)
    {
      print(paste("Processing sample", sample))
      result = cbind(result, signature$count)
      rownames(result) <- signature$context
      colnames(result)[ncol(result)] <- paste(sample,total_count,substr(cohort[i,]$biopsySite,1,8))
    }
    else
    {
      print(paste("No mutations available in sample", sample))
    }
  }
  return(result)
}

getCOSMICSignatures<-function(){
  sp_url = "http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"
  cancer_signatures = read.table(sp_url, sep = "\t", header = T)
  # reorder (to make the order of the trinucleotide changes the same)
  cancer_signatures = cancer_signatures[order(cancer_signatures[,1]),]
  # only signatures in matrix
  cancer_signatures = as.matrix(cancer_signatures[,4:33])
}

calculate_signatures<-function(dbConnect,cohort,cancer_signatures){
  mutation_matrix= cohort_signature(dbConnect, cohort[1:nrow(cohort),])
  fit_res = fit_to_signatures(mutation_matrix, cancer_signatures)
  fit_contribution<-fit_res$contribution[, order(colnames(fit_res$contribution),decreasing=F),drop=FALSE]
  fit_contribution[prop.table(fit_contribution, margin=2)<0.03 | fit_contribution<100]<-0
  orderVector<-colSums(fit_contribution)
  fit_contribution[, order(orderVector,decreasing=F),drop=FALSE]
}

plot_fitted_signatures<-function(dbConnect, fit_contribution,cancer_signatures,cancerType,chart_mode="absolute",writePDF=False){
  plots=list()
  if (cancerType=='Breast') {
    highMutationLoadCutOff=6000
  } else if (cancerType == 'Colorectal') {
    highMutationLoadCutOff=15000
  } else {
    highMutationLoadCutOff=10000  
  }
  
  veryHighMutationLoadCutOff=40000

  if (writePDF){
    pdf(file=paste(cancerType,"ALLRelative.pdf",sep=""),width=10)
  }
  selectRow = which(rowSums(fit_contribution)>0)
  if (any(colSums(fit_contribution)<highMutationLoadCutOff)) {
    p1<-plot_contribution(fit_contribution[selectRow,colSums(fit_contribution)<=highMutationLoadCutOff,drop=F], cancer_signatures[,selectRow],coord_flip = F, mode = chart_mode)+
      theme(axis.text.x = element_text(angle = 90, hjust = 1,size=7),legend.text=element_text(size=5),axis.title.y = element_text(size=6))+
      scale_fill_manual( values= myCOLORS[selectRow])+labs(fill="")+ggtitle(paste(cancerType,"# of mutations <",highMutationLoadCutOff))
    print(p1)
  }
  if (any(colSums(fit_contribution)>highMutationLoadCutOff&colSums(fit_contribution)<veryHighMutationLoadCutOff)) {
    p2<-plot_contribution(fit_contribution[selectRow,colSums(fit_contribution)>highMutationLoadCutOff&colSums(fit_contribution)<veryHighMutationLoadCutOff,drop=F], cancer_signatures[,selectRow],coord_flip = F, mode = chart_mode)+
      theme(axis.text.x = element_text(angle = 90, hjust = 1,size=7),legend.text=element_text(size=5),axis.title.y = element_text(size=6))+
      scale_fill_manual( values= myCOLORS[selectRow])+labs(fill="")+ggtitle(paste(cancerType,"# of mutations: ",highMutationLoadCutOff,"to",veryHighMutationLoadCutOff))
    print(p2)
  }
  if (any(colSums(fit_contribution)>veryHighMutationLoadCutOff)) {
    p3<-plot_contribution(fit_contribution[selectRow,colSums(fit_contribution)>=veryHighMutationLoadCutOff,drop=F], cancer_signatures[,selectRow],coord_flip = F, mode = chart_mode)+
      theme(axis.text.x = element_text(angle = 90, hjust = 1,size=7),legend.text=element_text(size=5),axis.title.y = element_text(size=6))+
      scale_fill_manual( values= myCOLORS[selectRow])+labs(fill="")+ggtitle(paste(cancerType,"# of mutations >",veryHighMutationLoadCutOff))
    print(p3)
  }
  if (writePDF){
    dev.off()
  }
}

clonal_vs_subclonal_signatures<-function(dbConnect, cohort,cancer_signatures,cancerType,chart_mode="absolute",writePDF=False){
  fit_contribution=list()
  i=1
  plots=list()
  for (clonality in c("CLONAL","SUBCLONAL")){
    mutation_matrix= cohort_signature(dbConnect, cohort[1:nrow(cohort),],clonality)
    fit_res = fit_to_signatures(mutation_matrix, cancer_signatures)
    fit_contribution[[i]]<-fit_res$contribution[, order(colnames(fit_res$contribution),decreasing=F),drop=FALSE]
    fit_contribution[[i]][prop.table(fit_contribution[[i]], margin=2)<0.03 | fit_contribution[[i]]<100]<-0
    i=i+1
  }
  selectRow = which(rowSums(fit_contribution[[1]])+rowSums(fit_contribution[[2]])>0)
  orderVector<-colSums(fit_contribution[[1]])#+colSums(fit_contribution[[2]])

  
  #myColors<-hcl(h = seq(15, 375, length =  31), l = 65, c = 100)[1:30]
  for (i in 1:2){
    fit_contribution[[i]]<-fit_contribution[[i]][, order(orderVector,decreasing=F),drop=FALSE]
    plots[[i]]<-plot_contribution(fit_contribution[[i]][selectRow,,drop=F], cancer_signatures[,select],coord_flip = F, mode = chart_mode)+
      theme(axis.text.x = element_text(angle = 90, hjust = 1,size=6),legend.text=element_text(size=5),axis.title.y = element_text(size=6))+
      scale_fill_manual( values= myCOLORS[selectRow])+labs(fill="")
  }
  if (writePDF){
    pdf(file=paste(cancerType,"subClonal.pdf",sep=""),width=10)
    multiplot(plots[[1]],plots[[2]])
    dev.off()
  }
  else {
    multiplot(plots[[1]],plots[[2]])
  }

}

###########################
########################MAIN CODE###################
# Create mutational matrix
cancer_signatures = getCOSMICSignatures()

# Retrieve cancer types
dbConnect = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
cancerTypes = select_cancer_types(dbConnect)
dbDisconnect(dbConnect)

# PLOT subclonal
dbConnect = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")


for (cancerType in cancerTypes$cancerType) {
  print(paste("running for type:",cancerType))
  cohort = select_cohort_with_subclones(dbConnect, cancerType,3000,0.10)
  if (length(cohort)>0){
    clonal_vs_subclonal_signatures(dbConnect,cohort,cancer_signatures,cancerType,"relative",TRUE)
    
  }
}
dbDisconnect(dbConnect)

cohort

# PLOT ALL
dbConnect = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis")
fitted_signatures=list()
#cohort = select_DRUP(dbConnect, "none")
for (cancerType in cancerTypes$cancerType) {
  print(paste("running for type:",cancerType))
  cohort = select_cohort(dbConnect, cancerType)
  if (length(cohort)>0){
    fitted_signatures[[cancerType]]<-calculate_signatures(dbConnect,cohort,cancer_signatures)
    
  }
}
for (cancerType in cancerTypes$cancerType) {
  plot_fitted_signatures(dbConnect,fitted_signatures[[cancerType]],cancer_signatures,cancerType,"relative",T)
}  

# PLOT 2 samples
dbConnect = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis")
fitted_signatures=list()
cancerType="none"
cohort = data.frame(sampleId=c("DRUP01050010T","CPCT02050172T"),biopsySite=c(cancerType,cancerType))
fitted_signatures[[cancerType]]<-calculate_signatures(dbConnect,cohort,cancer_signatures)
plot_fitted_signatures(dbConnect,fitted_signatures[[cancerType]],cancer_signatures,cancerType,"relative",T)

dbDisconnect(dbConnect)
######################################################################
#SINGLE SAMPLE
dbConnect = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysis")
cancerType="Brain"
cohort=data.frame(sampleId="CPCT02030265TII",biopsySite=cancerType)
fitted_signatures[[cancerType]]<-calculate_signatures(dbConnect,cohort,cancer_signatures)
fitted_signatures
plot_contribution(fitted_signatures, cancer_signatures,coord_flip = F, mode = chart_mode)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=7),legend.text=element_text(size=5),axis.title.y = element_text(size=6))+
  scale_fill_manual( values= myCOLORS[selectRow])+labs(fill="")
#MUT_MATRIX
mut_matrix=list()
for (cancerType in cancerTypes$cancerType) {
  print(paste("running for type:",cancerType))
  cohort = select_cohort(dbConnect, cancerType)
  if (length(cohort)>0){
    mut_matrix[[cancerType]]<-cohort_signature(dbConnect, cohort[1:nrow(cohort),])
    
  }
}


##### LOCAL STORAGE ##################
dataFile = "~/hmf/mutMatrix.RData"

# SAVE TO FILE
save(mut_matrix, file = dataFile)

# LOAD FROM FILE
load(dataFile)


######## RANDOM TESTING UNDER HERE #############



mut_mat <- mut_matrix + 0.00001
estimate <- nmf(mut_mat, rank=2:20, method="brunet", nrun=10, seed=123456)
plot(estimate)
nmf_res <- extract_signatures(mut_mat, rank = 20, nrun = 10)

colnames(nmf_res$signatures) <- 1:20
rownames(nmf_res$contribution) <- 1:20
plot_96_profile(nmf_res$signatures)
pc1 <- plot_contribution(nmf_res$contribution, nmf_res$signature,mode = "relative")
pc2 <- plot_contribution(nmf_res$contribution, nmf_res$signature, mode = "absolute")
grid.arrange(pc1, pc2)


ggplot(data=mutation_matrix), 
      aes(x=CPCT02020412T) + stat_ecdf(geom = "step") + labs(title="Test")


library('proxy') # Library of similarity/dissimilarity measures for 'dist()'
library(reshape2)
tMatrix=t(mutation_matrix)
sim_matrix = dist(tMatrix/rowSums(tMatrix),method="cosine",pairwise = TRUE,diag = TRUE,upper= TRUE,)
View(as.matrix(sim_matrix))
clusterOutput=hclust(sim_matrix,method="median")
par(cex=0.6, mar=c(5, 8, 4, 1))
plot(clusterOutput, xlab="", ylab="", main="", sub="", axes=FALSE)
par(cex=1)
axis(2)


#########################

load("~/Dropbox/HMF Australia team folder/RData/mutSignature2.RData")


# Fit the signature for all samples and clonalities
# For each cancerType:
# S

calculate_signatures<-function(dbConnect,cohort,cancer_signatures){
  mutation_matrix= cohort_signature(dbConnect, cohort[1:nrow(cohort),])
  fit_res = fit_to_signatures(mutation_matrix, cancer_signatures)
  fit_contribution<-fit_res$contribution[, order(colnames(fit_res$contribution),decreasing=F),drop=FALSE]
  fit_contribution[prop.table(fit_contribution, margin=2)<0.03 | fit_contribution<100]<-0
  orderVector<-colSums(fit_contribution)
  fit_contribution[, order(orderVector,decreasing=F),drop=FALSE]
}

mutation_matrix


?fit_to_signatures

