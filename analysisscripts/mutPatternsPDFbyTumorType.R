library(MutationalPatterns)
library(BSgenome)
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
library("NMF")

loadVCFsAndFilterforAutosomes<-function(vcfDir){
  vcf_files = list.files(vcfDir, pattern = ".vcf", full.names = TRUE)
  sample_names<-tools::file_path_sans_ext(basename(vcf_files))
  vcfs = read_vcfs_as_granges(vcf_files, sample_names, genome = "hg19")### GENOME??????########
  auto = extractSeqlevelsByGroup(species="Homo_sapiens",style="UCSC",group="auto")
  vcfs = lapply(vcfs, function(x) keepSeqlevels(x, auto))
  vcfs
}

getCOSMICSignatures<-function(vcfDir){
  sp_url = "http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"
  cancer_signatures = read.table(sp_url, sep = "\t", header = T)
  # reorder (to make the order of the trinucleotide changes the same)
  cancer_signatures = cancer_signatures[order(cancer_signatures[,1]),]
  # only signatures in matrix
  cancer_signatures = as.matrix(cancer_signatures[,4:33])
}

cancer_signatures = getCOSMICSignatures()

####### Find Data, Load VCFs and ###########
searchDir = "/data/cpct/runs"
fileNamePattern = "melted.vcf"
patientlistFile = "/home/peter/tumor_data.csv"
tumor_data <- read.csv(file=patientlistFile, header=TRUE, sep=",",strip.white=TRUE)[,1:2]
my_tumor_list <- split(tumor_data,tumor_data[2])
allVCFFiles<-list.files(path=searchDir,pattern=fileNamePattern,recursive=TRUE)
for (name in names(my_tumor_list)) {
  if (name != "") {
    vcf_files<-allVCFFiles[grepl((paste(my_tumor_list[[name]][,1],collapse="|")),allVCFFiles)]
    if (length(vcf_files)>0 & length(vcf_files)<10){
      print(paste(name,length(vcf_files)))
      #print(vcf_files)
      vcf_files<-paste(searchDir,vcf_files,sep="")
      vcfs = read_vcfs_as_granges(vcf_files, tools::file_path_sans_ext(basename(vcf_files)), genome = "hg19")### GENOME??????########
      auto = extractSeqlevelsByGroup(species="Homo_sapiens",style="UCSC",group="auto")
      vcfs = lapply(vcfs, function(x) keepSeqlevels(x, auto))
      
      # Fit & plot to pdf
      fit_res = fit_to_signatures(mut_matrix(vcf_list = vcfs, ref_genome = ref_genome), cancer_signatures)
      select = which(rowSums(fit_res$contribution) > 0.02 * sum(fit_res$contribution))  #2% contribution cutoff
      pc<-plot_contribution(fit_res$contribution[select,,drop=FALSE], cancer_signatures[,select], coord_flip = F, mode = "absolute")
      pdf(file=paste("test",name,".pdf",sep=""),width=10)#, width=10, height=2, pointsize=6, useDingbats=FALSE)
      print(pc)
      dev.off()
    }
  }
}