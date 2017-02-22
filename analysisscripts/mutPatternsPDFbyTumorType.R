library(MutationalPatterns)
library(BSgenome)
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

getCOSMICSignatures<-function(vcfDir){
  sp_url = "http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"
  cancer_signatures = read.table(sp_url, sep = "\t", header = T)
  # reorder (to make the order of the trinucleotide changes the same)
  cancer_signatures = cancer_signatures[order(cancer_signatures[,1]),]
  # only signatures in matrix
  cancer_signatures = as.matrix(cancer_signatures[,4:33])
}

cancer_signatures = getCOSMICSignatures()

searchDir = "/data/cpct/runs/"
fileNamePattern = "melted.vcf"
patientlistFile = "/home/peter/tumor_data.csv"

#Load patients and tumor types from ECRF
tumor_data <- read.csv(file=patientlistFile, header=TRUE, sep=",",strip.white=TRUE)[,1:2]
my_tumor_list <- split(tumor_data,tumor_data[2])

# Find all somatic vcfs in the health check runs
allVCFFiles<-list.files(path=searchDir,pattern=fileNamePattern,recursive=TRUE)

for (name in names(my_tumor_list)) {
  if (name != "") {
    # filter for vcfs with the  tumor type
    vcf_files<-allVCFFiles[grepl((paste(my_tumor_list[[name]][,1],collapse="|")),allVCFFiles)]
    print(paste(name,length(vcf_files)))
    if (length(vcf_files)>0 & length(vcf_files)<10){   ### TEMP CONDITION - ONLY DO SMALL COHORTS FOR TESTING
      
      # Load VCFs
      vcf_files<-paste(searchDir,vcf_files,sep="")
      vcfs = read_vcfs_as_granges(vcf_files, tools::file_path_sans_ext(basename(vcf_files)), genome = "hg19")### GENOME??????########
      auto = extractSeqlevelsByGroup(species="Homo_sapiens",style="UCSC",group="auto")
      vcfs = lapply(vcfs, function(x) keepSeqlevels(x, auto))
      
      # Fitting
      fit_res = fit_to_signatures(mut_matrix(vcf_list = vcfs, ref_genome = ref_genome), cancer_signatures)
      select = which(rowSums(fit_res$contribution) > 0.02 * sum(fit_res$contribution))  #2% contribution cutoff
      
      # Plot to PDF
      pc<-plot_contribution(fit_res$contribution[select,,drop=FALSE], cancer_signatures[,select], coord_flip = F, mode = "absolute")
      pdf(file=paste("test",name,".pdf",sep=""),width=10)#, width=10, height=2, pointsize=6, useDingbats=FALSE)
      print(pc)
      dev.off()
    }
  }
}