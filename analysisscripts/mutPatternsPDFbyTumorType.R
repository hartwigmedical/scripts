library(MutationalPatterns)
library(BSgenome)
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
library(VariantAnnotation)
library(ggplot2)

read_vcfs_as_granges <-function (vcf_files, sample_names, genome = "-", style = "UCSC") 
{
  if (length(vcf_files) != length(sample_names)) 
    stop("Provide the same number of sample names as VCF files")
  vcf_list <- GRangesList(lapply(vcf_files, function(file) {
    print(paste(file,"check1"))
    print(Sys.time())
    vcf <- rowRanges(readVcf(file, genome))
    seqlevelsStyle(vcf) <- style
    rem <- which(all(!(!is.na(match(vcf$ALT, DNA_BASES)) & 
                         !is.na(match(vcf$REF, DNA_BASES)) & (lengths(vcf$ALT) == 
                                                                1))))
    if (length(rem) > 0) {
      vcf = vcf[-rem]
    }
    print(paste(file,"check2"))
    print(Sys.time())
    return(vcf)
  }))
  names(vcf_list) <- sample_names
  return(vcf_list)
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

#searchDir = "~/hmf/analyses/mutPatternsTest/temp/"
searchDir = "/data/experiments/consensus_filtered/"
fileNamePattern = ".vcf"
#patientlistFile = "~/hmf/analyses/mutPatternsTest/tumor_data.csv"#/home/peter/tmp/ecrf_dump_for_patients.csv"
patientlistFile = "/home/peter/tmp/ecrf_dump_for_patients.csv"

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
    
    if (length(vcf_files)>0){   
      print(vcf_files)
      # Load VCFs
      vcf_files<-paste(searchDir,vcf_files,sep="")
      #print(Sys.time())
      vcfs = read_vcfs_as_granges(vcf_files, tools::file_path_sans_ext(basename(vcf_files)), genome = "hg19")### GENOME??????########
      #print(Sys.time())
      auto = extractSeqlevelsByGroup(species="Homo_sapiens",style="UCSC",group="auto")
      vcfs = lapply(vcfs, function(x) keepSeqlevels(x, auto))
      
      # Fitting
      mutMatrix<-mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
      dimnames(mutMatrix)[[2]]<-substr(dimnames(mutMatrix)[[2]],1,12)
      fit_res = fit_to_signatures(mutMatrix, cancer_signatures)
      select = which(rowSums(fit_res$contribution) > 0.02 * sum(fit_res$contribution))  #2% contribution cutoff
      
      # Plot to PDF
      pc<-plot_contribution(fit_res$contribution[select,,drop=FALSE], cancer_signatures[,select], coord_flip = F, mode = "absolute")+ggtitle(name)+theme(axis.text.x = element_text(angle = 90, hjust = 1))
      pdf(file=paste(name,".pdf",sep=""),width=10)#, width=10, height=2, pointsize=6, useDingbats=FALSE)
      print(pc)
      dev.off()
    }
  }
}

