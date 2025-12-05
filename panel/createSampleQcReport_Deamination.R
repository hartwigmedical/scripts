#!/usr/bin/env Rscript
### Gene Report

library(tidyr)
library(dplyr)
library(stringi)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(cowplot)
library(VariantAnnotation) # required for readVcf

# Register functions
check_file<-function(fileName)
{
  if(!file.exists(fileName))
  {
    print(sprintf('Missing file: %s', fileName))
    stop()
  }
}

load_file<-function(fileName,fileDesc,delim='\t')
{
  data = read.csv(fileName,sep=delim)
  print(sprintf('loaded %d %s records from %s', nrow(data), fileDesc, fileName))
  return (data)
}

vcf_data_frame<-function(vcf) 
{
  vcf.rowRanges = rowRanges(vcf)
  vcf.info = info(vcf)
  vcf.alt = CharacterList(alt(vcf))
  vcf.alt[elementNROWS(vcf.alt) > 1L ] <- lapply(vcf.alt[elementNROWS(vcf.alt) > 1L ], paste0, collapse=",")
  
  vcf.df = data.frame(
    Chromosome = seqnames(vcf), 
    Pos = start(vcf), 
    Ref = as.character(ref(vcf)), 
    Alt = as.character(vcf.alt),  
    Filter = as.character(vcf.rowRanges$FILTER))

  return (vcf.df)
}

runSampleQcDeamination<-function()
{


    # Parse and check inputs
    

#     Parse and check inputs
    args <- commandArgs(trailingOnly = TRUE)
    print(args)
    
    if(length(args) < 3)
    {
      print("Requires arguments: 1=SampleId, 2=DataDir, 3=RunDir, [4=ResourceDir (default=/data/resources/public)]")
      stop()
    }

    sampleId <- args[1]
    sampleDataDir <- args[2]
    runDir = args[3]

    # Optional argument with default
    resourceDir = ifelse(length(args) >= 4, args[4], "/data/resources/public")

    cohortMedianDepthFile = paste0(resourceDir, "/target_regions/38/target_regions_exon_relative_depth.38.csv")
    cobaltRegionsFile = paste0(resourceDir, "/target_regions/38/target_regions_normalisation.38.tsv")
    ensembl_gene_data = paste0(resourceDir, "/ensembl_data_cache/38/ensembl_gene_data.csv")
    driverGenePanel = paste0(resourceDir, "/gene_panel/38/DriverGenePanel.38.tsv")

    print(sampleId)

    #sampleDataDir = paste0(runDir, "Test/")
#    sampleDataDir = paste0(runDir, "230501-msi-valexp09-func-116/",sampleId,"/")
    sampleExonCoverageFile = paste0(sampleDataDir,'/sage_somatic/',sampleId,'.sage.exon.medians.tsv')
    sampleGeneCnFile = paste0(sampleDataDir,'/purple/',sampleId,'.purple.cnv.gene.tsv')
    samplePurityFile = paste0(sampleDataDir,'/purple/',sampleId,'.purple.purity.tsv')
    sampleSomaticVcf = paste0(sampleDataDir,'/purple/',sampleId,'.purple.somatic.vcf.gz')




    check_file(cohortMedianDepthFile)
    check_file(cobaltRegionsFile)
    check_file(ensembl_gene_data)
    check_file(sampleExonCoverageFile)
    check_file(sampleGeneCnFile)
    check_file(samplePurityFile)
    check_file(sampleSomaticVcf)
    check_file(driverGenePanel)

    # set plotting features
    defaultFontSize=8
    theme_set(theme_bw(base_size=defaultFontSize) +
                theme(axis.text=element_text(size=defaultFontSize),
                      axis.title=element_text(size=defaultFontSize),
                      axis.ticks.y=element_blank(),
                      legend.text=element_text(size=defaultFontSize),
                      legend.title=element_text(size=defaultFontSize),
                      panel.grid.minor.x=element_blank()))


    # DATA PREPARATION

    # load reference data
    cohortMedianDepth = load_file(cohortMedianDepthFile,"cohort median depth",',')
    colnames(cohortMedianDepth) = c('Gene','ExonRank','CohortRelDepth')

    cobaltRegions = load_file(cobaltRegionsFile,'cobalt region enrichment')
    colnames(cobaltRegions) = c('Chromosome','Position','RelativeEnrichment')

    ensemblGenes = load_file(ensembl_gene_data,'ensembl gene data', delim=',')
    driverGenePanel = load_file(driverGenePanel, 'driver gene panel')

    # arrange genes on alphabetical order of gene name
    # genesList = genesList %>% arrange(Gene) %>% mutate(GeneIndex=row_number())

    # sort ensembl on chromosomal order
    ensemblGenes = ensemblGenes %>% mutate(ChromIndex=as.integer(ensemblGenes$Chromosome))
    ensemblGenes$ChromIndex[which(ensemblGenes$Chromosome=='X')] <- 23
    ensemblGenes$ChromIndex[which(ensemblGenes$Chromosome=='Y')] <- 24
    ensemblGenes = ensemblGenes %>% arrange(ChromIndex, GeneStart)
    # check if correct: table(ensemblGenes$Chromosome, ensemblGenes$ChromIndex)

    # filter genes from driver gene panel - annotated with likelihoodType (TSG / ONCO)
    # driverGenes = ensemblGenes %>% filter(GeneName %in% driverGenePanel$gene)
    colnames(ensemblGenes)[2] <- 'gene'
    driverGenes = dplyr::inner_join(ensemblGenes, driverGenePanel, by="gene")

    # genesOfInterest = c('MPL','NRAS','ALK','IDH1','ERBB4','VHL','MLH1','CTNNB1','PIK3CA','FGFR3','PDGFRA','KIT','KDR','FBXW7','APC','CSF1R','NPM1',
    #                     'ROS1','EGFR','MET','SMO','BRAF','EZH2','FGFR1','JAK2','CDKN2A','GNAQ','ABL1','NOTCH1','RET','PTEN','FGFR2','HRAS',"ATM",
    #                     'KRAS','PTPN11','HNF1A','POLE','FLT3','RB1','AKT1','IDH2','CDH1','TP53','ERBB2','SMAD4','STK11','GNA11','JAK3','SRC','GNAS','SMARCB1')

    # filter genes of interest for OncoPanel. Check: "Proposed validation gene list" in Oncopanel > CNV validation
    # 15 genes for DELs
    reportGenesDeletions = c('PALB2', 'RAD51B', 'RAD51C', 'BRCA1', 'BRCA2', 'CDKN2A', 'TP53', 'CREBBP', 'EP300', 'ARID1A', 'PTEN', 'RB1','MTAP','KEAP1','SMAD4')
    # 17 genes for AMPS
    reportGenesAmplification = c('AR','CCNE1','EGFR','ERBB2','FGFR2','FLT4','KDR','KIT','MDM2','MET','MYC','PDGFRA','PDGFRB','PIK3CA','RAF1','ROS1','FGFR1')
    # notes: MST1R is removed from the list, since it is not in Hartwig's driver gene panel

    # these genes do not work: "CCNE1" "MDM2"  "MYC"
#    reportGenesAmplification = c('AR','EGFR','ERBB2','FGFR2','FLT4','KDR','KIT','MET','PDGFRA','PDGFRB','PIK3CA','RAF1','ROS1', 'FGFR1')
    # all reportAmplification=TRUE from driverGenePanel
    # reportAmplification = driverGenes %>% filter(reportAmplification=='true' & likelihoodType=='ONCO')
    genesOfInterest = c(reportGenesDeletions, reportGenesAmplification)

    genesList = driverGenes %>% filter(gene %in% genesOfInterest)
    #likelihoodType = genesList$likelihoodType
    likelihoodType = as.data.frame(cbind(genesList$gene, genesList$likelihoodType))
    colnames(likelihoodType) = c('Gene', 'likelihoodType')
    # TODO: color genes in legend by likelihoodType (TSG/ONCO)

    genesList = as.data.frame(genesList$gene)
    colnames(genesList) = c('Gene')
    genesList = genesList %>% mutate(GeneIndex=row_number())

    #genesList = as.data.frame(cbind(genesList$gene, genesList$likelihoodType))
    #colnames(genesList) = c('Gene', 'likelihoodType')
    #genesList = genesList %>% mutate(GeneIndex=row_number())

    exons = load_file(sampleExonCoverageFile,'sample exon coverage')
    colnames(exons) = c('Gene','Chromosome','PosStart','PosEnd','ExonRank','MedianDepth')

    exons = exons %>% filter(Gene %in% genesOfInterest)

    # genesOfInterest[genesOfInterest %in% unique(exons$Gene) == F]
    # TODO: WHY ARE THERE NO EXONS FOR: "CCNE1" "MDM2"  "MYC"

    exons = merge(exons,genesList,by='Gene',all.x=T)
    exons = merge(exons,cohortMedianDepth,by=c('Gene','ExonRank'),all.x=T)
    exons = exons %>% mutate(MedianSampleDepth=median(MedianDepth),
                             RelDepth=MedianDepth/MedianSampleDepth,
                             NormDepth=ifelse(CohortRelDepth>0,MedianDepth/CohortRelDepth,MedianDepth))
    exons = exons %>% mutate(Position=floor(PosEnd/1000)*1000+1)
    exons = merge(exons,cobaltRegions,by=c('Chromosome','Position'),all.x=T)

     geneCN = load_file(sampleGeneCnFile,'sample gene copy number') %>% filter(transcriptId!='ENST00000579755')
     purity = load_file(samplePurityFile,'sample purity')

    geneCN = geneCN %>% filter(gene %in% genesOfInterest) %>% dplyr::select(Gene=gene,MinCopyNumber=minCopyNumber,MaxCopyNumber=maxCopyNumber)
    geneCN = cbind(geneCN,purity %>% dplyr::select(Ploidy=ploidy))
    geneCN = geneCN %>% mutate(RelCopyNumber=pmin(pmax(MinCopyNumber/Ploidy*0.98,0.13),15.8),
                               MaxCopyNumber=pmax(pmin(MaxCopyNumber/Ploidy*1.02,15.9),0.135))

    geneCN = merge(geneCN,genesList,by='Gene',all.x=T)

    nonGeneMargins=c(0,0,0,-4)

    # Relative Copy Number
    geneCnPlot = ggplot(geneCN,aes(x=reorder(Gene,-GeneIndex),y=log(RelCopyNumber,2))) +
      geom_point(colour='blue',size=1) +
      geom_point(aes(x=Gene,y=log(MaxCopyNumber,2)),colour='blue',size=1) +
      geom_linerange(aes(ymin = log(RelCopyNumber,2), ymax = log(MaxCopyNumber,2)),colour='black',size=0.25)+
      scale_y_continuous(limits = c(-3.0,4.0)) +
      coord_flip() +
      theme(legend.position='none',axis.text.y=element_blank(),text = element_text(size = 12),axis.text.x = element_text(size = 12),plot.margin=unit(nonGeneMargins,'cm')) +
      labs(title='Gene Copy Number (log2)', x='',y='')

    # Exon Depth Normalised
    exons = exons %>% mutate(ValidNormDepth=(!is.na(RelativeEnrichment)))

    exonColours = c('blue','blue')

    exonNormPlot = ggplot(exons, aes(x=reorder(Gene,-GeneIndex), y=log(RelDepth,2))) + # y=pmax(NormDepth,0.1)
      scale_y_continuous(limits = c(-3,3)) +
      geom_jitter(width=0.05,size=0.1,alpha=0.5,aes(colour=ValidNormDepth)) +
      scale_color_manual(values = exonColours) +
      coord_flip() +
      theme(legend.position='none',axis.text.y=element_blank(),text = element_text(size = 12),axis.text.x = element_text(size = 12),plot.margin=unit(nonGeneMargins,'cm')) +
      labs(title='Exon Normalised Coverage (log2)', x='',y='')

    # medianDepth
    exonMedianPlot = ggplot(exons, aes(x=reorder(Gene,-GeneIndex), y=log(MedianDepth,2))) +
      scale_y_continuous(limits = c(0,16)) +
      geom_jitter(width=0.05,size=0.1,alpha=0.5,aes(colour=ValidNormDepth)) +
      scale_color_manual(values = exonColours) +
      coord_flip() +
      geom_hline(yintercept = log2(100), color="red") +
      theme(legend.position='none',axis.text.y=element_blank(),text = element_text(size = 12),axis.text.x = element_text(size = 12),plot.margin=unit(nonGeneMargins,'cm')) +
      labs(title='Exon Median Coverage (log2)', x='',y='',size=12)

    # SNVs by gene and coding effect

    somaticVcf = readVcf(sampleSomaticVcf)
    somaticVcf=somaticVcf[rowRanges(somaticVcf)$FILTER=="PASS",] #only consider variants annotated with PASS
    somaticVariants = vcf_data_frame(somaticVcf)

    # extract gene impact
    infoData = as.data.frame(info(somaticVcf))
    impactData = infoData %>% mutate(Impact=ifelse(is.na(IMPACT)|as.character(IMPACT)=='character(0)','',as.character(IMPACT))) %>% dplyr::select(Impact,Hotspot=HOTSPOT)
    colnames(impactData) = c('Impact','Hotspot')
    rownames(impactData) = NULL
    impactData = impactData %>% mutate(Impact=stri_replace_all_fixed(Impact,'c(',''))
    impactData = impactData %>% mutate(Impact=stri_replace_all_fixed(Impact,')',''))
    impactData = impactData %>% mutate(Impact=stri_replace_all_fixed(Impact,'"',''))

    impactFields = c('Gene','Transcript','CanonicalEffect','CanonicalCodingEffect','SpliceRegion','HgvsCodingImpact','HgvsProteinImpact','OtherReportableEffects','WorstCodingEffect','GenesAffected')
    impactData = impactData %>% separate(Impact,impactFields,sep=', ')
    impactData = impactData %>% mutate(Gene=ifelse(is.na(Gene),'',as.character(Gene)),CanonicalCodingEffect=ifelse(is.na(CanonicalCodingEffect),'NONE',as.character(CanonicalCodingEffect)))

    # extract AF + DP
    genotypeData = geno(somaticVcf)
    dpData = cbind(genotypeData$AF, genotypeData$DP)
    colnames(dpData) = c('AF', 'DP')
    rownames(dpData) = NULL
    afData <- dpData
    varData = cbind(somaticVariants,impactData %>% dplyr::select(Gene,CanonicalCodingEffect,Hotspot))
    varData = cbind(varData,afData)

    # deamination type = C>T/G>A
    varDataPlus = varData %>% mutate(Deamination=case_when((Ref=='C' & Alt=='T') ~ 'TRUE',
                                                           (Ref=='G' & Alt=='A') ~ 'TRUE',
                                                            TRUE ~ 'FALSE'))
    legendCodingEffectColours =c('red','black')
    varGeneData <- varDataPlus

    # Plot of all PASSED variants in the panel by variant type (Deamination)
    # so not filtered for Genes of interest / Reportability / Impact
    sampleTitle = paste0('QC: Variant Base Substitutions - ', sampleId, ' ( mTCP ', 100*(purity$purity), '% )' )

    variantPlot = ggplot(varGeneData, aes(x=log(DP,10),y=AF)) +
         geom_point(aes(color=Deamination)) +
         scale_color_manual(values=c("black", "red")) +
         theme(axis.ticks.y = element_blank()) +
         scale_x_continuous(limits = c(0.5, 3.5)) +
         scale_y_continuous(limits = c(0, 1.0)) +
         coord_flip() +
         theme(panel.grid.minor.x = element_blank(),text = element_text(size = 12)) +
         theme( axis.title=element_text(size=12), axis.text = element_text(size = 12)) +
         theme(legend.title = element_text(size = 12), legend.text = element_text(size = 12), legend.position='bottom') +
         labs(title=sampleTitle, y='VAF',x='Coverage (log10)', plot.title=element_text(hjust=0.5))

    # CREATE PLOTS
    outputDir = paste0(runDir, "sampleQcDeamination/")
    if (!dir.exists(outputDir)){ dir.create(outputDir)}

    # PDF creation
    outputFile = paste0(outputDir, sampleId, '.sampleQcDeamination.pdf')
    print(paste0("writing output to pdf file: ", outputFile))
    pdf(file = outputFile, height = 7, width = 10)
    par(mar = c(1, 1, 1, 1))
    print(variantPlot)
    dev.off()

    # PNG creation
    outputPNG = paste0(outputDir, sampleId, '.sampleQcDeamination.png')
    print(paste0("writing output to png file: ", outputPNG))
    png(file = outputPNG, height = 25, width = 50, units="cm", res=1200, pointsize=1)
    par(mar = c(1, 1, 1, 1))
    print(variantPlot)
    dev.off()

}



#     title = textGrob(paste0(sampleId, ' ( TCP ', 100*(purity$purity), '% )' ), gp = gpar(fontface = "bold", fontsize = 16))
#
#     plotWidth=3
#     gapWidth=0
#     grid.arrange(plot_grid(NULL,NULL,NULL, title, NULL,NULL,NULL,NULL,
#                            variantPlot,
#                            NULL,NULL,NULL,NULL,NULL,NULL,NULL,
#                            ncol=1, nrow=3,
#                            rel_heights=c(1,20,1),align='v',axis='l'))
#
#     grid.arrange(plot_grid(NULL,NULL,NULL, title, NULL,NULL,NULL,NULL,
#                            variantPlot,NULL,exonMedianPlot,NULL,exonNormPlot,NULL,geneCnPlot,NULL,
#                            NULL,NULL,NULL,NULL,NULL,NULL,NULL,
#                            ncol=8, nrow=3, rel_widths=c(plotWidth,gapWidth,plotWidth,gapWidth,plotWidth,gapWidth,plotWidth,0.1),
#                            rel_heights=c(1,20,1),align='v',axis='l'))




#Samples =read.table("/home/evandijk/TEST_QC_Plots/QC_rundir/sampleMSI.txt")



#for(S in Samples$V1){
  runSampleQcDeamination()
#}


