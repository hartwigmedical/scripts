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

runSampleQcReport<-function()
{



#     Parse and check inputs
    args <- commandArgs(trailingOnly = TRUE)
    print(args)

    if(length(args) < 3)
    {
      print("Requires arguments: 1=SampleId, 2=DataDir, 3=RunDir, [4=ResourceDir (default=/data/resources/public)], [5=ExcludedExonsFile (default=/data/resources/ops/panel/excludedExons.tsv)]")
      stop()
    }

    sampleId <- args[1]
    sampleDataDir <- args[2]
    runDir = args[3]

    # Optional arguments with defaults
    resourceDir = ifelse(length(args) >= 4, args[4], "/data/resources/public")
    excludedFile = ifelse(length(args) >= 5, args[5], "/data/resources/ops/panel/excludedExons.tsv")

    cohortMedianDepthFile = paste0(resourceDir, "/target_regions/38/target_regions_exon_relative_depth.38.csv")
    cobaltRegionsFile = paste0(resourceDir, "/target_regions/38/target_regions_normalisation.38.tsv")
    ensembl_gene_data = paste0(resourceDir, "/ensembl_data_cache/38/ensembl_gene_data.csv")
    driverGenePanel = paste0(resourceDir, "/gene_panel/38/DriverGenePanel.38.tsv")

    print(sampleId)

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
    check_file(excludedFile)
    print("files present")

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
    excludedExons = load_file(excludedFile,'excluded exons')




    # sort ensembl on chromosomal order
    ensemblGenes = ensemblGenes %>% mutate(ChromIndex=as.integer(ensemblGenes$Chromosome))
    ensemblGenes$ChromIndex[which(ensemblGenes$Chromosome=='X')] <- 23
    ensemblGenes$ChromIndex[which(ensemblGenes$Chromosome=='Y')] <- 24
    ensemblGenes = ensemblGenes %>% arrange(ChromIndex, GeneStart)

    colnames(ensemblGenes)[2] <- 'gene'
    driverGenes = dplyr::inner_join(ensemblGenes, driverGenePanel, by="gene")

    # filter genes of interest for OncoPanel. Check: "Proposed validation gene list" in Oncopanel > CNV validation
    # 15 genes for DELs
    reportGenesDeletions = c('PALB2', 'RAD51B', 'RAD51C', 'BRCA1', 'BRCA2', 'CDKN2A', 'TP53', 'CREBBP', 'EP300', 'ARID1A', 'PTEN', 'RB1','MTAP','KEAP1','SMAD4')
    # 17 genes for AMPS
    reportGenesAmplification = c('AR','CCNE1','EGFR','ERBB2','FGFR2','FLT4','KDR','KIT','MDM2','MET','MYC','PDGFRA','PDGFRB','PIK3CA','RAF1','ROS1','FGFR1')

    genesOfInterest = c(reportGenesDeletions, reportGenesAmplification)

    genesList = driverGenes %>% filter(gene %in% genesOfInterest)
    likelihoodType = as.data.frame(cbind(genesList$gene, genesList$likelihoodType))
    colnames(likelihoodType) = c('Gene', 'likelihoodType')
    # TODO: color genes in legend by likelihoodType (TSG/ONCO)

    genesList = as.data.frame(genesList$gene)
    colnames(genesList) = c('Gene')
    genesList = genesList %>% mutate(GeneIndex=row_number())


    exons = load_file(sampleExonCoverageFile,'sample exon coverage')

    colnames(exons) = c('Gene','Chromosome','PosStart','PosEnd','ExonRank','MedianDepth')
    print(paste(nrow(exons)," total exons in driver catalog"),sep="")

   # genesNotInDesign = c("FANCM","H3-3B","H3C13","LINC00290","LINC01001","OR11H1","OR4F21","OR4N2","RABAC1","SPATA31A7","U2AF1","SMARCE1") 
    genesNotInDesign = c("SPATA31A7","LINC01001","U2AF1")
    exons = exons %>% filter(!(Gene %in% genesNotInDesign))

    outputDir = paste0(runDir, "sampleQcReports/")
        if (!dir.exists(outputDir)){ dir.create(outputDir)}

    exons = merge(exons,genesList,by='Gene',all.x=T)
    exons = merge(exons,cohortMedianDepth,by=c('Gene','ExonRank'),all.x=T)
    exons = exons %>% mutate(MedianSampleDepth=median(MedianDepth),
                             RelDepth=MedianDepth/MedianSampleDepth,
                             NormDepth=ifelse(CohortRelDepth>0,MedianDepth/CohortRelDepth,MedianDepth))
    exons = exons %>% mutate(Position=floor(PosEnd/1000)*1000+1)
    exons = merge(exons,cobaltRegions,by=c('Chromosome','Position'),all.x=T)
    allExons = exons
    exons = exons %>% filter(Gene %in% genesOfInterest)

    allExons$excludedInValidation = rep(FALSE,nrow(allExons))
    for(i in 1:nrow(allExons)){
       if(!(allExons$Gene[i] %in% excludedExons$Gene)){
            allExons$excludedInValidation[i]=TRUE
       }
       else{
           excludedExonsGene = excludedExons[excludedExons$Gene == allExons$Gene[i],]$Trancript.ID.exons.not.included
           excludedExonsGene = strsplit(excludedExonsGene,split=",")[[1]]
       

           if(allExons$ExonRank[i] %in% excludedExonsGene){
               allExons$excludedInValidation[i]=TRUE
           }
       }
    }
    excluded=allExons[allExons$excludedInValidation==TRUE,]
    print(excluded[order(excluded$Gene),])
    checkedExons = allExons[allExons$excludedInValidation==FALSE,]
    print(paste(nrow(checkedExons)," exons checked"),sep="")

    checkedExons$excludedInValidation=NULL
    insufficientCoverageExons = checkedExons %>% filter(MedianDepth < 100)
    
    write.table( insufficientCoverageExons %>% filter(MedianDepth < 100),paste0(outputDir,sampleId,'.insufficientCoverage.tsv'),row.names=F,quote=F,sep="\t")
    entriesPerPage = 30



    print(paste(nrow( checkedExons %>% filter(MedianDepth < 100))," exons with insufficient coverage",sep=""))
    

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

    # extract AF
    genotypeData = geno(somaticVcf)
    afData = genotypeData$AF
    colnames(afData) = c('AF')
    rownames(afData) = NULL

    varData = cbind(somaticVariants,impactData %>% dplyr::select(Gene,CanonicalCodingEffect,Hotspot))
    varData = cbind(varData,afData)

    varGeneData = varData %>% filter(Gene %in% genesOfInterest)
    varGeneData = varData %>% filter(CanonicalCodingEffect!='NONE')

    varGeneData = merge(varGeneData,genesList,by='Gene',all.y=T)
    varGeneData = varGeneData %>% mutate(CanonicalCodingEffect=ifelse(is.na(CanonicalCodingEffect),'EMPTY',as.character(CanonicalCodingEffect)),
                                         IsHotspot=ifelse(is.na(Hotspot)|CanonicalCodingEffect=='EMPTY','EMPTY',ifelse(Hotspot,'HOTSPOT','NON-HOTSPOT')),
                                         AF=ifelse(is.na(AF),-1,AF))
    likelihoodType = likelihoodType %>% mutate(likelihoodTypeColours=ifelse(likelihoodType=='ONCO','tomato1','turquoise3'))
    likelihoodTypeColours = likelihoodType$likelihoodTypeColours

    legendCodingEffectColours =c('red','blue')
    codingEffectColours =c('yellow','red','blue')

    # plot to generate a legend for the variant types
    varCreateLegendPlot = ggplot(varGeneData, aes(x=reorder(Gene,-GeneIndex),y=AF)) +
      geom_point(aes(shape=CanonicalCodingEffect,color=IsHotspot)) +
      scale_color_manual(values = codingEffectColours) +
      theme(legend.title = element_blank()) +
      coord_flip() +
      theme(legend.position='bottom',legend.text = element_text(size = 12))

    varLegendPlot = as_ggplot(get_legend(varCreateLegendPlot))

    variantPlot = ggplot(varGeneData, aes(x=reorder(Gene,-GeneIndex),y=AF)) +
      geom_point(aes(shape=CanonicalCodingEffect,color=IsHotspot)) +
      theme(axis.ticks.y = element_blank()) +
      scale_y_continuous(limits = c(0, 1.0)) +
      theme(axis.text.y=element_text(colour=rev(likelihoodTypeColours))) +
      coord_flip() +
      scale_color_manual(values = codingEffectColours) +
      theme(panel.grid.minor.x = element_blank(),text = element_text(size = 12),axis.text = element_text(size = 12)) +
      theme(legend.position='none') +
      labs(title='Variant VAFs by Type', x='',y='')


    # CREATE PLOTS

    # PDF creation
    outputFile = paste0(outputDir, sampleId, '.sampleQcReport.pdf')
    print(paste0("writing output to pdf file: ", outputFile))
    pdf(file = outputFile, height = 14, width = 20)
    par(mar = c(1, 1, 1, 1))
    title = textGrob(paste0(sampleId, ' ( mTCP ', 100*(purity$purity), '%; Percentage_exon_100x ', format(100-100*(nrow(insufficientCoverageExons)/nrow(checkedExons)),digits=4), '% )'), gp = gpar(fontface = "bold", fontsize = 16))
    plotWidth=3
    gapWidth=0
    grid.arrange(plot_grid(NULL,NULL,NULL, title, NULL,NULL,NULL,NULL,
                           variantPlot,NULL,exonMedianPlot,NULL,exonNormPlot,NULL,geneCnPlot,NULL,
                           NULL,varLegendPlot,NULL,NULL,NULL,NULL,NULL,
                           ncol=8, nrow=3, rel_widths=c(plotWidth,gapWidth,plotWidth,gapWidth,plotWidth,gapWidth,plotWidth,0.1),
                           rel_heights=c(1,20,1),align='v',axis='l'))
    print(insufficientCoverageExons)
    insufficientCoverageExons$GeneIndex =NULL
    insufficientCoverageExons$CohortRelDepth =NULL
    insufficientCoverageExons$RelDepth =NULL
    insufficientCoverageExons$NormDepth =NULL
    insufficientCoverageExons$RelativeEnrichment =NULL
    
    if(nrow(insufficientCoverageExons)>0){
      rownames(insufficientCoverageExons) = 1:nrow(insufficientCoverageExons)
      grid.newpage()
      
      for(i in seq(1,nrow(insufficientCoverageExons),entriesPerPage)){
        grid.table(insufficientCoverageExons[i:min(nrow(insufficientCoverageExons),(i+entriesPerPage-1)),])
        if( i+entriesPerPage <= nrow(insufficientCoverageExons)){
          grid.newpage()
        }
      }
    }

    dev.off()

    # PNG creation
    outputPNG = paste0(outputDir, sampleId, '.sampleQcReport.png')
    print(paste0("writing output to png file: ", outputPNG))
    png(file = outputPNG, height = 25, width = 50, units="cm", res=1200, pointsize=1)
    par(mar = c(1, 1, 1, 1))
    title = textGrob(paste0(sampleId, ' ( mTCP ', 100*(purity$purity), '%; Percentage_exon_100x ', format(100-100*(nrow(insufficientCoverageExons)/nrow(checkedExons)),digits=4), '% )'), gp = gpar(fontface = "bold", fontsize = 16))
    plotWidth=3
    gapWidth=0
    grid.arrange(plot_grid(NULL,NULL,NULL, title, NULL,NULL,NULL,NULL,
                           variantPlot,NULL,exonMedianPlot,NULL,exonNormPlot,NULL,geneCnPlot,NULL,
                           NULL,varLegendPlot,NULL,NULL,NULL,NULL,NULL,
                           ncol=8, nrow=3, rel_widths=c(plotWidth,gapWidth,plotWidth,gapWidth,plotWidth,gapWidth,plotWidth,0.1),
                           rel_heights=c(1,20,1),align='v',axis='l'))
    dev.off()

}


runSampleQcReport()

