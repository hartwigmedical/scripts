library(RMySQL)
library(dplyr)
library(tidyr)
library(GenomicRanges)

outputDir = "~/garvan/RData/"
resourceDir = "~/hmf/resources/"

outputDir = "~/Documents/LKCGP_projects/RData/"
resourceDir = "/Users/marwo2/Documents/LKCGP_projects/RData/Resources/"

referenceDir = paste0(outputDir, "reference/")
processedDir = paste0(outputDir, "processed/")

garvanPanel = c('AADACL4','ABL1','ABL2','ACPP','ACTB','ACVR1','ACVR2A','ACVR2B','ACVRL1','ADD3','AGBL4','AKT1','AKT2','AKT3','ALK','AMER1','ARAF','ARID1A',
  'ARID2','ASXL1','ASXL2','ATF7IP','ATM','ATR','ATRX','AURKA','AURKB','AURKC','AUTS2','AXL','BAZ1A','BCL11B','BCL2','BCOR','BCORL1','BCR','BLK',
  'BPIFB1','BRAF','BRCA1','BRCA2','BRD2','BRD3','BRD4','BTG1','BTK','C10orf11','MALRD1','CARD11','CBFB','CBL','CCND1','CCND2','CCND3','CCNE1','CD200',
  'CDK2','CDK4','CDK6','CDK7','CDK9','CDKN1A','CDKN1B','CDKN2A','CDKN2B','CDKN2B-AS1','CEBPA','CHD4','CHEK1','CNOT3','COL1A1','CREBBP','CRKL','CRLF2',
  'CROCC','CSF1R','CSF3R','CTCF','CTLA4','CTNNB1','DDR1','DDR2','DDX31','DDX3X','DGCR8','DHX15','DKK1','DLEU1','DLL1','DLL3','DLL4','DNM2','DNMT1','DNMT3A',
  'DNMT3B','DROSHA','DYM','EBF1','EED','EGFR','ELF1','ENG','EP300','EPHA2','EPHA3','ERBB2','ERBB3','ERBB4','ERG','ERRFI1','ETV6','EYS','EZH2','FAM49A',
  'FAT2','FBN2','FBXO11','FBXO28','FBXW7','FGF1','FGF2','FGF3','FGF4','FGFR1','FGFR2','FGFR3','FGFR4','FGR','FLG','FLT1','FLT3','FLT4','FMR1','FOXO1',
  'FYN','GAK','GALNTL6','GATA1','GATA2','GATA3','GFI1B','GLI2','GLIS2','GMDS','GNA11','GNAQ','GNB1','GRB2','H3F3A','HCK','HDAC2','HDAC7','HDAC9','HGF','HIST1H3B',
  'HOXA10','HRAS','HSP90AA1','HSP90AB1','HSP90B1','ID3','IDH1','IDH2','IGF1R','IGHJ4','IGHV3-71','IGLL5','IGLV4-69','IKZF1','IKZF2','IKZF3','IL1B',
  'IL2','IL6','IL6R','IL7R','INO80','INSR','JAK1','JAK2','JAK3','JMJD1C','KBTBD4','KDM1A','KDM6A','KDR','KIAA1549','KIT','KMT2A','KMT2C','KMT2D','KRAS',
  'LAPTM4B','LCK','LEF1','LEMD3','LINC00511','LMO2','LPAR6','LPHN2','LRIG1','LYN','MAP2K1','MAP2K2','MAPK1','MAPK11','MAPK14','MAPK3','MAX','MCL1',
  'MDM2','MED12','MEF2C','MEF2D','MERTK','MET','MGA','MLLT1','MLLT10','MST1R','MTAP','MTOR','MYB','MYC','MYCN','MYCNOS','MYD88','MYH11','NCOR1','NDE1',
  'NF1','NIPBL','NKX2-1','NOTCH1','NOTCH2','NOTCH3','NOTCH4','NPLOC4','NPM1','NRAS','NRG1','NRG2','NRG3','NRG4','NRXN1','NTRK1','NTRK2','NTRK3','NUP214',
  'NUP98','OSM','OSMR','OTX2','OTX2-AS1','P2RY8','PAX5','PCBP1','PDGFA','PDGFB','PDGFRA','PDGFRB','PHF6','PHLPP1','PHLPP2','PIK3CA','PIK3CB','PIK3CD',
  'PIK3CG','PIK3R1','PIPOX','PLK1','PRKCA','PRKCB','PRKCD','PRKCE','PRKCH','PRKCQ','PTCH1','PTCHD4','PTEN','PTK2','PTPN11','PTPN2','PTPRD','PVT1','RAF1',
  'RAG1','RAG2','RASA1','RB1','RET','RGPD2','RHEB','RHOA','ROS1','RPL10','RPL22','RPL5','RPS6KB1','RPTOR','RUNX1','RUNX1T1','SELP','SERP2','SETBP1','SETD2',
  'SF3B1','SH2B3','SHANK2','SHH','SI','SIRPA','SIX1','SIX2','SLC6A18','SMARCA4','SMARCB1','SMC3','SMG8','SMO','SMYD3','SNCAIP','SOS1','SRC','STAG2','STAT3',
  'STAT5B','SUFU','SUZ12','SYK','TAL1','TAL2','TBL1XR1','TBR1','TCF3','TCF7','TERT','TET2','TFAP4','TGFB1','TGFB2','TGFB3','TGFBR1','TGFBR2','TLX1','TLX3',
  'TOX','TP53','TRAJ29','TRAP1','TRDD2','TRDD3','TRDV2','TSC1','TSC2','TSPYL2','TTYH1','TYRO3','U2AF1','UBA2','USP22','USP7','USP9X','VEGFA','VEGFB','WAC',
  'WDR64','WEE1','WHSC1','WT1','XBP1','XPO1','ZBTB7A','ZCCHC7','ZEB2','ZFHX3','ZFP36L2','ZFPM2','ZIC1','ZMIZ1','ZMYM3','ZNF217','ZNF384')
garvanPanel = data.frame(gene_name = garvanPanel, stringsAsFactors = F)


########################### Cosmic Genes
cosmicCurated = read.csv(paste0(resourceDir, "CosmicCurated.csv"), stringsAsFactors = F)
cosmicCurated$cosmicCurated <- TRUE
cosmicCurated = cosmicCurated[, c("Genes","cosmicCurated")]
colnames(cosmicCurated) <- c("gene_name", "cosmicCurated")

cosmicCensus = read.csv(paste0(resourceDir, "CosmicCensus.csv"), stringsAsFactors = F)
colnames(cosmicCensus) <- c("gene_name", "cosmicOncogene", "cosmicTsg")
cosmicCensus = cosmicCensus %>% filter(cosmicOncogene | cosmicTsg)

cosmicGenes = merge(cosmicCurated, cosmicCensus, by = "gene_name", all = T)
rm(cosmicCurated, cosmicCensus)

########################### KNOWN AMPS AND DELS
cgi = read.table(paste0(resourceDir, "cgi_biomarkers_per_variant.tsv"), stringsAsFactors = F, header = T, sep = "\t") %>%
    select(Biomarker) %>%
    separate(Biomarker, c("gene","alteration"), sep = " ") %>%
    filter(alteration %in% c("amplification", "deletion")) %>%
    distinct(gene, alteration)

onco = read.table(paste0(resourceDir, "oncoKb.tsv"), stringsAsFactors = F, header = T, sep = "\t") %>%
    select(gene = Gene, alteration = Alteration) %>%
    mutate(alteration = tolower(alteration)) %>%
    filter(alteration %in% c("amplification", "deletion")) %>%
    distinct(gene, alteration)

civic = read.table(paste0(resourceDir, "civic_variants.tsv"), stringsAsFactors = T, header = T, sep = "\t", quote = "") %>%
    select(gene, alteration = variant) %>%
    mutate(alteration = tolower(alteration)) %>%
    filter(alteration %in% c("amplification", "deletion")) %>%
    distinct(gene, alteration)

manual = data.frame(gene = "ZNF703", alteration = "amplification")
knownAmpsDels = bind_rows(cgi, onco) %>%
    bind_rows(civic) %>%
    bind_rows(manual) %>%
    select(gene_name = gene, alteration) %>%
    distinct(gene_name, alteration) %>%
    mutate(value = T) %>% spread(alteration, value)
rm(cgi, onco, civic, manual)


########################### Fragile Sites
fragileSites = read.csv(paste0(resourceDir, "FragileSite.csv"))
fragileSites$range = GRanges(fragileSites$chrom, IRanges(fragileSites$start, fragileSites$end))

load(paste0(referenceDir, "canonicalTranscripts.RData"))
canonicalTranscripts$range = GRanges(canonicalTranscripts$chromosome, IRanges(canonicalTranscripts$geneStart, canonicalTranscripts$geneEnd))
ol = as.matrix(findOverlaps(fragileSites$range, canonicalTranscripts$range))
fragileGenes = canonicalTranscripts[ol[,2], c("gene", "chromosome")]
fragileGenes$fragile <- TRUE
fragileGenes$chromosome <- NULL
colnames(fragileGenes) <- c("gene_name", "fragile")
rm(canonicalTranscripts, ol, fragileSites)
save(fragileGenes, file = paste0(processedDir,"fragileGenes.RData"))

########################### Gene Panel
load(file = paste0(processedDir,"HmfRefCDSCv.RData"))

sig = 0.05
hmfSignificant =  HmfRefCDSCv %>% filter(qglobal_cv < sig) %>% distinct(gene_name)
hmfSignificant$dnds <- TRUE

garvanPanel$garvan <- TRUE

genePanel = merge(garvanPanel, hmfSignificant, by = "gene_name", all = T)
genePanel = merge(genePanel, cosmicGenes, by = "gene_name", all = T)
genePanel = merge(genePanel, knownAmpsDels, by = "gene_name", all=T)
#genePanel = genePanel %>% filter(!gene_name %in% c("POM121L12","TRIM49B","LPCAT2"))
save(fragileGenes, file = paste0(processedDir,"genePanel.RData"))

oncoGenes = genePanel %>% filter(cosmicOncogene | is.na(cosmicTsg)) %>% mutate(classification = "onco")
tsGenes = genePanel %>% filter(cosmicTsg | is.na(cosmicOncogene)) %>% mutate(classification = "tsg")
save(tsGenes, oncoGenes, file = paste0(processedDir,"driverGenes.RData"))


