library(tidyr)
library(dplyr)
library(stringi)
# Base panel
bachelorGenePanel = read.csv('~/Dropbox/HMF Australia team folder/Bachelor/genePanel152.csv')

# Bachelor SNV / INDEL output
bachelorLatest = read.csv('~/Dropbox/HMF Australia team folder/Bachelor/bach_filtered_clinvar152.csv')
#bachelorLatest = read.csv('~/Dropbox/HMF Australia team folder/Bachelor/bach_germline_variants_drivers_152.csv')

#'SELECT * FROM somaticVariant WHERE worstCodingEffect <> 'NONE' AND gene IN ('ALK','APC','ATR','AXIN2','BAP1','BMPR1A','BRAF','BRCA1','CBL','CDC73','CDH1','CDK4','CDKN1C','CDKN2A','CEBPA','CHEK2','CTR9','CYLD','DDB2','DICER1','DROSHA','EGFR','ELANE','EPCAM','ETV6','EXT1','EXT2','FANCM','FLCN','GATA2','GJB2','HMBS','HRAS','JMJD1C','KIT','KRAS','LMO1','MAP2K1','MAP2K2','MAX','MEN1','MET','MITF','MPL','MTAP','NF1','NF2','NRAS','PAX5','PDGFRA','PHOX2B','POLD1','POT1','PRDM9','PRKAR1A','PRSS1','PTCH1','PTEN','PTPN11','RAD51D','RAF1','RB1','RECQL','RET','RHBDF2','RUNX1','SDHAF2','SDHB','SDHC','SDHD','SETBP1','SHOC2','SMAD4','SMARCA4','SMARCB1','SMARCE1','SOS1','STAT3','STK11','SUFU','TGFBR1','TMEM127','TNFRSF6','TP53','TSC1','TSC2','TSHR','VHL','WT1','ERCC1','FANCD2','FANCE','FANCF','FANCI','FANCL','HNF1A','NHP2','NOP10','NTHL1','PMS1','PRF1','SH2B3','SPRTN','ABCB11','BLM','BUB1B','DIS3L2','DOCK8','ERCC2','ERCC3','ERCC4','ERCC5','FAH','FANCA','FANCC','FANCG','GBA','HFE','ITK','MUTYH','NBN','POLH','RECQL4','RMRP','SBDS','SERPINA1','SLC25A13','TRIM37 ','WRN','XPA','XPC','ATM','BRCA2','BRIP1','CDKN1B','COL7A1','FH','MLH1','MSH2','MSH6','PALB2','PMS2','POLE','RAD51C','SDHA','TERT','UROD','DKC1','GPC3','SH2D1A','WAS','SRY');
bachelorSomatic = read.csv('~/Dropbox/HMF Australia team folder/Bachelor/somaticVariants152.csv')  ### please replace with drivers from driver catalog.

### SQL QUERY FOR GERMLINE Copy NUMBER EVENTS  #### Please replace with direct query from R based on genePanel152 file
bachelorGermlineCN= read.csv('~/Dropbox/HMF Australia team folder/Bachelor/germlineHettoHom152.csv') 
# SELECT g.* FROM geneCopyNumber g INNER JOIN purity p ON p.sampleId = g.sampleId WHERE qcstatus = 'PASS' AND STATUS <> 'NO_TUMOR' AND purity > 0.195 AND (germlineHetRegions>0 OR germlineHomRegions>0) 
#AND gene IN ('ALK','APC','ATR','AXIN2','BAP1','BMPR1A','BRAF','BRCA1','CBL','CDC73','CDH1','CDK4','CDKN1C','CDKN2A','CEBPA','CHEK2','CTR9','CYLD','DDB2','DICER1','DROSHA','EGFR','ELANE','EPCAM','ETV6',
#'EXT1','EXT2','FANCM','FLCN','GATA2','GJB2','HMBS','HRAS','JMJD1C','KIT','KRAS','LMO1','MAP2K1','MAP2K2','MAX','MEN1','MET','MITF','MPL','MTAP','NF1','NF2','NRAS','PAX5','PDGFRA','PHOX2B','POLD1','POT1',
#PRDM9','PRKAR1A','PRSS1','PTCH1','PTEN','PTPN11','RAD51D','RAF1','RB1','RECQL','RET','RHBDF2','RUNX1','SDHAF2','SDHB','SDHC','SDHD','SETBP1','SHOC2','SMAD4','SMARCA4','SMARCB1','SMARCE1','SOS1','STAT3',
#STK11','SUFU','TGFBR1','TMEM127','TNFRSF6','TP53','TSC1','TSC2','TSHR','VHL','WT1','ERCC1','FANCD2','FANCE','FANCF','FANCI','FANCL','HNF1A','NHP2','NOP10','NTHL1','PMS1','PRF1','SH2B3','SPRTN','ABCB11',
#'BLM','BUB1B','DIS3L2','DOCK8','ERCC2','ERCC3','ERCC4','ERCC5','FAH','FANCA','FANCC','FANCG','GBA','HFE','ITK','MUTYH','NBN','POLH','RECQL4','RMRP','SBDS','SERPINA1','SLC25A13','TRIM37 ','WRN','XPA','XPC',
#'ATM','BRCA2','BRIP1','CDKN1B','COL7A1','FH','MLH1','MSH2','MSH6','PALB2','PMS2','POLE','RAD51C','SDHA','TERT','UROD','DKC1','GPC3','SH2D1A','WAS','SRY');

#info only please ignore
bachelorTranscripts= read.csv('~/Dropbox/HMF Australia team folder/Bachelor/trancripts152.csv')   
load('~/hmf/RData/Processed/highestPurityCohortSummary.RData')

includeList = c('BAP1','CDKN2A','ATM','CHEK2','APC','BMPR1A','BRCA1','BRCA2','MEN1','MLH1','MSH2','MSH6','MUTYH','NF2','PMS2','PTEN','RB1','RET','SDHAF2','SDHB','SDHC','SDHD','SMAD4','STK11','TP53','TSC1','TSC2','VHL','WT1')
#View(includeList)

# Germline CN deletions
View(bachelorGermlineCN)# %>% filter((gene %in% includeList)))

filteredBachelor = merge(bachelorLatest %>% mutate(refSampleVaf=GermlineAltCount/GermlineReadDepth) %>% 
  group_by(SampleId,Gene) %>% mutate(netFS=sum(stri_length(Ref)-stri_length(Alt)) %% 3) %>% ungroup() %>%  #net frameshift per sample
  group_by(Gene,Chromosome,Position,Ref,Alt) %>% mutate(medianRefSampleVaf=median(refSampleVaf),medianNetFS=median(netFS),medianRefReadDepth=median(GermlineReadDepth)) %>% ungroup() %>%  #median vaf per sample
  filter(!(ClinvarSignificance %in% c('Benign/Likely_benign','Benign','Likely_benign')),   #Known Benign NONSENSE, FRAMESHIFT or SPLICE
  !Effects=='missense_variant'|!(ClinvarDiagnosis %in% c('PI_Z|PI_Z(AUGSBURG)|PI_Z(TUN)|Inborn_genetic_diseases|Alpha-1-antitrypsin_deficiency|FRAXE|not_provided','PI_M(PROCIDA)|Alpha-1-antitrypsin_deficiency','Hereditary_pancreatitis|not_specified')),  #Non cancer related pathogenic
  !(Gene=='GJB2' & HgvsProtein %in% c('p.Leu90Pro')),   #Associated with deafness. Unlikely to be pathogenic
  !(Gene=='TSHR' & HgvsProtein %in% c('p.Pro68Ser')),   #Associated with Hypothyroidism, congenital, nongoitrous. Unlikely to be pathogenic in cancer
  WorstCodingEffect!='SPLICE'|Type!='INDEL',  #remove microsatellite expansions in INDELs adjacent to splice acceptor sites
  Filter!='ARTEFACT',  #remove samples which are found by bachelor to have a negative implied copy number in the tumor
  medianRefSampleVaf>0.2&medianRefSampleVaf<0.8,  #remove variants consistently found at very low or very high Vaf in ref
  !grepl('frameshift',Effects)|medianNetFS!=0),   #remove frameshifts which are offset by another frameshift in 50% or more of samples of which is found
  highestPurityCohortSummary %>% select(sampleId,primaryTumorLocation),by.x='SampleId',by.y='sampleId',all.x=F) %>% 
  mutate(include=Gene %in% includeList)

### Match with somatics - can just do this in the catalog
filteredBachelor=merge(filteredBachelor,bachelorSomatic %>% filter(filter=='PASS',canonicalCodingEffect!='SYNONYMOUS') %>% 
             group_by(sampleId,gene) %>% summarise(somaticCount=n(),somaticPloidy=sum(adjustedVaf*adjustedCopyNumber),somaticMinCopyNumber=min(adjustedCopyNumber)),by.x=c('Gene','SampleId'),by.y=c('gene','sampleId'),all.x=T) %>% mutate(somaticHit=Biallelic=='true'|!is.na(somaticCount))

View(filteredBachelor)
# By  Biallleic
View(filteredBachelor %>%
             group_by(include,somaticHit) %>% count() %>% 
             spread(somaticHit,n,fill=0) %>% mutate(total=`TRUE`+`FALSE`,LOHProportion=round(`TRUE`/total,3)) %>%
             arrange(-total))

# By Biallleic/LOH
View(filteredBachelor %>% mutate(LOH=AdjustedVaf*AdjCopyNumber<0.5) %>% unite(LOH_BI,LOH,somaticHit) %>%
       group_by(include,LOH_BI) %>% count() %>% 
       spread(LOH_BI,n,fill=0))#

# By Gene & Biallleic
View(merge(filteredBachelor %>%
       group_by(include,Gene,somaticHit) %>% count() %>% 
       spread(somaticHit,n,fill=0) %>% mutate(total=`TRUE`+`FALSE`,LOHProportion=round(`TRUE`/total,2)) %>%
       arrange(include,-total),bachelorTranscripts %>% select(gene,chromosome,codingStart,codingBases,exonBases),by.x='Gene',by.y='gene'))

# By Gene & Biallleic/LOH
View(filteredBachelor %>% mutate(LOH=AdjustedVaf*AdjCopyNumber<0.5) %>% unite(LOH_BI,LOH,somaticHit) %>%
       group_by(include,Gene,LOH_BI) %>% count() %>% 
       spread(LOH_BI,n,fill=0))# 

# By Gene & PrimaryTumorLocation
View(filteredBachelor %>%
       group_by(Biallelic,Gene,primaryTumorLocation) %>% count() %>% 
       arrange(-n)  %>% spread(primaryTumorLocation,n))

temp = HmfRefCDSCv %>% filter(cancerType=='All')

temp = temp %>% filter(wmis_cv>0,wnon_cv>0) %>% 
  mutate(d_mis=n_mis*pmin(1,(wmis_cv-1)/wmis_cv)) %>%
  mutate(d_ind=n_ind*pmin(1,(wind_cv-1)/wind_cv)) %>%
  mutate(d_spl=n_spl*pmin(1,(wspl_cv-1)/wspl_cv)) %>%
  mutate(d_non=n_non*pmin(1,(wnon_cv-1)/wnon_cv))
temp %>% filter(!is.na(d_mis)) %>% summarise(sum(d_mis),sum(n_mis),sum(d_non),sum(n_non),sum(d_spl),sum(n_spl))


View(temp %>% select(gene_name,d_mis,n_mis))
