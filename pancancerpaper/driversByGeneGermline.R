library(tidyr)
library(dplyr)
library(purple)

#### COMMON DATA
path = "~/hmf/resources/Bachelor/"
load(file = "~/hmf/RData/processed/driverGenes.RData")
includeList = c('BAP1','CDKN2A','ATM','CHEK2','APC','BMPR1A','BRCA1','BRCA2','MEN1','MLH1','MSH2','MSH6','MUTYH','NF2','PMS2','PTEN','RB1','RET','SDHAF2','SDHB','SDHC','SDHD','SMAD4','STK11','TP53','TSC1','TSC2','VHL','WT1')

genomic_coordinates <- function(chromosome, position, ref, alt) {
  return (paste0(chromosome, ":", position, ref, ">", alt))
}


##### GERMLINE CNA
load(file = "~/hmf/RData/Processed/fragileGenes.RData")
load(file = "~/hmf/RData/reference/hpcGeneCopyNumberDeletes.RData")
bachelorGenes = read.csv(paste0(path,'genePanel152.csv')) %>% select(Gene)
germlineDeletes = hpcGeneCopyNumberDeletes %>% 
  filter(germlineHetRegions > 0 | germlineHomRegions > 0) %>% 
  filter(gene %in% bachelorGenes$Gene) %>% left_join(fragileGenes, by = c("gene" = "gene_name")) %>%
  mutate(category = "CNA", fragile = ifelse(is.na(fragile),F, T), driver = "Del", partial = somaticRegions > 1) %>%
  select(sampleId, gene,  driver, category)
rm(bachelorGenes, fragileGenes, hpcGeneCopyNumberDeletes)

###### GERMLINE SOMATICS
bachelorLatest = read.csv(paste0(path,'bach_filtered_clinvar152.csv'), stringsAsFactors = F) 
load(file = "~/hmf/RData/Processed/hpcDndsMutations.RData")
hpcMutations = hpcMutations %>% filter(impact != 'Synonymous')
somaticHits = hpcMutations %>% group_by(sampleId, gene) %>% summarise(somaticVariants = n())

germlineSomatics = bachelorLatest %>% mutate(refSampleVaf=GermlineAltCount/GermlineReadDepth) %>%
  group_by(SampleId,Gene) %>% mutate(netFS=sum(nchar(Ref)-nchar(Alt)) %% 3) %>% ungroup() %>%  #net frameshift per sample
  group_by(Gene,Chromosome,Position,Ref,Alt) %>% mutate(medianRefSampleVaf=median(refSampleVaf),medianNetFS=median(netFS),medianRefReadDepth=median(GermlineReadDepth)) %>% ungroup() %>%  #median vaf per sample
  filter(
    !(ClinvarSignificance %in% c('Benign/Likely_benign','Benign','Likely_benign')),   #Known Benign NONSENSE, FRAMESHIFT or SPLICE
    !Effects=='missense_variant'|!(ClinvarDiagnosis %in% c('PI_Z|PI_Z(AUGSBURG)|PI_Z(TUN)|Inborn_genetic_diseases|Alpha-1-antitrypsin_deficiency|FRAXE|not_provided','PI_M(PROCIDA)|Alpha-1-antitrypsin_deficiency','Hereditary_pancreatitis|not_specified')),  #Non cancer related pathogenic
    !(Gene=='GJB2' & HgvsProtein %in% c('p.Leu90Pro')),   #Associated with deafness. Unlikely to be pathogenic
    !(Gene=='TSHR' & HgvsProtein %in% c('p.Pro68Ser')),   #Associated with Hypothyroidism, congenital, nongoitrous. Unlikely to be pathogenic in cancer
    WorstCodingEffect!='SPLICE'|Type!='INDEL',  #remove microsatellite expansions in INDELs adjacent to splice acceptor sites
    Filter!='ARTEFACT',  #remove samples which are found by bachelor to have a negative implied copy number in the tumor
    medianRefSampleVaf>0.2&medianRefSampleVaf<0.8,  #remove variants consistently found at very low or very high Vaf in ref
    !grepl('frameshift',Effects)|medianNetFS!=0)   #remove frameshifts which are offset by another frameshift in 50% or more of samples of which is found
   
germlineSomatics = germlineSomatics %>%
  mutate(
    coordinate = genomic_coordinates(Chromosome, Position, Ref, Alt), 
    category = ifelse(Type == 'SNP', "SNV", Type),
    ploidy = AdjCopyNumber * AdjustedVaf, 
    variantLostInTumor = ploidy < 0.5,
    driver = ifelse(category == 'INDEL' & WorstCodingEffect == 'NONSENSE_OR_FRAMESHIFT', "Frameshift", WorstCodingEffect),
    driver = ifelse(category == 'INDEL' & WorstCodingEffect == 'MISSENSE', "Inframe", driver),
    driver = ifelse(category == 'SNV' & WorstCodingEffect == 'MISSENSE', "Missense", driver),
    driver = ifelse(category == 'SNV' & WorstCodingEffect == 'SPLICE', "Splice", driver),
    driver = ifelse(category == 'SNV' & WorstCodingEffect == 'NONSENSE_OR_FRAMESHIFT', "Nonsense", driver),
    driver = ifelse(category == 'SNV' & WorstCodingEffect == 'SYNONYMOUS', "Synonymous", driver),
    impact = driver,
    gene = Gene,
    biallelic = ifelse(Biallelic == "false", F, T)
    ) %>%
  select(sampleId = SampleId, coordinate, gene, category, driver, impact, pHGVS = HgvsProtein, germlineGenotype = GermlineStatus, biallelic, variantLostInTumor ) 

germlineSomatics = germlineSomatics %>% group_by(sampleId, gene) %>%
  summarise(
    n = n(), 
    driver = ifelse(n > 1, "Multihit", driver), 
    category = paste0(impact, collapse =","),
    coordinate = paste0(coordinate, collapse =","),
    pHGVS = paste0(pHGVS, collapse =","), 
    germlineGenotype = paste0(germlineGenotype, collapse =","), 
    biallellicInTumor = any(biallelic),
    variantLostInTumor = all(variantLostInTumor),
    germlineVariants = n) %>%
  select(-n)

germlineSomatics = germlineSomatics %>% 
  left_join(somaticHits, by = c("sampleId", "gene")) %>% 
  mutate(somaticVariants = ifelse(is.na(somaticVariants), 0, somaticVariants), wildTypeLostInTumor = somaticVariants > 0 | biallellicInTumor)

load(file = "~/hmf/RData/Reference/highestPurityCohort.RData")
germlineDriverCatalog = bind_rows(germlineSomatics, germlineDeletes) %>% 
  mutate(highConfidenceGenePanel = gene %in% includeList) %>%
  left_join(highestPurityCohort %>% select(sampleId, cancerType), by = "sampleId") %>%
  select(sampleId, cancerType, gene, coordinate, category, mutationClass = driver, pHGVS, biallellicInTumor, variantLostInTumor, wildTypeLostInTumor, germlineVariants, somaticVariants, highConfidenceGenePanel) 
save(germlineDriverCatalog, file = "~/hmf/RData/Processed/germlineCatalog.RData")

sampleIdMap = read.csv(file = "/Users/jon/hmf/secure/SampleIdMap.csv", stringsAsFactors = F)
germlineDriverCatalogSupp = germlineDriverCatalog %>% ungroup() %>%
  mutate(
    pHGVS = ifelse(pHGVS == ",","",pHGVS)) %>% 
  left_join(sampleIdMap, by = "sampleId") %>%
  select(-sampleId) %>%
  select(sampleId = hmfSampleId, everything()) 

write.csv(germlineDriverCatalogSupp, file = "~/hmf/RData/Supp/Supplementary Table 6_GermlinePredispositionVariants.csv", row.names = F) 
