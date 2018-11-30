library(stringi)
bachelorLatest = read.csv('~/Dropbox/HMF Australia team folder/Bachelor/bach_filtered_result_152.csv')
load('~/hmf/RData/Processed/highestPurityCohortSummary.RData')

558/

filteredBachelor = merge(bachelorLatest %>% mutate(refSampleVaf=RefCount/RefReadDepth) %>% 
                           group_by(SampleId,Gene) %>% mutate(netFS=sum(stri_length(Ref)-stri_length(Alt)) %% 3) %>% ungroup() %>%  #check Frameshift per sample
                           group_by(Gene,Chromosome,Position,Ref,Alt) %>% mutate(medianRefSampleVaf=median(refSampleVaf),medianNetFS=median(netFS)) %>% ungroup() %>%  #check median vaf
                           filter(!(ClinvarSignificance %in% c('Benign/Likely_benign','Benign','Likely_benign')),   #Known Benign NONSENSE, FRAMESHIFT or SPLICE
                                  !Effects=='missense_variant'|!(ClinvarDiagnosis %in% c('PI_Z|PI_Z(AUGSBURG)|PI_Z(TUN)|Inborn_genetic_diseases|Alpha-1-antitrypsin_deficiency|FRAXE|not_provided','PI_M(PROCIDA)|Alpha-1-antitrypsin_deficiency')),  #Non cancer related pathogenic
                                  !(Gene=='GJB2' & HgvsProtein %in% c('p.Leu90Pro')),   #Associated with deafness. Unlikely to be pathogenic
                                  !(Gene=='TSHR' & HgvsProtein %in% c('p.Pro68Ser')),   #Associated with Hypothyroidism, congenital, nongoitrous. Unlikely to be pathogenic in cancer
                                  !DbsnpId=='rs111033566',   #very common SNP in PRSS1 (47% of ExAc). unlikely to be pathogenic
                                  WorstCodingEffect!='SPLICE'|Type!='INDEL',  #remove microsatellite expansions in INDELs adjacent to splice acceptor sites
                                  Filter!='ARTEFACT',  #remove samples which are found by bachelor to have a negative implied copy number in the tumor
                                  medianRefSampleVaf>0.2,medianRefSampleVaf<0.8,  #remove variants consistently found at very low or very high Vaf in ref
                                  !grepl('frameshift',Effects)|medianNetFS!=0),   #remove frameshifts which are offset by another frameshift in 50% or more of samples of which is found
                         highestPurityCohortSummary %>% select(sampleId,primaryTumorLocation),by.x='SampleId',by.y='sampleId',all.x=F) 

View(filteredBachelor)
write.csv(filteredBachelor,'~/filteredBachelor.csv')

#Biallelic 
View(filteredBachelor %>%
       group_by(Biallelic) %>% count())

# By Variant & Biallleic
View(filteredBachelor %>%
       group_by(Gene,CodonInfo,ClinvarDiagnosis,ClinvarSignificance,Chromosome,Position,Effects,Ref,Alt,HgvsProtein,DbsnpId,CosmicId,medianRefSampleVaf,Biallelic) %>% count() %>% 
       arrange(-n)  %>% spread(Biallelic,n))

# By Gene & Biallleic
View(filteredBachelor %>%
       group_by(Gene,Biallelic) %>% count() %>% 
       spread(Biallelic,n) %>% mutate(total=true+false,LOHProportion=round(true/total,2)) %>%
       arrange(-total))

# By Gene & PrimaryTumorLocation
View(filteredBachelor %>%
       group_by(Gene,primaryTumorLocation) %>% count() %>% 
       arrange(-n)  %>% spread(primaryTumorLocation,n))

View(filteredBachelor %>% group_by(SampleId) %>% count())
