library(stringi)
bachelorLatest = read.csv('~/Dropbox/HMF Australia team folder/Bachelor/bach_filtered_result_152.csv')
load('~/hmf/RData/Processed/highestPurityCohortSummary.RData')

### SUSPICIOUS BUT SUGGEST TO KEEP #####
#   GJB2 p.Gly12fs rs80338939 (33) - present in 0.9% of germline, gene seems to be mainly associated with hearing.  BUT frequently biallleic  => Keep?
#   BLM  15:91306246 & 15:91306246 (STOP - 3 instances each)
#   COL7A1 3:48621017 - Enriched in Prostate cancer
#   DOCK8 9:271626 Splice Acceptor (5 instances) - but 2 biallelic

#### PROBABLY FILTER ############
#   RECQL 12:21623219 (STOP GAINED - 5 instances) - Near end of gene (620/649 codon) never bialllelic   => probably filter?
#   NBN 8:90955525 (STOP GAINED - 4 instances) - Near end of gene (714/754 codon) only once bialllelic => probably filter?
#   KRAS  12:25368442 (STOP GAINED 2 instances) - (codon 168/189) => probably filter?
#   TSHR p.Pro68Ser (4 instances) - assocaited with Hypothyroidism, congenital, nongoitrous => probably filter?

View(merge(bachelorLatest %>% mutate(refSampleVaf=RefCount/RefReadDepth) %>% 
       group_by(SampleId,Gene) %>% mutate(netFS=sum(stri_length(Ref)-stri_length(Alt)) %% 3) %>% ungroup() %>%  #check Frameshift per sample
       group_by(Gene,Chromosome,Position,Ref,Alt) %>% mutate(medianRefSampleVaf=median(refSampleVaf),medianNetFS=median(netFS)) %>% ungroup() %>%  #check median vaf
       filter(!(ClinvarSignificance %in% c('Benign/Likely_benign','Benign','Likely_benign')),
              !Effects=='missense_variant'|!(ClinvarDiagnosis %in% c('PI_Z|PI_Z(AUGSBURG)|PI_Z(TUN)|Inborn_genetic_diseases|Alpha-1-antitrypsin_deficiency|FRAXE|not_provided','PI_M(PROCIDA)|Alpha-1-antitrypsin_deficiency')),  #Non cancer related pathogenic
              !(Gene=='GJB2' & HgvsProtein %in% c('p.Leu90Pro')),   #Associated with deafness
              !DbsnpId=='rs111033566',   #very common SNP in PRSS1 (47% of ExAc)
              WorstCodingEffect!='SPLICE'|Type!='INDEL',
              Filter!='ARTEFACT',
              medianRefSampleVaf>0.2,medianRefSampleVaf<0.8,
              !grepl('frameshift',Effects)|medianNetFS!=0),
       highestPurityCohortSummary %>% select(sampleId,primaryTumorLocation),by.x='SampleId',by.y='sampleId',all.x=T) %>%
       group_by(Gene,CodonInfo,ClinvarDiagnosis,ClinvarSignificance,Chromosome,Position,Effects,Ref,Alt,HgvsProtein,DbsnpId,CosmicId,medianRefSampleVaf,Biallelic) %>% count() %>% 
       arrange(-n)  %>% spread(Biallelic,n))

  
