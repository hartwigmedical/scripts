library(ggplot2)

load(file = "~/hmf/analysis/cohort/processed/genePanel.RData")
genePanel = genePanel %>% filter(!is.na(reportablePointMutation) | !is.na(reportableAmp) | !is.na(reportableDel)) %>% pull(gene)

purity = data.frame(
SampleId = c("A-CPCT02010944T", "U-CPCT02010963T", "P-CPCT02020378T", "A-CPCT02020629T","A-CPCT02020787T", "U-CPCT02020834T", "P-CPCT02080180T", "P-CPCT02100137T", "U-CPCT02360014T", "P-DRUP01010040T"),
Purity = c(0.58,0.93,	0.63,	0.2,	0.95,	1,	0.88,	0.94,	0.9,	0.88), stringsAsFactors = F)

ensemblTransExonData = read.csv('~/hmf/analysis/RNA/ensembl_trans_exon_data.csv')
ensemblGeneData = read.csv('~/hmf/analysis/RNA/ensembl_gene_data.csv')
allTranscripts = ensemblTransExonData %>% distinct(GeneId,Trans)

salmonTest = rbind(
read.csv('/Users/jon/hmf/analysis/salmon/CPCT02010944T.salmon_quant/quant.sf',sep='\t') %>% mutate(SampleId='A-CPCT02010944T'),
read.csv('/Users/jon/hmf/analysis/salmon/CPCT02020629T.salmon_quant/quant.sf',sep='\t') %>% mutate(SampleId='A-CPCT02020629T'),
read.csv('/Users/jon/hmf/analysis/salmon/CPCT02020787T.salmon_quant/quant.sf',sep='\t') %>% mutate(SampleId='A-CPCT02020787T'),
read.csv('/Users/jon/hmf/analysis/salmon/CPCT02360014T.salmon_quant/quant.sf',sep='\t') %>% mutate(SampleId='U-CPCT02360014T'),
read.csv('/Users/jon/hmf/analysis/salmon/CPCT02010963T.salmon_quant/quant.sf',sep='\t') %>% mutate(SampleId='U-CPCT02010963T'),
read.csv('/Users/jon/hmf/analysis/salmon/CPCT02020834T.salmon_quant/quant.sf',sep='\t') %>% mutate(SampleId='U-CPCT02020834T'),
read.csv('/Users/jon/hmf/analysis/salmon/CPCT02020378T.salmon_quant/quant.sf',sep='\t') %>% mutate(SampleId='P-CPCT02020378T'),
read.csv('/Users/jon/hmf/analysis/salmon/CPCT02080180T.salmon_quant/quant.sf',sep='\t') %>% mutate(SampleId='P-CPCT02080180T'),
read.csv('/Users/jon/hmf/analysis/salmon/CPCT02100137T.salmon_quant/quant.sf',sep='\t') %>% mutate(SampleId='P-CPCT02100137T'),
read.csv('/Users/jon/hmf/analysis/salmon/DRUP01010040T.salmon_quant/quant.sf',sep='\t') %>% mutate(SampleId='P-DRUP01010040T'))

save(ensemblTransExonData, ensemblGeneData, allTranscripts, salmonTest, genePanel, purity, file = "/Users/jon/hmf/analysis/salmon/salmonTest.RData")
rm(list=setdiff(ls(), c("ensemblTransExonData","allTranscripts","ensemblGeneData","salmonTest", "genePanel")))

unique(filteredSalmonTest$SampleId)



load(file = "/Users/jon/hmf/analysis/salmon/salmonTest.RData")

filteredSalmonTest = salmonTest %>%
    group_by(SampleId) %>%
    left_join(purity, by = "SampleId") %>%
    mutate(CancerType = substr(SampleId, 1, 1)) %>%
    mutate(normalisationFactor = sum(ifelse(TPM > 5000,0, NumReads) / EffectiveLength)) %>%
    ungroup() %>%
    mutate(TPM = NumReads / EffectiveLength / normalisationFactor * 1000000) %>%
    select(-normalisationFactor) %>%
    left_join(allTranscripts %>% select(Name = Trans, GeneId), by = c("Name")) %>%
    left_join(ensemblGeneData %>% select(GeneId, GeneName), by = "GeneId") %>%
    mutate( GenePanel = GeneName %in% genePanel) %>%
    filter(GenePanel)

geneData = filteredSalmonTest %>% group_by(SampleId, GeneName, Purity, CancerType) %>% summarise(TPM = sum(TPM)) %>% ungroup()
ggplot(geneData %>% filter(GeneName == "B2M"), aes(Purity, TPM)) + geom_point(aes(color = CancerType))

urinaryData = geneData %>% filter(grepl("^A", SampleId) | grepl("^U", SampleId), !is.na(GeneName))
prostateData = geneData %>% filter(grepl("^P", SampleId), !is.na(GeneName))


summariseAndCombine <- function(trainingData, testData) {
    trainingDataSummary = trainingData %>%
        filter(!is.na(GeneName)) %>%
        group_by(GeneName) %>%
        summarise(meanTrainingTPM = mean(TPM), sdTrainingTPM = sd(TPM))

    testDataSummary = testData %>%
        filter(!is.na(GeneName)) %>%
        group_by(GeneName) %>%
        summarise(meanTestTPM = mean(TPM), sdTestTPM = sd(TPM))

    result = full_join(trainingDataSummary, testDataSummary, by = "GeneName") %>%
        mutate(diff = abs(meanTestTPM - meanTrainingTPM), foldChange = meanTestTPM / meanTrainingTPM, z = diff / sdTrainingTPM) %>%
        filter(meanTrainingTPM > 0) %>%
        arrange(-z)
    return (result)
}

sampleOutliers <- function(sample, data) {
    trainingData = data %>% filter(SampleId != sample)
    testData = data  %>% filter(SampleId == sample)

    result = summariseAndCombine(trainingData, testData) %>%
    filter(meanTrainingTPM > 0, z > 2, diff > 10)
    return (result)
}

urinarySamples = unique(urinaryData$SampleId)
urinaryOutliers1 = sampleOutliers(urinarySamples[1], urinaryData)
urinaryOutliers2 = sampleOutliers(urinarySamples[2], urinaryData)
urinaryOutliers3 = sampleOutliers(urinarySamples[3], urinaryData)
urinaryOutliers4 = sampleOutliers(urinarySamples[4], urinaryData)
urinaryOutliers5 = sampleOutliers(urinarySamples[5], urinaryData)
urinaryOutliers6 = sampleOutliers(urinarySamples[6], urinaryData)

prostateSamples = unique(prostateData$SampleId)
prostateOutliers1 = sampleOutliers(prostateSamples[1], prostateData)
prostateOutliers2 = sampleOutliers(prostateSamples[2], prostateData)
prostateOutliers3 = sampleOutliers(prostateSamples[3], prostateData)
prostateOutliers4 = sampleOutliers(prostateSamples[4], prostateData)

cohortOutliers = summariseAndCombine(urinaryData, prostateData) %>%
filter(z > 3, diff > 10)


ggplot(cohortOutliers, aes(meanTrainingTPM, meanTestTPM)) +
    geom_point() +
    geom_text(aes(label = GeneName), size = 3, nudge_x = .01, color = "red") +
#xlim(0, 50) + ylim(0, 50) +
    xlab("Prostate TPM") + ylab("Urinary TPM")  + scale_x_log10() + scale_y_log10()



ggplot(urinaryOutliers1 %>% filter(diff > 30), aes(meanTrainingTPM, meanTestTPM)) +
    geom_point() +
    geom_text(aes(label = GeneName), size = 3, nudge_x = 0.01, color = "red") +
    xlab("Cohort TPM") + ylab("Sample TPM") + scale_x_log10() + scale_y_log10()









salmonTPM = filteredSalmonTest %>%  mutate(TPM=round(TPM,1)) %>% spread(SampleId,TPM)
salmonTPMGene = filteredSalmonTest %>% filter(!is.na(GeneName)) %>% mutate(TPM=round(TPM,1)) %>% group_by(GeneName,SampleId) %>% summarise(TPM=sum(TPM)) %>% spread(SampleId,TPM)

View(salmonTPMGene %>% filter(GeneName=='AHR'))

ggplot(data=salmonTPMGene %>% filter(`A-CPCT02010944T`>100|`A-CPCT02020629T`>100),aes(`A-CPCT02010944T`,`A-CPCT02020629T`))+geom_point() + geom_text(aes(label = GeneName), size = 3, nudge_x = 100, color = "red") #+scale_x_log10() +scale_y_log10()
ggplot(data=salmonTPMGene %>% filter(`A-CPCT02010944T`>100|`P-CPCT02020378T`>100),aes(`A-CPCT02010944T`,`P-CPCT02020378T`))+geom_point()+ geom_text(aes(label = GeneName), size = 3, nudge_x = 100, color = "red") #+scale_x_log10() +scale_y_log10()
ggplot(data=salmonTPMGene %>% filter(`A-CPCT02010944T`>100|`U-CPCT02020834T`>100),aes(`A-CPCT02010944T`,`U-CPCT02020834T`))+geom_point() + geom_text(aes(label = GeneName), size = 3, nudge_x = 100, color = "red") #+scale_x_log10() +scale_y_log10()


# High Purity urinary tract vs NO_TUMOR
ggplot(data=salmonTPMGene %>% filter(`U-CPCT02010963T`+`U-CPCT02020834T`>100),aes(`U-CPCT02010963T`,`U-CPCT02020834T`))+geom_point() + geom_text(aes(label = GeneName), size = 3, nudge_x = 100, color = "red") #+scale_x_log10() +scale_y_log10()

# High Purity prostate vs NO_TUMOR
ggplot(data=salmonTPMGene %>% filter(`P-CPCT02100137T`+`U-CPCT02020834T`>100),aes(`P-CPCT02100137T`,`U-CPCT02020834T`))+geom_point() + geom_text(aes(label = GeneName), size = 3, nudge_x = 0.1, color = "red") +scale_x_log10() +scale_y_log10()
