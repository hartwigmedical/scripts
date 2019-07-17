library(RMySQL)
library(dplyr)
library(tidyr)
library(GenomicRanges)

########################### Cosmic Genes
cosmicCurated = read.csv("~/hmf/resources/CosmicCurated.csv", stringsAsFactors = F)
cosmicCurated$cosmicCurated <- TRUE
cosmicCurated = cosmicCurated[, c("Genes","cosmicCurated")]
colnames(cosmicCurated) <- c("gene", "cosmicCurated")

cosmicCensus = read.csv("~/hmf/resources/CosmicCensus.csv", stringsAsFactors = F)
colnames(cosmicCensus) <- c("gene", "cosmicOncogene", "cosmicTsg")
cosmicCensus = cosmicCensus %>% filter(cosmicOncogene | cosmicTsg)

cosmicGenes = merge(cosmicCurated, cosmicCensus, by = "gene", all = T)
rm(cosmicCurated, cosmicCensus)

########################### Actionable Genes
actionableGenes = read.table("~/hmf/analysis/actionable/actionablePanel.tsv", header = T, stringsAsFactors = F) 
colnames(actionableGenes) <- c("gene", "actionableAmplification", "actionableDeletion", "actionableVariant", "actionableDrup", "actionableResponse", "actionableResponseSource", "actionableResistance", "actionableResistanceSource")
actionableGenes$actionableAmplification <- ifelse(actionableGenes$actionableAmplification  == "true", T, NA)
actionableGenes$actionableDeletion <- ifelse(actionableGenes$actionableDeletion  == "true", T, NA)
actionableGenes$actionableVariant <- ifelse(actionableGenes$actionableVariant  == "true", T, NA)
actionableGenes$actionableDrup <- ifelse(actionableGenes$actionableDrup  == "true", T, NA)

########################### KNOWN AMPS AND DELS
cgi = read.table(file = "~/hmf/resources/cgi_biomarkers_per_variant_181127.tsv", stringsAsFactors = F, header = T, sep = "\t") %>%
    select(Biomarker) %>%
    separate(Biomarker, c("gene","alteration"), sep = " ") %>%
    filter(alteration %in% c("amplification", "deletion")) %>%
    distinct(gene, alteration)

onco = read.table(file = "~/hmf/resources/onco_allAnnotatedVariants_181127.csv", stringsAsFactors = F, header = T, sep = "\t", quote = "") %>%
    select(gene = Gene, alteration = Alteration) %>%
    mutate(alteration = tolower(alteration)) %>%
    filter(alteration %in% c("amplification", "deletion")) %>%
    distinct(gene, alteration)

civic = read.table(file = "~/hmf/resources/civic_01-Nov-2018-VariantSummaries_181127.tsv", stringsAsFactors = T, header = T, sep = "\t", quote = "", fill = T) %>%
  select(gene, alteration = variant) %>%
    mutate(alteration = tolower(alteration)) %>%
    filter(alteration %in% c("amplification", "deletion")) %>%
    distinct(gene, alteration)


knownAmpsDels = bind_rows(cgi, onco) %>%
    bind_rows(civic) %>%
    select(gene_name = gene, alteration) %>%
    distinct(gene_name, alteration) %>%
    mutate(value = T) %>% spread(alteration, value)
rm(cgi, onco, civic)
colnames(knownAmpsDels) <- c("gene", "knownAmplification", "knownDeletion")


########################### Gene Panel
load(file="~/hmf/dnds/PcawgRefCDSCv.RData")
load(file="~/hmf/analysis/cohort/processed/HmfRefCDSCv.RData")
sig = 0.01
hmfSignificant =  HmfRefCDSCv %>% filter(qglobal_cv < sig) %>% distinct(gene_name) %>% select(gene = gene_name)
hmfSignificant$hmfDnds <- TRUE
martincorenaSignificant =  PcawgRefCDSCv %>% filter(qglobal < sig) %>% distinct(gene_name) %>% select(gene = gene_name)
martincorenaSignificant$martincorenaDnds <- TRUE

genePanel = merge(martincorenaSignificant, hmfSignificant, by = "gene", all = T)
genePanel = merge(genePanel, cosmicGenes, by = "gene", all = T)
genePanel = merge(genePanel, actionableGenes, by = "gene", all = T)
genePanel = merge(genePanel, knownAmpsDels, by = "gene", all=T)

genePanelInitial = genePanel %>% filter(!gene %in% c("POM121L12","TRIM49B","LPCAT2"))
save(genePanelInitial, file="~/hmf/analysis/cohort/processed/genePanelInitial.RData")


########################### PART 2 - GOTO ampsDelsTarget


########################### PART 3 - GOTO dndsClassification


####### WRITE TO REPOSITORY AND DB
load(file = "~/hmf/analysis/cohort/processed/genePanel.RData")

reportableDels = genePanel %>% 
  filter(reportableDel) %>% 
  mutate(
    loci = coalesce(centromere, telomere), 
    loci = gsub("_", " ", loci),
    loci = gsub("[0-9]{1,2}", "", loci),
    loci = gsub("[X|Y]", "", loci)) %>% 
  select(gene, loci) %>% arrange(is.na(loci), gene)

write.table(genePanel %>% filter(reportableAmp) %>% select(gene), file = "/Users/jon/hmf/repos/hmftools/hmf-common/src/main/resources/cna/AmplificationTargets.tsv", quote = F, row.names = F, col.names = F)
write.table(reportableDels, file = "/Users/jon/hmf/repos/hmftools/hmf-common/src/main/resources/cna/DeletionTargets.tsv", quote = F, row.names = F, col.names = F, sep = "\t")

genePanelDB = genePanel %>%
  mutate(
    actionableResponse = ifelse(is.na(actionableResponse), "", actionableResponse),
    actionableResponseSource = ifelse(is.na(actionableResponseSource), "", actionableResponseSource),
    actionableResistance = ifelse(is.na(actionableResistance), "", actionableResistance),
    actionableResistanceSource = ifelse(is.na(actionableResistanceSource), "", actionableResistanceSource),
    reportablePointMutation = ifelse(is.na(reportablePointMutation), "", reportablePointMutation),
    armEndLocus = coalesce(telomere, centromere, "")
    ) %>%
  select(-telomere, -centromere) %>%
  select(gene,martincorenaDnds,hmfDnds,cosmicCurated,cosmicOncogene,cosmicTsg,actionableAmplification,actionableDeletion,actionableVariant,actionableDrup,actionableResponse,actionableResponseSource,
         actionableResistance,actionableResistanceSource,knownAmplification,knownDeletion,hmfAmplification,hmfDeletion,armEndLocus,reportablePointMutation,reportableAmp, reportableDel) %>%
  filter(is.na(reportablePointMutation) | reportableDel | reportableAmp | cosmicOncogene | cosmicTsg | actionableDrup | actionableVariant | knownAmplification | knownDeletion )

genePanelDB[is.na(genePanelDB)] <- 0

dbLocal = dbConnect(MySQL(), dbname='hmfpatients_pilot', user = "build", password = "build", host = "localhost")
dbLocal = dbConnect(MySQL(), dbname='hmfpatients_pilot', groups="RAnalysisWrite", host = "127.0.0.1")
dbLocal = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysisWrite", host = "127.0.0.1")
dbWriteTable(dbLocal, value = genePanelDB, name = "genePanel", overwrite = T, row.names = F )
dbDisconnect(dbLocal)
rm(dbLocal)


########################### PART 4 Check with hotspots

load("~/hmf/analysis/cohort/reference/canonicalTranscripts.RData")
load("~/hmf/analysis/cohort/processed/genePanel.RData")
cgi = read.csv("~/hmf/resources/cgi_variant_list", sep ="\t", stringsAsFactors = F)
cgi = cgi %>% 
  select(gene = GENE, chromosome = CHROMOSOME, position = POSITION, ref = REF, alt = ALT) %>%
  group_by(gene, chromosome, position, ref, alt) %>% distinct(gene, chromosome, position, ref, alt)

onco = read.csv("~/hmf/resources/oncoKb_variant_list", sep ="\t", stringsAsFactors = F)
onco = onco %>% 
  filter(INFO %in% c("Likely Oncogenic", "Oncogenic")) %>%
  select(gene = GENE, chromosome = CHROMOSOME, position = POSITION, ref = REF, alt = ALT) %>%
  group_by(gene, chromosome, position, ref, alt) %>% distinct(gene, chromosome, position, ref, alt)

civic = read.csv("~/hmf/resources/civic_variant_list", sep ="\t", stringsAsFactors = F)
civic = civic %>% 
  filter(INFO == "true") %>%
  select(gene = GENE, chromosome = CHROMOSOME, position = POSITION, ref = REF, alt = ALT) %>%
  group_by(gene, chromosome, position, ref, alt) %>% distinct(gene, chromosome, position, ref, alt) 


hotspots = full_join(cgi, onco, by = c("gene", "chromosome", "position", "ref", "alt")) %>% full_join(civic, by = c("gene", "chromosome", "position", "ref", "alt")) %>% 
  left_join(genePanel %>% select(gene) %>% mutate(inGenePanel = T), by = "gene")

hotspotsMissing = hotspots %>% filter(is.na(inGenePanel))
length(unique(hotspotsMissing$gene))

hotspotgenes = full_join(cgi, onco, by = "gene") %>% full_join(civic, by = "gene") %>% left_join(genePanel %>% select(gene) %>% mutate(inGenePanel = T), by = "gene")
missingHotspotgenes = hotspotgenes %>% filter(is.na(inGenePanel))
canonicalTranscripts %>% filter(gene %in% hotspotsMissing$gene)

head(hpcExonicSomatics)

load(file = "~/hmf/analysis/cohort/reference/hpcExonicSomatics.RData")
hotspotSomatics = hpcExonicSomatics %>% inner_join(hotspotsMissing, by = c("chromosome", "position", "ref", "alt")) %>% group_by(gene.x) %>% count()
  
hotspotSomatics = hpcExonicSomatics %>% inner_join(hotspotsMissing, by = c("chromosome", "position", "ref", "alt")) %>% group_by(gene.x)
  
missingHotspotCv = hotspotSomatics %>% left_join(HmfRefCDSCv, by = c("gene.x" = "gene_name")) %>% filter(cancerType == 'All')
missingHotspotCv$n






########################### PART Check actionable variants

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO

actionableVariant = read.csv(file = "~/hmf/analysis/actionable/actionableVariants.tsv", sep = "\t") %>% select(2,3,4,5)
names(actionableVariant) <- c("CHROM", "POS", "REF", "ALT")

actionableVariant = actionableVariant %>% mutate(ID = ".", QUAL = "100", FILTER = ".") %>%
  select(CHROM,	POS,	ID,	REF,	ALT,	QUAL,	FILTER) %>% mutate(CHROM = factor(CHROM, levels = c(1:22, 'X', 'Y'))) %>% arrange(CHROM, POS) %>% distinct()

write.table(actionableVariant, file = "~/hmf/analysis/actionable/actionableVariants.vcf", sep = "\t", row.names = F, quote = F)

load("/Users/jon/hmf/analysis/cohort/processed/genePanel.RData")
result = read.table(file = "/Users/jon/hmf/analysis/cohort/processed/actionableVariants.ann.vcf.info.summary", sep  = "|")
names(result) <- c("impact", "gene")

actionalVariantCategories = result %>%
  group_by(gene, impact) %>% count() %>% 
  spread(impact, n) %>% 
  left_join(genePanel %>% select(gene, reportablePointMutation))
save(actionalVariantCategories, file = "~/hmf/analysis/actionable/actionalVariantCategories.RData")
