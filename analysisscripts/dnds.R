library(RMySQL)
library(dndscv)
library(IRanges)
detach("package:purple", unload=TRUE); library(purple);

dbProd = dbConnect(MySQL(), dbname='hmfpatients', groups="RAnalysis")
cat("Querying purple")
rawCohort = purple::query_purity(dbProd)
cohort = rawCohort

# PatientIds
cat("Mapping samples to patients")
patientIdLookups = query_patient_id_lookup(dbProd)
patientIds = purple::apply_to_cohort(cohort, function(x) {purple::sample_to_patient_id(x$sampleId, patientIdLookups)})
cohort$patientId <- patientIds$V1

#Clinical Data
cat("Querying clinical data")
clinicalData = purple::query_clinical_data(dbProd)
cohort = left_join(cohort, clinicalData[, c("sampleId", "cancerType")])

# Cohort
cohort = purple::highest_purity_patients(cohort)
save(cohort, file="/Users/jon/hmf/RData/dnds/dndsCohort.RData")

# Somatics 
cat("Querying somatics")
rawSomatics = purple::query_somatic_variants(dbProd, cohort, filterEmptyGenes = TRUE)
save(rawSomatics, file="/Users/jon/hmf/RData/dnds/dndsSomatics.RData")
somatics = rawSomatics[rawSomatics$type == "INDEL" | rawSomatics$type == "SNP", c("sampleId", "chromosome", "position", "ref", "alt")]
colnames(somatics) <- c("sampleId", "chr", "pos", "ref", "alt")

# Takes AGES to annotate with cancer type... don't bother!
#somatics$cancerType <- "UNKNOWN"
#for (sampleId in unique(somatics$sampleId)) {
#  matchedCancerType = cohort[match(sampleId, cohort$sampleId), c("cancerType")]
#  cat("SampleId:", sampleId, ", cancerType:", matchedCancerType, "\n")
#  somatics[somatics$sampleId == sampleId, ]$cancerType <- matchedCancerType
#}

# Clean up DB Connection
dbDisconnect(dbProd)
rm(dbProd)
rm(patientIdLookups)
rm(patientIds)
rm(clinicalData)

cancerTypes = unique(cohort$cancerType)
cancerTypes = cancerTypes[!is.na(cancerTypes)]

dndsResults = list()
for (cancerType in cancerTypes[!is.na(cancerTypes)]) {
  cat("Processing", cancerType)
  cancerTypeSampleIds = cohort[!is.na(cohort$cancerType) & cohort$cancerType == cancerType, c("sampleId")]
  input = somatics[somatics$sampleId %in% cancerTypeSampleIds, c("sampleId", "chr", "pos", "ref", "alt")]
  output = dndscv(input)
  dndsResults[[cancerType]] <- output$sel_cv
  save(dndsResults, file="/Users/jon/hmf/RData/dnds/dndsSelCV.RData")
}

output = dndscv(somatics)

output$sel_cv
dndsResults[["Pan Cancer"]] <- output$sel_cv
names(dndsResults)

panCancer <-  output$sel_cv

### WORKING
pcawgRaw = read.csv("/Users/jon/hmf/pcawg/PCAWG_counts.txt", sep = '\t')
pcawgPanCancer = read.csv("/Users/jon/hmf/pcawg/dNdScv_output_PANCANCER.txt", sep = "\t")


pcawgPanCancer

pcawgPanCancer = read.csv("/Users/jon/hmf/pcawg/dNdScv_output_PANCANCER.txt", sep = "\t")
hmfPanCancer = dndsResults[["Pan Cancer"]]
panCancerCombined = merge(pcawgPanCancer, hmfPanCancer, all.X = TRUE, all.Y = TRUE, by = "gene_name", suffixes = c(".pcawg", ".hmf"))
exclusiveToHMF = panCancerCombined[panCancerCombined$qglobal > 0.05 & panCancerCombined$qglobal_cv < 0.05, ]
exclusiveToPCAWG = panCancerCombined[panCancerCombined$qglobal < 0.05 & panCancerCombined$qglobal_cv > 0.05, ]
inBoth = panCancerCombined[panCancerCombined$qglobal < 0.05 & panCancerCombined$qglobal_cv < 0.05, ]

pcawgSignificant = panCancerCombined[panCancerCombined$qglobal < 0.05, ] 
hmfSignificant = panCancerCombined[panCancerCombined$qglobal_cv < 0.05, ] 

dt = data.table(cohort)
dt[, .N, by = cancerType][order(N)]


names(hmfPanCancer)
write.table(hmfPanCancer[1:100, ], file = "/Users/jon/hmf/repos/hmftools/hmf-common/src/test/resources/dndscv/dndscv_long.tsv", sep="\t", quote=F, row.names = F)


##### DRIVERS
hotspots = read.csv(file="/Users/jon/hmf/hotspot/Hotspot.tsv", header=FALSE, sep="\t")
colnames(hotspots) <- c("chromosome", "position", "ref", "alt")

isHotspot<-function(variant, hotspots) {
  result = hotspots[hotspots$chromosome == variant[1] & hotspots$position == variant[2] & hotspots$ref == variant[3] & hotspots$alt == variant[4],  ]
  return (nrow(result) != 0)
}



ACVR2A
head(rawSomatics)

tp53Somatics = purple::query_somatic_variants(dbProd, cohort, filterEmptyGenes = TRUE)
tp53Somatics[tp53Somatics$chromosome == 17 & tp53Somatics$position == 7578406, ]
hotspots[hotspots$chromosome == 17 & hotspots$position == 7578406, ]


brafSomatics = purple::query_somatic_variants(dbProd, cohort, gene = "ACVR2A")
acvr2aSomatics = purple::query_somatic_variants(dbProd, cohort, gene = "ACVR2A")
jon = acvr2aSomatics[acvr2aSomatics$type == "INDEL" & acvr2aSomatics$effect == "frameshift variant", ]

save(brafSomatics, tp53Somatics, file="/Users/jon/hmf/RData/dnds/dndsBRAFSomatics.RData")

variant = brafSomatics[1, ]
result = hotspots[hotspots$chromosome == variant[1] & hotspots$position == variant[2] & hotspots$ref == variant[3] & hotspots$alt == variant,  ]
nrow(result)

brafSomatics$hotspot = apply(brafSomatics[, c("chromosome", "position", "ref", "alt")], 1, function (x) {isHotspot(x, hotspots)})
tp53Somatics$hotspot = apply(tp53Somatics[, c("chromosome", "position", "ref", "alt")], 1, function (x) {isHotspot(x, hotspots)})

dt = data.table(brafSomatics)
brafStats = dt[,
   list(
    biallelicHotspot = sum(.SD$hotspot == TRUE & .SD$loh == 1),  
    diploidHotspot = sum(.SD$hotspot == TRUE & .SD$loh == 0),  
    biallelic = sum(.SD$hotspot == FALSE & .SD$loh == 1),  
    remainder = sum(.SD$hotspot == FALSE & .SD$loh == 0),  
    total = .N), 
  by=.(sampleId, type)]

dt[sampleId == "CPCT02030213T", 
     print(.SD), 
   by=.(sampleId, type)]

dt[sampleId == 'CPCT02030213T',  ]

dt[hotspot == TRUE, .N, by=sampleId ][order(-N)]

dcast(dt, sampleId + hotspot ~ type)


any(jon)

apply(brafSomatics[1:3, ], 1, function (x) {str(x[1])})
apply(brafSomatics[1:3, ], 1, function (x) {str(x)})

ifelse(isHotspot(brafSomatics, hotspots), TRUE, FALSE)

?apply()
isHotspot(brafSomatics[1, ], hotspots)

sum
unique(brafSomatics$effect)
?grep
grep("stop", "splice region variant; intron variant")

standardEffect<-function(variant) {
  effect = variant$effect
  
  if (variant$type == "INDEL") {
    if (grep("stop gained", effect) || grep("frameshift", effect)) {
      return ("IndelNon")
    }
    return ()
  }
  
}


load(file="/Users/jon/hmf/RData/dnds/dndsSomatics.RData")
jon2 = rawSomatics[rawSomatics$chromosome == 2 & rawSomatics$position == 148683685,]
nrow(jon2)
