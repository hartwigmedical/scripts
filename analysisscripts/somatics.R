library(ggplot2)


############# Load PURPLE DATA
load("~/hmf/somatics.RData")
purpleRaw = allSomatics
rm(allSomatics)

# TRANSFORM PURPLE into usable
purpleSomatics = purpleRaw[!is.na(purpleRaw$cancerType), c("snv", "indels", "cancerType")]
purpleSomatics$source <- "PURPLE"
colnames(purpleSomatics) <- c("snvs", "indels", "cancerType", "source")

# PURPLE CancerTypes
unique(purpleSomatics$cancerType)

#############Load PCAWG Data
pcawgRaw = read.csv("/Users/jon/hmf/pcawg/PCAWG_counts.txt", sep = '\t')
pcawgRaw$cancerType <- NA

# Mapping helpers
unique(purpleSomatics$cancerType)
unique(pcawgRaw[is.na(pcawgRaw$cancerType), c('organ_system')])
unique(pcawgRaw[is.na(pcawgRaw$cancerType), c('tumour_type_abbr')])

### Map PCAWG cancerTypes to PURPLE cancerType
pcawgRaw$cancerType <- ifelse(pcawgRaw$tumour_type_abbr == "Melanoma", "Melanoma", pcawgRaw$cancerType)
pcawgRaw$cancerType <- ifelse(pcawgRaw$tumour_type_abbr == "Sarcoma_soft_tissue", "Sarcoma", pcawgRaw$cancerType)
pcawgRaw$cancerType <- ifelse(pcawgRaw$organ_system == "BREAST", "Breast", pcawgRaw$cancerType)
pcawgRaw$cancerType <- ifelse(pcawgRaw$organ_system == "LIVER", "Liver", pcawgRaw$cancerType)
pcawgRaw$cancerType <- ifelse(pcawgRaw$organ_system == "KIDNEY", "Kidney", pcawgRaw$cancerType)
pcawgRaw$cancerType <- ifelse(pcawgRaw$organ_system == "OVARY", "Ovary", pcawgRaw$cancerType)
pcawgRaw$cancerType <- ifelse(pcawgRaw$organ_system == "PROSTATE GLAND", "Prostate", pcawgRaw$cancerType)
pcawgRaw$cancerType <- ifelse(pcawgRaw$organ_system == "STOMACH", "Stomach", pcawgRaw$cancerType)
pcawgRaw$cancerType <- ifelse(pcawgRaw$organ_system == "PANCREAS", "Pancreas", pcawgRaw$cancerType)
pcawgRaw$cancerType <- ifelse(pcawgRaw$organ_system == "UTERUS, NOS", "Uterus", pcawgRaw$cancerType)
pcawgRaw$cancerType <- ifelse(pcawgRaw$organ_system == "BRAIN, & CRANIAL NERVES, & SPINAL CORD, (EXCL. VENTRICLE, CEREBELLUM)", "Brain", pcawgRaw$cancerType)
pcawgRaw$cancerType <- ifelse(pcawgRaw$organ_system == "LUNG & BRONCHUS", "Lung", pcawgRaw$cancerType)
pcawgRaw$cancerType <- ifelse(pcawgRaw$organ_system == "GALLBLADDER & EXTRAHEPATIC BILE DUCTS", "Gall bladder", pcawgRaw$cancerType)

# Mapped successfully
dim(pcawgRaw[!is.na(pcawgRaw$cancerType),])

# Missing
missing = pcawgRaw[is.na(pcawgRaw$cancerType),]

# TRANSFORM PCAWG into usable
pcawgSomatics = pcawgRaw[!is.na(pcawgRaw$cancerType), c("all.SNVs", "all.Indels", "cancerType")]
pcawgSomatics$source <- "PCAWG"
colnames(pcawgSomatics) <- c("snvs", "indels", "cancerType", "source")

#############  Combine into MEGA data frame
purpleIndels = purpleSomatics[, c("indels", "cancerType", "source")]
colnames(purpleIndels) <- c("value", "cancerType", "source")
purpleIndels$type <- "Indels"

purpleSnvs = purpleSomatics[, c("snvs", "cancerType", "source")]
colnames(purpleSnvs) <- c("value", "cancerType", "source")
purpleSnvs$type <- "SNVs"

pcawgIndels = pcawgSomatics[, c("indels", "cancerType", "source")]
colnames(pcawgIndels) <- c("value", "cancerType", "source")
pcawgIndels$type <- "Indels"

pcawgSnvs = pcawgSomatics[, c("snvs", "cancerType", "source")]
colnames(pcawgSnvs) <- c("value", "cancerType", "source")
pcawgSnvs$type <- "SNVs"

combined = rbind(purpleIndels, purpleSnvs, pcawgIndels, pcawgSnvs)
rm(purpleIndels)
rm(purpleSnvs)
rm(pcawgIndels)
rm(pcawgSnvs)

#############  PLOT

mappedCancerTypes = unique(pcawgSomatics[!is.na(pcawgSomatics$cancerType), c("cancerType")])
ggplot(data=subset(combined, cancerType %in% mappedCancerTypes), 
       aes(x=value)) + stat_ecdf(geom = "step", pad = FALSE, aes(color = source))  + facet_grid(cancerType ~ type, scales = "free")


#############  EXPERIMENTAL RUBBISH BELOW

ggplot(data=subset(combined, cancerType %in% c('Breast', 'Melanoma')), 
       aes(x=value)) + stat_ecdf(geom = "step", pad = FALSE, aes(color = source)) + labs(title="Breast")  + facet_wrap( ~type, scales = "free")


#############  PLOT INDELS
combinedSomatics = rbind(pcawgSomatics, purpleSomatics)

ggplot(data=subset(combinedSomatics, cancerType %in% c('Breast')), 
       aes(x=indels)) + stat_ecdf(geom = "step", pad = FALSE, aes(color = source)) + labs(title="Breast Indels")

ggplot(data=subset(combinedSomatics, cancerType %in% c('Melanoma')), 
       aes(x=indels)) + stat_ecdf(geom = "step", pad = FALSE, aes(color = source)) + labs(title="Melanoma Indels")



ggplot(data=subset(combinedSomatics, cancerType %in% c('Breast')), 
       aes(x=snvs)) + stat_ecdf(geom = "step", pad = FALSE, aes(color = source)) + labs(title="Breast SNVs")

ggplot(data=subset(combinedSomatics, cancerType %in% c('Melanoma')), 
       aes(x=snvs)) + stat_ecdf(geom = "step", pad = FALSE, aes(color = source)) + labs(title="Melanoma SNVs")

aes(chCopyNumber)
ggplot(aes(indels), data=allSomatics)+ stat_ecdf(geom = "step", pad = FALSE) + facet_wrap( ~cancerType )


purpleBreast = allSomatics[allSomatics$cancerType == 'Breast',]
purpleLung = allSomatics[allSomatics$cancerType == 'Lung',]

ggplot(aes(indels), data=purpleBreast) + stat_ecdf(geom = "step", pad = FALSE)
ggplot(aes(indels), data=purpleLung) + stat_ecdf(geom = "step", pad = FALSE)



ggplot(data=subset(allSomatics, cancerType %in% c('Breast', 'Lung')), 
       aes(x=indels, color=cancerType))+ stat_ecdf(geom = "step", pad = FALSE) 

ggplot(data=subset(allSomatics, cancerType %in% c('Breast', 'Lung')), 
       aes(x=indels))+ stat_ecdf(geom = "step", pad = FALSE, aes(color = cancerType)) 


