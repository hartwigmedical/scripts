library(ggplot2)

############# Load PURPLE DATA
load("~/hmf/purple.RData")
cohort$age <- ifelse(is.na(cohort$biopsyDate), NA, as.numeric(substr(cohort$biopsyDate,0,4)) - cohort$birthYear)
purpleRaw = cohort
rm(allSamples)
rm(cohort)
rm(backupSamples)

# TRANSFORM PURPLE into usable
purpleVariants = purpleRaw[!is.na(purpleRaw$cancerType), c("snv", "indels", "svTotal", "age", "cancerType")]
purpleVariants$source <- "PURPLE"
colnames(purpleVariants) <- c("snvs", "indels", "svs", "age", "cancerType", "source")

# PURPLE CancerTypes
unique(purpleVariants$cancerType)

#############Load PCAWG Data
pcawgRaw = read.csv("/Users/jon/hmf/pcawg/PCAWG_counts.txt", sep = '\t')
pcawgRaw$cancerType <- NA

# Mapping helpers
#unique(purpleVariants$cancerType)
#unique(pcawgRaw[is.na(pcawgRaw$cancerType), c('organ_system')])
#unique(pcawgRaw[is.na(pcawgRaw$cancerType), c('tumour_type_abbr')])

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
#dim(pcawgRaw[!is.na(pcawgRaw$cancerType),])

# Missing
#pcawgRaw[is.na(pcawgRaw$cancerType),]

# TRANSFORM PCAWG into usable
pcawgSomatics = pcawgRaw[!is.na(pcawgRaw$cancerType), c("all.SNVs", "all.Indels", "SV.events", "age", "cancerType")]
pcawgSomatics$source <- "PCAWG"
colnames(pcawgSomatics) <- c("snvs", "indels", "svs", "age", "cancerType", "source")



#############  Combine into MEGA data frame
purpleIndels = purpleVariants[, c("indels", "cancerType", "age", "source")]
colnames(purpleIndels) <- c("value", "cancerType", "age", "source")
purpleIndels$type <- "Purple Indels"

purpleSnvs = purpleVariants[, c("snvs", "cancerType", "age", "source")]
colnames(purpleSnvs) <- c("value", "cancerType", "age", "source")
purpleSnvs$type <- "Purple SNVs"

purpleSvs = purpleVariants[, c("svs", "cancerType", "age", "source")]
colnames(purpleSvs) <- c("value", "cancerType", "age", "source")
purpleSvs$type <- "Purple SVs"

pcawgIndels = pcawgSomatics[, c("indels", "cancerType", "age", "source")]
colnames(pcawgIndels) <- c("value", "cancerType", "age", "source")
pcawgIndels$type <- "PCAWG Indels"

pcawgSnvs = pcawgSomatics[, c("snvs", "cancerType", "age", "source")]
colnames(pcawgSnvs) <- c("value", "cancerType", "age", "source")
pcawgSnvs$type <- "PCAWG SNVs"

pcawgSvs = pcawgSomatics[, c("svs", "cancerType", "age", "source")]
colnames(pcawgSvs) <- c("value", "cancerType", "age", "source")
pcawgSvs$type <- "PCAWG SVs"

combined = rbind(purpleSvs, purpleIndels, purpleSnvs, pcawgSvs, pcawgIndels, pcawgSnvs)

rm(purpleIndels)
rm(purpleSnvs)
rm(purpleSvs)
rm(pcawgIndels)
rm(pcawgSnvs)
rm(pcawgSvs)

#############  Indel, SNV and SV CDFs
mappedCancerTypes = unique(pcawgSomatics[!is.na(pcawgSomatics$cancerType), c("cancerType")])

ggplot(data=subset(combined, cancerType %in% mappedCancerTypes), aes(x=value, color = type)) + 
  stat_ecdf(geom = "step", pad = FALSE) + 
  facet_grid(cancerType~., scales = "free") + 
  scale_x_log10()

#############  Indels by age
ggplot(data=subset(combined, cancerType %in% mappedCancerTypes & type %in% c('PCAWG Indels', "Purple Indels")),  aes(x=age, y=value, color = type)) + 
  geom_point() + 
  facet_grid(cancerType~.) + 
  scale_y_log10()  + 
  labs(title="Indel count by age")

#ggplot(data=subset(combined, cancerType %in% mappedCancerTypes), 
#       aes(x=value)) + stat_ecdf(geom = "step", pad = FALSE, aes(color = source))  + facet_grid(cancerType ~ type, scales = "free") + scale_x_log10()


#############  EXPERIMENTAL RUBBISH BELOW

#ggplot(data=subset(combined, cancerType %in% c('Breast', 'Melanoma')), 
#       aes(x=value)) + stat_ecdf(geom = "step", pad = FALSE, aes(color = source)) + labs(title="Breast")  + facet_wrap( ~type, scales = "free")


#############  PLOT INDELS
#combinedSomatics = rbind(pcawgSomatics, purpleVariants)

#ggplot(data=subset(combinedSomatics, cancerType %in% c('Breast')), 
#      aes(x=indels)) + stat_ecdf(geom = "step", pad = FALSE, aes(color = source)) + labs(title="Breast Indels")
#
#ggplot(data=subset(combinedSomatics, cancerType %in% c('Melanoma')), 
#       aes(x=indels)) + stat_ecdf(geom = "step", pad = FALSE, aes(color = source)) + labs(title="Melanoma Indels")



#ggplot(data=subset(combinedSomatics, cancerType %in% c('Breast')), 
#       aes(x=snvs)) + stat_ecdf(geom = "step", pad = FALSE, aes(color = source)) + labs(title="Breast SNVs")

#ggplot(data=subset(combinedSomatics, cancerType %in% c('Melanoma')), 
#       aes(x=snvs)) + stat_ecdf(geom = "step", pad = FALSE, aes(color = source)) + labs(title="Melanoma SNVs")

#aes(chCopyNumber)
#ggplot(aes(indels), data=allSomatics)+ stat_ecdf(geom = "step", pad = FALSE) + facet_wrap( ~cancerType )


#purpleBreast = allSomatics[allSomatics$cancerType == 'Breast',]
#purpleLung = allSomatics[allSomatics$cancerType == 'Lung',]

#ggplot(aes(indels), data=purpleBreast) + stat_ecdf(geom = "step", pad = FALSE)
#ggplot(aes(indels), data=purpleLung) + stat_ecdf(geom = "step", pad = FALSE)



#gplot(data=subset(allSomatics, cancerType %in% c('Breast', 'Lung')), 
#       aes(x=indels, color=cancerType))+ stat_ecdf(geom = "step", pad = FALSE) 

#ggplot(data=subset(allSomatics, cancerType %in% c('Breast', 'Lung')), 
#       aes(x=indels))+ stat_ecdf(geom = "step", pad = FALSE, aes(color = cancerType)) 


