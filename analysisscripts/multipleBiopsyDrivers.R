library(dplyr)
detach("package:purple", unload=TRUE);
library(purple)

load(file = "~/hmf/RData/reference/multipleBiopsyCohort.RData")
load(file = "~/hmf/RData/reference/mbExonicSomatics.RData")
load("~/hmf/RData/processed/excessRates.RData")
load(file = "~/hmf/RData/processed/dndsUnfilteredMultipleBiopsyMutations.RData")
dndsUnfilteredMultipleBiopsyMutations$pid <- NULL

somatics = mbExonicSomatics %>% 
  select(sampleId, chromosome, position, ref, alt, type, worstCodingEffect, canonicalCodingEffect, hotspot, biallelic, clonality, scope)

mutations = dnds_annotate_somatics(dndsUnfilteredMultipleBiopsyMutations, somatics)

load(file = "~/hmf/RData/processed/genePanel.RData")
genePanel = genePanel %>% filter(martincorena | hmf | cosmicCurated)
mutations = mutations %>% filter(gene %in% genePanel$gene_name, impact != "")

load("~/hmf/RData/processed/excessRates.RData")
load(file = "~/hmf/RData/processed/driverGenes.RData")

mbTsgDrivers = tsg_mutations(mutations %>% filter(gene %in% tsGenes$gene_name))
mbTsgDrivers$driverLikelihood <- tsg_driver_likelihood(mbTsgDrivers, excessTsgRates)
mbTsgDrivers = mbTsgDrivers %>% filter(driverLikelihood > 0)
save(mbTsgDrivers, file = "~/hmf/RData/processed/mbTsgDrivers.RData")

mbOncoDrivers = onco_mutations(mutations %>% filter(gene %in% oncoGenes$gene_name))
mbOncoDrivers$driverLikelihood <- onco_driver_likelihood(mbOncoDrivers, excessOncoRates)
mbOncoDrivers = mbOncoDrivers %>% filter(driverLikelihood > 0)
save(mbOncoDrivers, file = "~/hmf/RData/processed/mbOncoDrivers.RData")
