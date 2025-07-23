# HIGH LEVEL EVIDENCE

#characteristics
select sourceEvent, type, cutOffType, cutoff, source, sourceUrls, treatment, treatmentApproachesDrugClass, treatmentApproachesTherapy, SUBSTRING_INDEX(indication,"(", 1) as cancerType, REPLACE(SUBSTRING_INDEX(indication,"(", -1),")", "") as doid, evidenceLevel, evidenceDirection
from efficacyEvidence
inner join actionableCharacteristic on actionableCharacteristic.molecularCriteriumId = efficacyEvidence.molecularCriteriumId;

#codons
select  sourceEvent, gene, chromosome, start, end, applicableMutationType, source, sourceUrls, treatment, treatmentApproachesDrugClass, treatmentApproachesTherapy, SUBSTRING_INDEX(indication,"(", 1) as cancerType, REPLACE(SUBSTRING_INDEX(indication,"(", -1),")", "") as doid, evidenceLevel, evidenceDirection
from efficacyEvidence
inner join actionableCodon on actionableCodon.molecularCriteriumId = efficacyEvidence.molecularCriteriumId;

# exons
select  sourceEvent, gene, chromosome, start, end, applicableMutationType, source, sourceUrls, treatment, treatmentApproachesDrugClass, treatmentApproachesTherapy, SUBSTRING_INDEX(indication,"(", 1) as cancerType, REPLACE(SUBSTRING_INDEX(indication,"(", -1),")", "") as doid, evidenceLevel, evidenceDirection
from efficacyEvidence
inner join actionableExon on actionableExon.molecularCriteriumId = efficacyEvidence.molecularCriteriumId;

#fusions
select  sourceEvent, geneUp, minExonUp, maxExonUp, geneDown, minExonDown, maxExonDown source, sourceUrls, treatment, treatmentApproachesDrugClass, treatmentApproachesTherapy, SUBSTRING_INDEX(indication,"(", 1) as cancerType, REPLACE(SUBSTRING_INDEX(indication,"(", -1),")", "") as doid, evidenceLevel, evidenceDirection
from efficacyEvidence
inner join actionableFusion on actionableFusion.molecularCriteriumId = efficacyEvidence.molecularCriteriumId;

#genes
select  sourceEvent, gene, event, source, sourceUrls, treatment, treatmentApproachesDrugClass, treatmentApproachesTherapy, SUBSTRING_INDEX(indication,"(", 1) as cancerType, REPLACE(SUBSTRING_INDEX(indication,"(", -1),")", "") as doid, evidenceLevel, evidenceDirection
from efficacyEvidence
inner join actionableGene on actionableGene.molecularCriteriumId = efficacyEvidence.molecularCriteriumId;

#hla
select  sourceEvent, hlaAllele, treatment, treatmentApproachesDrugClass, treatmentApproachesTherapy, SUBSTRING_INDEX(indication,"(", 1) as cancerType, REPLACE(SUBSTRING_INDEX(indication,"(", -1),")", "") as doid, evidenceLevel, evidenceDirection
from efficacyEvidence
inner join actionableHla on actionableHla.molecularCriteriumId = efficacyEvidence.molecularCriteriumId;


# CLINICAL STUDIES

#characteristics
select sourceEvent, type, cutOffType, cutoff, source, nctId, title, acronym, countriesAndCities, therapynames, genderCriterium, SUBSTRING_INDEX(indications,"(", 1) as cancerType, REPLACE(SUBSTRING_INDEX(indications,"(", -1),")", "") as doid  from actionableTrial
inner join trialMolecularCriterium on trialMolecularCriterium.actionableTrialid = actionableTrial.id
inner join actionableCharacteristic on actionableCharacteristic.molecularCriteriumId = trialMolecularCriterium.molecularCriteriumId;

#codons
#exons
#fusions
#genes
#hla
