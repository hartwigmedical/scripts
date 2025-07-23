# HIGH LEVEL EVIDENCE

#characteristics
select sourceEvent, type, cutOffType, cutoff, source, sourceUrls, treatment, treatmentApproachesDrugClass, treatmentApproachesTherapy, SUBSTRING_INDEX(indication,"(", 1) as cancerType, REPLACE(SUBSTRING_INDEX(indication,"(", -1),")", "") as doid, evidenceLevel, evidenceDirection
from efficacyEvidence
inner join actionableCharacteristic on actionableCharacteristic.molecularCriteriumId = efficacyEvidence.molecularCriteriumId;

#codons
#exons
#fusions
#genes
#hla


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
