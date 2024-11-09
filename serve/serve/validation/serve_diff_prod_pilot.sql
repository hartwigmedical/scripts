#eventInterpretation
SELECT version, sourceEvent, interpretedGene, interpretedEvent, interpretedEventType
FROM (
(SELECT 'OnlyInTruth' AS version, sourceEvent, interpretedGene, interpretedEvent, interpretedEventType FROM serve_production.eventInterpretation
where source = "CKB")
UNION
(SELECT 'OnlyInNew' AS version,
sourceEvent, interpretedGene, interpretedEvent, interpretedEventType FROM serve_pilot.eventInterpretation)) as a
group by  sourceEvent, interpretedGene, interpretedEvent, interpretedEventType
having count(version) = 1
order by version, sourceEvent, interpretedGene, interpretedEvent, interpretedEventType desc;

#ActionableCharacteristics
SELECT version, type,cutoffType,cutoff,sourceEvent,sourceUrls,treatment,treatmentApproachesDrugClass, applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction,evidenceUrls
FROM (
(SELECT 'OnlyInTruth' AS version, type,cutoffType,cutoff,sourceEvent,sourceUrls,treatment,case when sourceTreatmentApproach = " " then null else sourceTreatmentApproach end as treatmentApproachesDrugClass,applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction,evidenceUrls FROM serve_production.actionableCharacteristic
where source = "CKB")
UNION
(SELECT 'OnlyInNew' AS version,
type,cutoffType,cutoff,sourceEvent,sourceUrls,treatment,treatmentApproachesDrugClass, applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction,evidenceUrls FROM serve_pilot.actionableCharacteristic where source = "CKB_EVIDENCE")) as a
group by  type,cutoffType,cutoff,sourceEvent,sourceUrls,treatment,treatmentApproachesDrugClass, applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction,evidenceUrls
having count(version) = 1
order by type,cutoffType,cutoff,sourceEvent,sourceUrls,treatment,treatmentApproachesDrugClass, applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction,evidenceUrls desc;

#ActionableCodon
SELECT version, gene,chromosome,start,end,applicableMutationType,sourceEvent,sourceUrls,treatment,treatmentApproachesDrugClass,applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction, evidenceUrls
FROM (
(SELECT 'OnlyInTruth' AS version, gene,chromosome,start,end,applicableMutationType,sourceEvent,sourceUrls,treatment,case when sourceTreatmentApproach = " " then null else sourceTreatmentApproach end as treatmentApproachesDrugClass,applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction, evidenceUrls FROM serve_production.actionableCodon
where source = "CKB")
UNION
(SELECT 'OnlyInNew' AS version,
gene,chromosome,start,end,applicableMutationType,sourceEvent,sourceUrls,treatment,treatmentApproachesDrugClass,applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction, evidenceUrls FROM serve_pilot.actionableCodon where source = "CKB_EVIDENCE")) as a
group by gene,chromosome,start,end,applicableMutationType,sourceEvent,sourceUrls,treatment,treatmentApproachesDrugClass,applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction, evidenceUrls
having count(version) = 1
order by gene,chromosome,start,end,applicableMutationType,sourceEvent,sourceUrls,treatment,treatmentApproachesDrugClass,applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction, evidenceUrls desc;

#actionableExon
SELECT version, gene,chromosome,start,end,applicableMutationType,sourceEvent,sourceUrls,treatment,treatmentApproachesDrugClass,applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction, evidenceUrls
FROM (
(SELECT 'OnlyInTruth' AS version, gene,chromosome,start,end,applicableMutationType,sourceEvent,sourceUrls,treatment,case when sourceTreatmentApproach = " " then null else sourceTreatmentApproach end as treatmentApproachesDrugClass, applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction, evidenceUrls FROM serve_production.actionableExon
where source = "CKB")
UNION
(SELECT 'OnlyInNew' AS version,
gene,chromosome,start,end,applicableMutationType,sourceEvent,sourceUrls,treatment,treatmentApproachesDrugClass,applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction, evidenceUrls FROM serve_pilot.actionableExon where source = "CKB_EVIDENCE")) as a
group by gene,chromosome,start,end,applicableMutationType,sourceEvent,sourceUrls,treatment,treatmentApproachesDrugClass,applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction, evidenceUrls
having count(version) = 1
order by gene,chromosome,start,end,applicableMutationType,sourceEvent,sourceUrls,treatment,treatmentApproachesDrugClass,applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction, evidenceUrls desc;

#actionableFusion
SELECT version, geneUp,minExonUp,maxExonUp,geneDown,minExonDown,maxExonDown,sourceEvent,sourceUrls,treatment,treatmentApproachesDrugClass,applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction,evidenceUrls
FROM (
(SELECT 'OnlyInTruth' AS version, geneUp,minExonUp,maxExonUp,geneDown,minExonDown,maxExonDown,sourceEvent,sourceUrls,treatment,case when sourceTreatmentApproach = " " then null else sourceTreatmentApproach end as treatmentApproachesDrugClass, applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction,evidenceUrls FROM serve_production.actionableFusion
where source = "CKB")
UNION
(SELECT 'OnlyInNew' AS version,
geneUp,minExonUp,maxExonUp,geneDown,minExonDown,maxExonDown,sourceEvent,sourceUrls,treatment,treatmentApproachesDrugClass,applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction,evidenceUrls FROM serve_pilot.actionableFusion where source = "CKB_EVIDENCE")) as a
group by geneUp,minExonUp,maxExonUp,geneDown,minExonDown,maxExonDown,sourceEvent,sourceUrls,treatment,treatmentApproachesDrugClass,applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction,evidenceUrls
having count(version) = 1
order by geneUp,minExonUp,maxExonUp,geneDown,minExonDown,maxExonDown,sourceEvent,sourceUrls,treatment,treatmentApproachesDrugClass,applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction,evidenceUrls	desc;

#actionableGene
SELECT version, gene,event,sourceEvent,sourceUrls, treatment,treatmentApproachesDrugClass,applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction, evidenceUrls
FROM (
(SELECT 'OnlyInTruth' AS version, gene,event,sourceEvent,sourceUrls, treatment,case when sourceTreatmentApproach = " " then null else sourceTreatmentApproach end as treatmentApproachesDrugClass, applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction, evidenceUrls FROM serve_production.actionableGene
where source = "CKB")
UNION
(SELECT 'OnlyInNew' AS version,
gene,event,sourceEvent,sourceUrls, treatment,treatmentApproachesDrugClass,applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction, evidenceUrls	 FROM serve_pilot.actionableGene where source = "CKB_EVIDENCE")) as a
group by gene,event,sourceEvent,sourceUrls, treatment,treatmentApproachesDrugClass,applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction, evidenceUrls
having count(version) = 1
order by gene,event,sourceEvent,sourceUrls, treatment,treatmentApproachesDrugClass,applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction, evidenceUrls desc;

#actionableHla
SELECT version, hlaAllele,sourceEvent,sourceUrls,treatment,treatmentApproachesDrugClass,applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction,evidenceUrls
FROM (
(SELECT 'OnlyInTruth' AS version, hlaAllele,sourceEvent,sourceUrls,treatment,case when sourceTreatmentApproach = " " then null else sourceTreatmentApproach end as treatmentApproachesDrugClass, applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction,evidenceUrls FROM serve_production.actionableHla
where source = "CKB")
UNION
(SELECT 'OnlyInNew' AS version,
hlaAllele,sourceEvent,sourceUrls,treatment,treatmentApproachesDrugClass,applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction,evidenceUrls FROM serve_pilot.actionableHla where source = "CKB_EVIDENCE")) as a
group by hlaAllele,sourceEvent,sourceUrls,treatment,treatmentApproachesDrugClass,applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction,evidenceUrls
having count(version) = 1
order by hlaAllele,sourceEvent,sourceUrls,treatment,treatmentApproachesDrugClass,applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction,evidenceUrls desc;

#actionableHotspot
SELECT version, gene,chromosome,position,ref,alt,sourceEvent,sourceUrls,treatment,treatmentApproachesDrugClass,applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction, evidenceUrls
FROM (
(SELECT 'OnlyInTruth' AS version, gene,chromosome,position,ref,alt,sourceEvent,sourceUrls,treatment,case when sourceTreatmentApproach = " " then null else sourceTreatmentApproach end as treatmentApproachesDrugClass, applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction, evidenceUrls FROM serve_production.actionableHotspot
where source = "CKB")
UNION
(SELECT 'OnlyInNew' AS version,
gene,chromosome,position,ref,alt,sourceEvent,sourceUrls,treatment,treatmentApproachesDrugClass,applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction, evidenceUrls FROM serve_pilot.actionableHotspot where source = "CKB_EVIDENCE")) as a
group by gene,chromosome,position,ref,alt,sourceEvent,sourceUrls,treatment,treatmentApproachesDrugClass,applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction, evidenceUrls
having count(version) = 1
order by gene,chromosome,position,ref,alt,sourceEvent,sourceUrls,treatment,treatmentApproachesDrugClass,applicableCancerType,applicableDoid,blacklistCancerTypes,level,direction, evidenceUrls desc;