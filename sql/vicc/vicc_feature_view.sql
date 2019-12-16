CREATE OR REPLACE VIEW viccFeature AS

SELECT feature.viccEntryId, feature.name, feature.biomarkerType,
	feature.referenceName, feature.chromosome, feature.start, feature.end, feature.ref, feature.alt,
    feature.provenanceRule, provenanceNames, feature.geneSymbol, feature.entrezId, feature.description,
	sequenceOntology.soid, sequenceOntology.parentSoid, soidHierarchy,
    sequenceOntology.name AS sequenceOntologyName, sequenceOntology.parentName AS sequenceOntologyParentName,
    synonymNames, linkNames, featureInfo.germlineOrSomatic,
    featureAttribute.aminoAcidChange AS attrAminoAcidChange,
    featureAttribute.germline AS attrGermline,
    featureAttribute.partnerGene AS attrPartnerGene,
    featureAttribute.description AS attrDescription,
    featureAttribute.exons AS attrExons,
    featureAttribute.notes AS attrNotes,
    featureAttribute.cosmic AS attrCosmic,
    featureAttribute.effect AS attrEffect,
    featureAttribute.cnvType AS attrCnvType,
    featureAttribute.featureAttributeId AS attrId,
    featureAttribute.cytoband AS attrCytoband,
    featureAttribute.variantType AS attrVariantType,
    featureAttribute.dnaChange AS attrDnaChange,
    featureAttribute.codons AS attrCodons,
    featureAttribute.chromosomeBasedCnv AS attrChromosomeBandCnv,
    featureAttribute.transcript AS attrTranscript,
    featureAttribute.descriptionType AS attrDescriptionType,
    featureAttribute.chromosome AS attrChromosome
FROM feature
LEFT JOIN sequenceOntology ON sequenceOntology.featureId = feature.id
LEFT JOIN featureInfo ON featureInfo.featureId = feature.id
LEFT JOIN featureAttribute ON featureAttribute.featureId = feature.id
LEFT JOIN
	(SELECT sequenceOntologyId, GROUP_CONCAT(hierarchyName SEPARATOR ",") AS soidHierarchy FROM hierarchy GROUP BY 1)
	hierarchies ON hierarchies.sequenceOntologyId = sequenceOntology.id
LEFT JOIN
	(SELECT featureId, GROUP_CONCAT(provenanceName SEPARATOR ",") AS provenanceNames FROM provenance GROUP BY 1)
	provenances ON provenances.featureId = feature.id
LEFT JOIN
	(SELECT featureId, GROUP_CONCAT(synonymName SEPARATOR ",") AS synonymNames FROM synonym GROUP BY 1)
	synonyms ON synonyms.featureId = feature.id
LEFT JOIN
	(SELECT featureId, GROUP_CONCAT(linkName SEPARATOR ",") AS linkNames FROM link GROUP BY 1)
	links ON links.featureId = feature.id;