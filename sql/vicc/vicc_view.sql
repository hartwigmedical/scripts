CREATE OR REPLACE VIEW vicc AS

SELECT
	viccEntry.id, viccEntry.source,
	geneNames.genes, associationVariants.variants, featureNames.features,
    geneIdentifiers.geneSymbols, geneIdentifiers.entrezIds, geneIdentifiers.ensemblGeneIds,
	association.evidenceLevel, association.evidenceLabel, association.responseType, association.drugLabels,
    association.sourceLink, association.description, association.oncogenic,
    evidence.description AS evidenceDescription, publications, publicationUrls,
    evidenceType.sourceName AS evidenceSourceName, evidenceType.idEvidenceType AS evidenceSourceIdentifier,
    phenotype.description AS phenotypeDescription, phenotype.family AS phenotypeFamily, phenotype.idPhenotype AS phenotypeIdentifier,
    phenotypeType.source AS phenotypeSource, phenotypeType.term AS phenotypeTerm, phenotypeType.idPhenotypeType AS phenotypeTypeIdentifier,
    tagNames.tags, devTagNames.devTags
FROM viccEntry
INNER JOIN
	(SELECT viccEntryId, GROUP_CONCAT(geneName SEPARATOR ",") AS genes FROM gene GROUP BY 1)
	geneNames ON geneNames.viccEntryId = viccEntry.id
LEFT JOIN
	(select viccEntryId, GROUP_CONCAT(symbol SEPARATOR ",") AS geneSymbols, GROUP_CONCAT(entrezId SEPARATOR ",") AS entrezIds,
		GROUP_CONCAT(ensemblGeneId SEPARATOR ",") AS ensemblGeneIds FROM geneIdentifier GROUP BY 1)
	geneIdentifiers ON geneIdentifiers.viccEntryId = viccEntry.id
LEFT JOIN
	(select viccEntryId, GROUP_CONCAT(nameOfFeature SEPARATOR ",") AS features FROM featureName GROUP BY 1)
	featureNames ON featureNames.viccEntryId = viccEntry.id
LEFT JOIN
	(select viccEntryId, GROUP_CONCAT(tagName SEPARATOR ",") AS tags FROM tag GROUP BY 1)
	tagNames ON tagNames.viccEntryId = viccEntry.id
LEFT JOIN
	(select viccEntryId, GROUP_CONCAT(devTagName SEPARATOR ",") AS devTags FROM devTag GROUP BY 1)
	devTagNames ON devTagNames.viccEntryId = viccEntry.id
INNER JOIN association ON association.viccEntryId = viccEntry.id
INNER JOIN evidence ON evidence.associationId = association.id
INNER JOIN evidenceType ON evidenceType.evidenceId = evidence.id
LEFT JOIN phenotype ON phenotype.associationId = association.id
LEFT JOIN phenotypeType ON phenotypeType.phenotypeId = phenotype.id
LEFT JOIN
	(SELECT evidenceId, GROUP_CONCAT(publication SEPARATOR ",") AS publications FROM evidenceInfo GROUP BY 1)
	evidenceInfos ON evidenceInfos.evidenceId = evidence.id
LEFT JOIN
	(SELECT associationId, GROUP_CONCAT(urlOfPublication SEPARATOR ",") AS publicationUrls FROM publicationUrl GROUP BY 1)
	publicationUrls ON publicationUrls.associationId = association.id
LEFT JOIN
	(SELECT associationId, GROUP_CONCAT(variantName SEPARATOR ",") AS variants FROM associationVariant GROUP BY 1)
	associationVariants ON associationVariants.associationId = association.id;