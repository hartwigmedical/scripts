select association.*, evidence.*, evidenceType.*, publications, variantNames, publicationUrls, phenotype.*, phenotypeType.* from association
inner join evidence on evidence.associationId = association.id
inner join evidenceType on evidenceType.evidenceId = evidence.id
left join phenotype on phenotype.associationId = association.id
left join phenotypeType on phenotypeType.phenotypeId = phenotype.id
left join
	(select evidenceId, GROUP_CONCAT(publication SEPARATOR ",") as publications from evidenceInfo group by 1)
	evidenceInfos on evidenceInfos.evidenceId = evidence.id
left join
	(select associationId, GROUP_CONCAT(urlOfPublication SEPARATOR ",") as publicationUrls from publicationUrl group by 1)
	publicationUrls on publicationUrls.associationId = association.id
left join
	(select associationId, GROUP_CONCAT(variantName SEPARATOR ",") as variantNames from associationVariant group by 1)
	associationVariants on associationVariants.associationId = association.id;