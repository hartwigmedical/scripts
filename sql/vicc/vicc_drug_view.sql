CREATE OR REPLACE VIEW viccDrug AS

SELECT
	association.viccEntryId,
	environmentalContext.term, environmentalContext.description, environmentalContext.source, environmentalContext.usanStem,
    environmentalContext.toxicity, environmentalContext.idEnvironmentalContext AS identifier, approvedCountries,
    taxonomy.kingdom, taxonomy.directParent, taxonomy.class, taxonomy.subClass, taxonomy.superClass
FROM association
INNER JOIN environmentalContext ON environmentalContext.associationId = association.id
LEFT JOIN taxonomy ON environmentalContext.id = taxonomy.environmentalContextId
LEFT JOIN
	(SELECT environmentalContextId, GROUP_CONCAT(approvedCountryName SEPARATOR ",") AS approvedCountries FROM approvedCountry GROUP BY 1)
	approvedCountries ON approvedCountries.environmentalContextId = environmentalContext.id;