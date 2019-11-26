select viccEntry.*, geneNames, symbols, entrezIds, ensemblGeneIds, featureNames, tagNames, devTagNames from viccEntry
inner join
	(select viccEntryId, GROUP_CONCAT(geneName SEPARATOR ",") as geneNames from gene group by 1)
	geneNames on geneNames.viccEntryId = viccEntry.id
left join
	(select viccEntryId, GROUP_CONCAT(symbol SEPARATOR ",") as symbols, GROUP_CONCAT(entrezId SEPARATOR ",") as entrezIds,
		GROUP_CONCAT(ensemblGeneId SEPARATOR ",") as ensemblGeneIds from geneIdentifier group by 1)
	geneIdentifiers on geneIdentifiers.viccEntryId = viccEntry.id
left join
	(select viccEntryId, GROUP_CONCAT(nameOfFeature SEPARATOR ",") as featureNames from featureName group by 1)
	featureNames on featureNames.viccEntryId = viccEntry.id
left join
	(select viccEntryId, GROUP_CONCAT(tagName SEPARATOR ",") as tagNames from tag group by 1)
	tags on tags.viccEntryId = viccEntry.id
left join
	(select viccEntryId, GROUP_CONCAT(devTagName SEPARATOR ",") as devTagNames from devTag group by 1)
	devTags on devTags.viccEntryId = viccEntry.id;