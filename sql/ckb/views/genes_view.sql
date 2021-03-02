CREATE OR REPLACE VIEW genes AS (

SELECT DISTINCT
	createDate, updateDate, geneSymbol,
    geneRole, entrezId, chromosome, mapLocation, canonicalTranscript, description,
    group_concat(DISTINCT pubmedId) as pubmeds
FROM gene
LEFT JOIN geneReference on geneReference.geneId = gene.id AND NOT(isnull(pubmedId))
GROUP BY 1,2,3,4,5,6,7,8,9);