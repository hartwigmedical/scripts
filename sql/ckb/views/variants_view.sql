CREATE OR REPLACE VIEW variants AS (

SELECT
	createDate, updateDate, variant.fullName, variant, variant.impact, variant.proteinEffect, type, description,
	group_concat(DISTINCT variantPath) AS variantPaths,
	group_concat(DISTINCT pubmedId) AS pubmeds
FROM variant
LEFT JOIN categoryVariantPath ON categoryVariantPath.variantId = variant.id
LEFT JOIN variantReference ON variantReference.variantId = variant.id AND NOT(isnull(pubmedId))
GROUP BY 1,2,3,4,5,6,7,8);