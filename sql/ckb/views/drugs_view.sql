CREATE OR REPLACE VIEW drugs AS (

SELECT
    drug.createDate, drugName, tradeName, group_concat(DISTINCT term) AS terms, group_concat(DISTINCT synonym) AS synonyms,
    group_concat(DISTINCT drugClass) AS drugClasses,
    casRegistryNum, ncitId, description, group_concat(DISTINCT pubmedId) AS pubmeds
FROM drug
LEFT JOIN drugClass ON drugClass.drugId = drug.id
LEFT JOIN drugTerm ON drugTerm.drugId = drug.id
LEFT JOIN drugSynonym ON drugSynonym.drugId = drug.id
LEFT JOIN drugReference ON drugReference.drugId = drug.id AND NOT(isnull(pubmedId))
GROUP BY 1,2,3,7,8,9);