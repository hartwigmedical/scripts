CREATE OR REPLACE VIEW drugs AS (

SELECT drug.createDate, drugName, tradeName, casRegistryNum, ncitId, group_concat(DISTINCT drugClass) AS drugClasses
FROM drug LEFT JOIN drugClass ON drugClass.drugId = drug.id
GROUP BY 1,2,3,4,5);