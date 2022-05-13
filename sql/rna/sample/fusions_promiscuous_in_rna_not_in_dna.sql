SELECT primaryTumorLocation, cuppaTumorLocation, b.reportedType, dnaFusionCount, svFusion.reported AS reportedDna, a.*
FROM
(SELECT * FROM rnaFusion WHERE
(name LIKE 'SLC45A3%' OR
name LIKE 'TMPRSS2%' OR
name LIKE 'FGFR3%' OR
name LIKE 'FGFR1%' OR
name LIKE 'FGFR2%' OR
name LIKE 'HMGA2%' OR
name LIKE 'BCOR%' OR
name LIKE 'BCR%' OR
name LIKE 'CD74%' OR
name LIKE 'CLTC%' OR
name LIKE 'ETV6%' OR
name LIKE 'EWSR1%' OR
name LIKE 'FUS%' OR
name LIKE 'KAT6A%' OR
name LIKE 'KIF5B%' OR
name LIKE 'KMT2A%' OR
name LIKE 'NPM1%' OR
name LIKE 'NUP98%' OR
name LIKE 'PAX3%' OR
name LIKE 'PAX5%' OR
name LIKE 'PAX8%' OR
name LIKE 'RUNX1%' OR
name LIKE 'SS18%' OR
name LIKE 'STRN%' OR
name LIKE 'TFG%' OR
name LIKE 'TPM3%' OR
name LIKE 'TPR%' OR
name LIKE 'YAP1%' OR
name LIKE 'YWHAE%' OR
name LIKE 'CSF1%' OR
name LIKE 'SRF%')
AND RIGHT(name,1) !="_" AND (svType IN ('BND','INV','INS') OR abs(positionUp-positionDown)>1e6)
OR
(name LIKE '%BRAF' OR
name LIKE '%ETV4' OR
name LIKE '%RET' OR
name LIKE '%ROS1' OR
name LIKE '%ALK' OR
name LIKE '%ETV1' OR
name LIKE '%MET' OR
name LIKE '%NRG1' OR
name LIKE '%NTRK1' OR
name LIKE '%NTRK2' OR
name LIKE '%NTRK3' OR
name LIKE '%ABL1' OR
name LIKE '%CRLF2' OR
name LIKE '%ERG' OR
name LIKE '%ETV6' OR
name LIKE '%FGFR1' OR
name LIKE '%JAK2' OR
name LIKE '%MAML2' OR
name LIKE '%MYC' OR
name LIKE '%NCOA2' OR
name LIKE '%NFIB' OR
name LIKE '%NR4A3' OR
name LIKE '%NUTM1' OR
name LIKE '%PDGFRA' OR
name LIKE '%PDGFRB' OR
name LIKE '%PHF1' OR
name LIKE '%PLAG1' OR
name LIKE '%RAF1' OR
name LIKE '%RARA' OR
name LIKE '%TFE3' OR
name LIKE '%TFEB' OR
name LIKE '%USP6' OR
name LIKE '%NCOA1') AND LEFT(NAME,1) !="_" AND (svType IN ('BND','INV','INS') OR abs(positionUp-positionDown)>1e6)
)
AS a
LEFT JOIN (SELECT DISTINCT NAME, reportedType, count(*) AS dnaFusionCount FROM svFusion GROUP BY 1,2) AS b
    ON a.name=b.name
LEFT JOIN svFusion
    ON a.sampleId=svFusion.sampleId AND b.name=svFusion.name
LEFT JOIN clinical
    ON a.sampleId=clinical.sampleId
LEFT JOIN cuppa
    ON a.sampleId=cuppa.sampleId
WHERE (svFusion.sampleId IS NULL OR svFusion.reported != 1) AND a.sampleId IN ('XXX') AND (b.reportedType IS NULL OR b.reportedType != 'KNOWN_PAIR')
ORDER BY b.reportedType, a.name, a.sampleId;