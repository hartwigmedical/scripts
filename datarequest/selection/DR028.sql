SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE
    primaryTumorLocation = "Bone/Soft tissue" and primaryTumorType not like "%Chondro%" and primaryTumorType not like "%Osteo%"
ORDER BY 1;