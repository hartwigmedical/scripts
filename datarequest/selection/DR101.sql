SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE primaryTumorLocation = 'Prostate' AND sampleId IN
(        SELECT sampleId FROM samplesPanCancerPaper
     UNION
        SELECT sampleId FROM samplesProstatePaperDessel
)
ORDER BY 1;