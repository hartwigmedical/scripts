SELECT
	DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE
    primaryTumorLocation = 'Colon/Rectum'
AND (treatment like '%folfiri%' or treatment like '%folfox%' or treatment like '%capox%'
or treatment like '%capiri%' or treatment like '%irinotecan%' or treatment like '%capecitabine%'
or concatenatedTreatmentType like '%Chemo%')
ORDER BY 1;
