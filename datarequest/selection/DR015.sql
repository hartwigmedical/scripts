SELECT
	DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE
    primaryTumorLocation = 'Colon/Rectum'
AND (treatment rlike '%folfiri|folfox|capox|capiri|irinotecan|capecitabine%')
ORDER BY 1;