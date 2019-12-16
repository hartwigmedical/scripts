SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE hmfPatientId IN (SELECT hmfPatientId FROM datarequest GROUP BY 1 HAVING count(*) >= 2)
ORDER BY 1;