SELECT DISTINCT(patientId) AS '#patientId'
FROM datarequest
WHERE sampleId IN
(SELECT DISTINCT(svBreakend.sampleId) FROM svBreakend
WHERE gene = 'FGFR2');