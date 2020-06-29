SELECT DISTINCT patientId as '#patientId'
FROM (
    SELECT DISTINCT sampleId
        FROM datarequest
        WHERE hasRNA=1
    UNION
    SELECT DISTINCT(svBreakend.sampleId)
        FROM svBreakend
        WHERE gene = 'FGFR2')
    AS tmp
 INNER JOIN datarequest ON tmp.sampleId = datarequest.sampleId;