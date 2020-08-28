(SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest WHERE hasRNA = 0
ORDER BY 1 limit 4)
UNION
(SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest WHERE hasRNA = 1
ORDER BY 1 limit 4)
UNION
(SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest WHERE patientId in (SELECT patientId FROM datarequest GROUP BY patientId HAVING count(*) > 1)
    ORDER BY 1 limit 4)
UNION
(SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest WHERE hasRNA = 1 AND patientId in (SELECT patientId FROM datarequest GROUP BY patientId HAVING count(*) > 1)
    ORDER BY 1 limit 4);