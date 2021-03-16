SELECT
    DISTINCT patientId AS '#patientId'
FROM
    datarequest
WHERE SUBSTRING_INDEX(doids, ',', 1) IN (1793, 1798, 4074) OR
        SUBSTRING_INDEX(SUBSTRING_INDEX(doids, ',', 2),',',-1) IN (1793, 1798, 4074) OR
        SUBSTRING_INDEX(SUBSTRING_INDEX(doids, ',', 3),',',-1) IN (1793, 1798, 4074)
ORDER BY 1;