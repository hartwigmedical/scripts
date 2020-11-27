SELECT
DISTINCT patientId AS '#patientId'
FROM
datarequest where preTreatmentsMechanism like "%PI3K%"
				or  preTreatmentsMechanism like "%mTOR%"
				or  preTreatmentsMechanism like "%PARP%"
                or concatenatedTreatmentMechanism like "%PI3K%"
                or concatenatedTreatmentMechanism like "%mTOR%"
                or concatenatedTreatmentMechanism like "%PARP%"
ORDER BY 1;