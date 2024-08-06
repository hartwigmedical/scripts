USE actin_paper;

DROP TABLE IF EXISTS eligibleCohorts_addition;
CREATE TABLE eligibleCohorts_addition (
patientId varchar(50) NOT NULL,
trialId varchar(50) NOT NULL,
trialAcronym varchar(50) NOT NULL,
cohortDescription varchar(500) NOT NULL,
event varchar(50)
);

DROP TABLE IF EXISTS paperSamples;
CREATE TABLE paperSamples (
sampleId varchar(50) NOT NULL,
patientId varchar(50) NOT NULL
);

CREATE OR REPLACE VIEW molecularDrivers
AS (
SELECT * FROM (
	SELECT sampleId, event, gene, driverLikelihood, "a_variant" as category
    FROM variant
    WHERE isReportable
		UNION
	SELECT sampleId, event, gene, driverLikelihood, "b_copy_number" as category
    FROM copyNumber
    WHERE isReportable
		UNION
	SELECT sampleId, event, gene, driverLikelihood, "e_hom_disruptions" as category
    FROM homozygousDisruption
    WHERE isReportable
		UNION
	SELECT sampleId, event, gene, driverLikelihood, "f_het_disruptions" as category
    FROM disruption
    WHERE isReportable
		UNION
	SELECT sampleId, event, "NA" as gene, driverLikelihood, "c_fusions" as category
    FROM fusion
    WHERE isReportable
		UNION
	SELECT sampleId, event, "NA" as gene, driverLikelihood, "d_viruses" as category
    FROM virus
    WHERE isReportable)
    AS a
ORDER BY 1,4,3,2
);