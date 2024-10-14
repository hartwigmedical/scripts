USE actin_paper;

DROP TABLE IF EXISTS eligibleCohorts_addition;
CREATE TABLE eligibleCohorts_addition (
patientId varchar(50) NOT NULL,
trialId varchar(50) NOT NULL,
trialAcronym varchar(50) NOT NULL,
cohortDescription varchar(500) NOT NULL,
event varchar(50)
);

DROP TABLE IF EXISTS variant_addition;
CREATE TABLE variant_addition (
sampleId varchar(50) NOT NULL,
isReportable boolean NOT NULL,
gene varchar(50) NOT NULL,
event varchar(500) NOT NULL,
driverLikelihood varchar(50) NOT NULL
);

DROP TABLE IF EXISTS copyNumber_addition;
CREATE TABLE copyNumber_addition (
sampleId varchar(50) NOT NULL,
isReportable boolean NOT NULL,
gene varchar(50) NOT NULL,
event varchar(500) NOT NULL,
driverLikelihood varchar(50) NOT NULL
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
	SELECT sampleId, event, gene, driverLikelihood, "a_variant" as category
    FROM variant_addition
    WHERE isReportable
		UNION
	SELECT sampleId, event, gene, driverLikelihood, "b_copy_number" as category
    FROM copyNumber
    WHERE isReportable
		UNION
	SELECT sampleId, event, gene, driverLikelihood, "b_copy_number" as category
    FROM copyNumber_addition
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
ORDER BY 1,5,4,3,2
);