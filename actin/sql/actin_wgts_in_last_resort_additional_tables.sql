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