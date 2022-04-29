SELECT DISTINCT sampleId, trialId, trialAcronym, cohortDescription
FROM trialEvaluation
WHERE sampleId IN ('XXX')
AND ((isEligibleTrial AND NOT trialHasCohorts AND trialOpen) OR (isEligibleTrial AND isEligibleCohort AND cohortOpen AND NOT cohortBlacklist));