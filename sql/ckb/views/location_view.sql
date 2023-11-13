CREATE OR REPLACE VIEW locations AS (

SELECT DISTINCT
    nctId, group_concat(distinct country) as countries from location
    group by clinicalTrialId);