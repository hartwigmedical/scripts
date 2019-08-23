CREATE OR REPLACE VIEW datarequestFiltered AS
SELECT 
    datarequestBase.*
FROM
    datarequestBase
WHERE
    informedConsentDate > '2016-04-20'
        AND purpleQC = 'PASS'
        AND purpleStatus <> 'NO_TUMOR'
        AND purplePurity > 0.195;
