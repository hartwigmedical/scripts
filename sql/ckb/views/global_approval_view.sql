CREATE OR REPLACE VIEW globalApproval AS (

SELECT
    DISTINCT therapyName AS treatment, cancerType, tumorProfile, approvalStatus, approvalAuthority
FROM globalApprovalStatus
INNER JOIN therapy ON therapy.id = globalApprovalStatus.therapyId
INNER JOIN (SELECT DISTINCT ckbIndicationId, name AS cancerType FROM indication)
    AS cancerType ON cancerType.ckbIndicationId = globalApprovalStatus.ckbIndicationId
INNER JOIN (SELECT DISTINCT ckbProfileId, profileName AS tumorProfile FROM ckbEntry)
    AS profile ON profile.ckbProfileId = globalApprovalStatus.ckbProfileId);