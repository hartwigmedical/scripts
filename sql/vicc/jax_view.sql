CREATE OR REPLACE VIEW jaxView AS

SELECT
    viccEntryId, profileName AS profile, therapyName AS therapy, responseType, evidenceType, approvalStatus, efficacyEvidence,
    source AS indicationSource, name AS indicationName,
    group_concat(url SEPARATOR ",") AS referenceUrls, group_concat(title SEPARATOR ",") AS referenceTitles
FROM jax
INNER JOIN jaxMolecularProfile ON jax.id = jaxMolecularProfile.jaxId
INNER JOIN jaxTherapy ON jax.id = jaxTherapy.jaxId
INNER JOIN jaxIndication ON jax.id = jaxIndication.jaxId
INNER JOIN jaxReference ON jax.id = jaxReference.jaxId
GROUP BY 1,2,3,4,5,6,7,8,9;