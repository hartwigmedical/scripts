CREATE OR REPLACE VIEW datarequest_CB01 AS

SELECT * FROM hmfpatients.datarequest_all where allowExternalUseWithCheck = 'true';