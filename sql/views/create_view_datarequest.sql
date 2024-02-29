# TODO (KD) Could potentially reuse clinical view once consents are integrated into sample table.

CREATE OR REPLACE VIEW datarequest AS

SELECT * FROM hmfpatients.datarequest_all where allowExternalUseWithoutCheck = 'True';