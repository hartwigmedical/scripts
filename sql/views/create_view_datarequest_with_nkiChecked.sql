CREATE OR REPLACE VIEW datarequest_nkiChecked AS

SELECT * FROM hmfpatients.datarequest_all where AllowExternalUseIRBchecked = 1;