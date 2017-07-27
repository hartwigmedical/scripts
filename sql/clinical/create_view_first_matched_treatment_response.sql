USE hmfpatients;

DROP VIEW IF EXISTS firstMatchedTreatmentResponse;

CREATE VIEW firstMatchedTreatmentResponse AS 
	SELECT *
	FROM treatmentResponse as tr
	WHERE NOT EXISTS (
		SELECT * 
		FROM treatmentResponse AS tr1 
		WHERE tr1.treatmentId = tr.treatmentId
		AND tr1.responseDate <= tr.responseDate
        AND tr1.id != tr.id
	)
	AND NOT(isnull(treatmentId));