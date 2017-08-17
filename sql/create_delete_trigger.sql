USE hmfpatients;
DROP TRIGGER IF EXISTS trigger_sv_deletion;
CREATE TRIGGER trigger_sv_deletion AFTER DELETE ON structuralVariant
FOR EACH ROW
DELETE FROM sv_annotation2 WHERE sv_id = old.id;