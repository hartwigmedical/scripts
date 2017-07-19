CREATE TRIGGER trigger_sv_deletion BEFORE DELETE ON structuralVariant
FOR EACH ROW
DELETE FROM sv_annotation2 WHERE sv_id = old.id;