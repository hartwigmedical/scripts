CREATE OR REPLACE VIEW datarequest_all AS

select * from datarequest
union all
select * from datarequest_CB01;