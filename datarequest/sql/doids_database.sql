SET SESSION group_concat_max_len = @@max_allowed_packet;
SELECT GROUP_CONCAT(distinct REPLACE(doids, ' ' , '' ) SEPARATOR ',') as doids
FROM  datarequest where doids is not null;