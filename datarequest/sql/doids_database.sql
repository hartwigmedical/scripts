SELECT GROUP_CONCAT(distinct REPLACE(doids, ' ' , '' ) SEPARATOR ', ') as doids
FROM  datarequest where doids is not null;