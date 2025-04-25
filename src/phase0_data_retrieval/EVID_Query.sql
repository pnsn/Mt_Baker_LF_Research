\copy (
    SELECT e.evid,o.orid,o.lat,o.lon,o.depth,o.datetime,e.etype,
        wheres.separation_km(
            o.lat::double precision,o.lon::double precision,
            48.7745::double precision,-121.8172::double precision)
    FROM event e INNER JOIN origin o ON e.prefor=o.orid 
    WHERE wheres.separation_km(
            o.lat::double precision,o.lon::double precision,
            48.7745::double precision,-121.8172::double precision) <= 30. 
        AND e.selectflag=1 ORDER BY o.datetime
    ) 
TO 'MtBakerQuery.csv' CSV HEADER;