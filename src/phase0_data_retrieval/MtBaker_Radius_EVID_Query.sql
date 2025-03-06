\copy SELECT e.evid
    FROM event e
        INNER JOIN origin o ON e.prefor = o.orid
    WHERE (
        6371. * acos(
            cos(radians(48.7745)) * cos(radians(o.lat)) *
            cos(radians(o.lon) - radians(-121.8172)) + 
            sin(radians(48.7745)) * sin(radians(o.lat))
            )
        ) <= :radius
        AND e.selectflag = 1 AND e.evid > :levid
    ORDER BY o.datetime TO :ocsv CSV HEADER;