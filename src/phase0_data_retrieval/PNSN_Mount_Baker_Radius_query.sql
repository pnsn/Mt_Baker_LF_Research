SELECT e.evid, e.etype, "uw"e.evid
    FROM event e
        INNER JOIN origin o ON e.prefor = o.orid
    WHERE (
        6371. * acos(
            cos(radians(48.7745)) * cos(radians(o.lat)) *
            cos(radians(o.lon) - radians(-121.8172)) + 
            sin(radians(48.7745)) * sin(radians(o.lat))
            )
        ) <= 50.
        AND e.selectflag = 1 AND e.evid > 62061782
    ORDER BY o.datetime; TO :ocsv CSV HEADER;