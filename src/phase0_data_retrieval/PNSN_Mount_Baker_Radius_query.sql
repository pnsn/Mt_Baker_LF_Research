\copy SELECT e.evid
FROM event e
    INNER JOIN origin o ON e.prefor = o.orid
WHERE (
    6371. * acos(
        cos(radians(48.7745)) * cos(radians(o.lat)) *
        cos(radians(o.lon) - radians(-121.8172)) + 
        sin(radians(48.7745)) * sin(radians(o.lat))
        )
    ) <= 50.
    AND e.selectflag = 1
ORDER BY o.datetime TO "MtBaker_50km_evids_supp.csv" CSV HEADER;