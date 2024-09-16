-- Query composed by Alex Hutko (2024 09 04)
-- Pulls all event types 'le' and 'lf' from the PNSN postprocessing archive
-- Provides event, origin, magnitude, and arrival information for each associated arrival to one of these events preferred solutions
\copy ( 
SELECT 
    e.evid, e.etype, to_timestamp(o.datetime),
    o.lat, o.lon, o.depth, o.distance, o.wrms, o.erhor, o.sdep, o.rflag,
    to_timestamp(a.datetime), a.net, a.sta, a.location, a.seedchan, a.iphase, a.quality, a.deltim,
    n.magnitude, n.magtype, n.nsta, n.nobs, n.uncertainty
FROM event e
    INNER JOIN origin o ON e.prefor = o.orid
    INNER JOIN netmag n ON e.prefmag = n.magid
    INNER JOIN assocaro ac ON ac.orid = o.orid
    INNER JOIN arrival a ON a.arid = ac.arid WHERE (e.etype = 'le' or e.etype = 'lf')
        AND e.selectflag = 1
ORDER BY a.datetime)
TO 'le_and_lf_full_catalog_query.csv' WITH CSV HEADER; -- Nate added this line to output to file