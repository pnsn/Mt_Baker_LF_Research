-- This query is used on the PNSN archive database to fetch all cataloged seismic events
-- In a bounding box roughly centered on the summit of Mt. Baker
-- Adapted by N. Stevens (ntsteven@uw.edu) from a query composed by A. Hutko (ahutko@uw.edu)
\copy (
    SELECT  e.evid, e.prefor, e.prefmag, e.etype,
            to_timestamp(o.datetime), o.lat, o.lon, o.depth,
            n.magnitude, n.magtype, n.nsta, n.nobs, n.uncertainty 
    FROM 
        event e 
        INNER JOIN origin o ON e.prefor = o.orid 
        INNER JOIN netmag n ON e.prefmag = n.magid 
    WHERE 
        o.lat > 48.60785 AND
        o.lat < 48.93909 AND
        o.lon > -122.0514 AND
        o.lon < -121.48836
    ORDER BY o.datetime; 
) 
TO 'Mount_Bake_Region_Full_Catalog_Event_Summary.csv'
WITH CSV HEADER;