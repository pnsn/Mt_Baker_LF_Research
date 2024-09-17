-- script: PNNS_catalog_MB_origin_query.sql
-- author: Nathan T. Stevens (ntsteven@uw.edu) + ChatGPT (see attribution)
-- org: Pacific Northwest Seismic Network
-- license: GPLv3
--
-- Attribution: This script was based on a prompt to ChatGPT Auto on 17. Sept. 2024
-- and independently tested/modified by the author to verify expected performance
--
-- Purpose:
-- Use the Haversine function to get kilometer distance estimates between
-- input lat/lon points from a PostgreSQL table and a reference lat/lon
-- (the summit of Mt. Baker in this example) and filter by a reference distance
-- (20 km from Mt. Baker's summit in this example) and then provide the following
--      (e) Event ID and metadata
--      (r) Comment from analyst/import
--      (o) Origin best-fit solution
--     (oe) Origin Error statistics
--   (calc) Calculated distance in km from Mt. Baker summit
--      (n) Network Magnitude
-- Assumes a spherical earth with radius of 6371 km, which is reasonable for local
-- earthquake/reference point scales (a few degrees)

-- NOTE: To output to CSV this needs to be wrapped with \copy (<sql command>) TO <output_file_name> WITH CSV HEADER;

SELECT 
    e.evid, e.etype, r.remark, o.auth,
    (6371. * acos(
        cos(radians(48.7745)) * cos(radians(o.lat)) *
        cos(radians(o.lon) - radians(-121.8172)) + 
        sin(radians(48.7745)) * sin(radians(o.lat))
        )) AS calculated_distance_km,
    to_timestamp(o.datetime), o.lat, o.lon, o.depth,
        o.totalarr, o.ndef, o.nbs, o.quality, o.wrms,
        o.rflag, o.fdepth, o.fepi, o.ftime, o.bogusflag,
    n.magnitude, n.magtype, n.nsta, n.nobs, n.uncertainty,
    oe.*,
    e.subsource, e.version, e.lddate,
    
FROM event e
    INNER JOIN origin o ON e.prefor = o.orid
    INNER JOIN netmag n ON e.prefmag = n.magid
    INNER JOIN remark r ON e.commid = r.commid
    INNER JOIN origin_error oe ON o.orid = oe.orid

WHERE (
    6371. * acos(
        cos(radians(48.7745)) * cos(radians(o.lat)) *
        cos(radians(o.lon) - radians(-121.8172)) + 
        sin(radians(48.7745)) * sin(radians(o.lat))
        )
    ) <= 50.
ORDER BY o.datetime