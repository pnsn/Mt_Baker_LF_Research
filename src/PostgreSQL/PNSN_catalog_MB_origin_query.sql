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
--     (x) Origin Error statistics
--   (calc) Calculated distance in km from Mt. Baker summit
--      (n) Network Magnitude
-- Assumes a spherical earth with radius of 6371 km, which is reasonable for local
-- earthquake/reference point scales (a few degrees)

-- NOTE: To output to CSV this needs to be wrapped with \copy (<sql command>) TO <output_file_name> WITH CSV HEADER;

\copy (SELECT 
    e.evid, e.etype, e.selectflag, e.version,
    r.commid, r.remark,
    n.magid, n.magnitude, n.magtype, 
        n.uncertainty, n.nobs, n.nsta,
    o.orid, o.auth, o.subsource,
        to_timestamp(o.datetime), o.lat, o.lon, o.depth,
        (6371. * acos(
        cos(radians(48.7745)) * cos(radians(o.lat)) *
        cos(radians(o.lon) - radians(-121.8172)) + 
        sin(radians(48.7745)) * sin(radians(o.lat))
        )) AS MBS_distance_km,
        o.quality, o.wrms, o.erhor, o.sdep, o.erlat, o.erlon,
        o.bogusflag, o.fdepth, o.fepi, o.ftime,
    x.stt, x.stx, x.sty, x.stz,
           x.sxx, x.sxy, x.sxz,
                  x.syy, x.syz,
                         x.szz,
        x.maglarge, x.azilarge, x.diplarge,
        x.maginter, x.aziinter, x.dipinter,
        x.magsmall, x.azismall, x.dipsmall
FROM event e
    INNER JOIN origin o ON e.prefor = o.orid
    LEFT JOIN netmag n ON e.prefmag = n.magid
    LEFT JOIN remark r ON e.commid = r.commid
    LEFT JOIN origin_error x ON o.orid = x.orid

WHERE (
    6371. * acos(
        cos(radians(48.7745)) * cos(radians(o.lat)) *
        cos(radians(o.lon) - radians(-121.8172)) + 
        sin(radians(48.7745)) * sin(radians(o.lat))
        )
    ) <= 50.
    AND e.selectflag = 1
ORDER BY o.datetime) TO 'MtBaker_50km_radius_origins.csv' WITH CSV HEADER;