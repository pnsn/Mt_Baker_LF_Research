-- script: haversine_query.sql
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
-- (20 km from Mt. Baker's summit in this example)
--
-- Assumes a spherical earth with radius of 6371 km, which is reasonable for local
-- earthquake/reference point scales (a few degrees)


SELECT e.evid, e.etype,
    o.orid, to_timestamp(o.datetime), o.lat, o.lon, o.depth,
    o.wrms, o.erhor, o.sdep, 
FROM event e
    INNER JOIN origin o ON e.prefor = o.orid
WHERE (
    6371. * acos(
        cos(radians(48.7745)) * cos(radians(o.lat)) *
        cos(radians(o.lon) - radians(-121.8172)) + 
        sin(radians(48.7745)) * sin(radians(o.lat))
    )
) <= 20.
ORDER BY o.datetime;