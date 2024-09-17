-- script: PNNS_catalog_MB_main_query.sql
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
-- information for all EVID s selected
--      (e) Event - ID and Classification (etype)
--      (o) Origin - ID, hypocentral parameters, fix-status on solutions, author, 
--      (oe) Origin_Error - covariance matrix upper triangle
--      (n) NetMag - ID, magnitude value, type, uncertainty, and contributing observations
--      (r) remark - ID, comment
-- Assumes a spherical earth with radius of 6371 km, which is reasonable for local
-- earthquake/reference point scales (a few degrees)

\copy (
    SELECT 
        e.evid, e.etype,
        o.orid, to_timestamp(o.datetime), o.lat, o.lon, o.depth,
            o.totalarr, o.ndef, o.nbs, o.quality, s.wrms,
            o.rflag, o.fdepth, o.fepi, o.ftime,
        n.magnitude, n.magtype, n.nsta, n.nobs, n.uncertainty,
        oe.sxx, oe.sxy, oe.sxz, oe.sxt, 
                oe.syy, oe.syz, oe.syt,
                        oe.szz, oe.szt,
                                oe.stt,
        r.remark

    FROM event e
        INNER JOIN origin o ON e.prefor = o.orid
        INNER JOIN origin_error oe on o.orid = oe.orid
        INNER JOIN netmag n on e.prefmag = n.magid
        INNER JOIN remark_origin ro on o.orid = ro.orid
        INNER JOIN remark r on ro.commid = r.commid
    WHERE (
        6371. * acos(
            cos(radians(48.7745)) * cos(radians(o.lat)) *
            cos(radians(o.lon) - radians(-121.8172)) + 
            sin(radians(48.7745)) * sin(radians(o.lat))
            )
        ) <= 20.
        AND e.selectflag= 1

    ORDER BY o.datetime
)