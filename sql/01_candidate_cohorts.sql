-- READ ONLY. No CREATE/DROP/UPDATE against the database.
-- _tmp_exclusion_chain is the source cohort table (94,458 ICU stays).
-- Retain the source neurological-exclusion flag for exact attribution; do not silently
-- substitute a differently parenthesised historical ICD query.
WITH base AS (
    SELECT * FROM analysis._tmp_exclusion_chain
    WHERE NOT excl_not_first_icu AND NOT excl_underage AND NOT excl_no_mv
), timings AS (
    SELECT b.*, s.first_sedative_start, r.first_rass_time
    FROM base b
    LEFT JOIN LATERAL (
        SELECT MIN(infusion_start) AS first_sedative_start
        FROM analysis._tmp_sedative_raw x
        WHERE x.stay_id=b.stay_id AND x.infusion_start>=b.first_mv_start
          AND x.infusion_start<=b.first_mv_start+INTERVAL '12 hours'
          AND (b.icu_outtime IS NULL OR x.infusion_start<=b.icu_outtime)
          AND x.infusion_end>x.infusion_start AND x.rate>0
    ) s ON TRUE
    LEFT JOIN LATERAL (
        SELECT MIN(rass_time) AS first_rass_time
        FROM analysis._tmp_rass_raw x
        WHERE x.stay_id=b.stay_id AND x.rass_time>=b.first_mv_start
          AND x.rass_time<=b.first_mv_start+INTERVAL '24 hours'
          AND (b.icu_outtime IS NULL OR x.rass_time<=b.icu_outtime)
          AND x.rass_score BETWEEN -5 AND 4
    ) r ON TRUE
), timezero AS (
    SELECT *, CASE WHEN first_sedative_start IS NOT NULL AND first_rass_time IS NOT NULL
        THEN GREATEST(first_mv_start,first_sedative_start,first_rass_time) END AS time_zero
    FROM timings
), service_at_zero AS (
    SELECT DISTINCT ON (t.stay_id) t.stay_id,s.curr_service,s.transfertime
    FROM timezero t JOIN mimiciv_hosp.services s
      ON s.hadm_id=t.hadm_id AND s.transfertime<=t.time_zero
    ORDER BY t.stay_id,s.transfertime DESC,s.curr_service
)
SELECT t.*, (old.stay_id IS NOT NULL) AS original_5086,
    old.time_zero AS original_time_zero,
    old.first_rass_time AS original_first_rass_time,
    old.first_sedative_start AS original_first_sedative_start,
    v.on_mv_at_t0, srv.curr_service AS service_at_t0,
    srv.transfertime AS service_timestamp,
    (t.hospital_deathtime IS NULL OR t.hospital_deathtime>t.time_zero) AS alive_at_t0,
    (t.icu_intime<=t.time_zero AND (t.icu_outtime IS NULL OR t.icu_outtime>t.time_zero)) AS in_icu_at_t0,
    p.anchor_year_group,
    EXTRACT(EPOCH FROM (t.time_zero-t.first_mv_start))/3600.0 AS hours_mv_to_t0,
    EXTRACT(EPOCH FROM (t.first_rass_time-t.first_sedative_start))/3600.0 AS first_rass_delay_from_sed_hours
FROM timezero t
LEFT JOIN analysis.tte_eligible old USING (stay_id)
JOIN mimiciv_hosp.patients p ON p.subject_id=t.subject_id
LEFT JOIN LATERAL (
    SELECT EXISTS (SELECT 1 FROM analysis._tmp_mv_raw m
        WHERE m.stay_id=t.stay_id AND m.mv_start<=t.time_zero
          AND m.mv_end>t.time_zero) AS on_mv_at_t0
) v ON TRUE
LEFT JOIN service_at_zero srv ON srv.stay_id=t.stay_id
ORDER BY t.stay_id;
