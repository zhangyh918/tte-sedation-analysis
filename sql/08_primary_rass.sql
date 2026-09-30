-- Read-only; primary-cohort IDs and time zero supplied at run time.
WITH e AS (
  SELECT * FROM json_to_recordset($1::json)
    AS q(stay_id integer, time_zero timestamp)
)
SELECT e.stay_id,
  EXTRACT(EPOCH FROM (r.rass_time-e.time_zero))/3600.0 AS time_h,
  r.rass_score
FROM e JOIN analysis._tmp_rass_raw r USING(stay_id)
WHERE r.rass_time >= e.time_zero-INTERVAL '24 hours'
  AND r.rass_time < e.time_zero+INTERVAL '168 hours'
ORDER BY e.stay_id,r.rass_time,r.rass_score;
