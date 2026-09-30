-- Same agents, dose thresholds, units and source tables as the primary cohort.
WITH e AS (
  SELECT * FROM json_to_recordset($1::json)
    AS q(stay_id integer, time_zero timestamp)
)
SELECT e.stay_id,
  EXTRACT(EPOCH FROM (s.infusion_start-e.time_zero))/3600.0 AS start_h,
  EXTRACT(EPOCH FROM (s.infusion_end-e.time_zero))/3600.0 AS end_h,
  s.drug_class,s.rate,s.rateuom
FROM e JOIN analysis._tmp_sedative_raw s USING(stay_id)
WHERE s.infusion_end > s.infusion_start
  AND s.infusion_start < e.time_zero+INTERVAL '168 hours'
  AND s.infusion_end > e.time_zero
ORDER BY e.stay_id,s.infusion_start,s.infusion_end,s.drug_class,s.rate;
