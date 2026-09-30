-- sensitivity cohort without future duration criteria; read only, no new eligibility definition.
WITH e AS (
  SELECT * FROM json_to_recordset($1::json)
    AS q(stay_id integer, time_zero timestamp)
)
SELECT e.stay_id,
  EXTRACT(EPOCH FROM (i.intime-e.time_zero))/3600.0 AS icu_in_h,
  EXTRACT(EPOCH FROM (i.outtime-e.time_zero))/3600.0 AS icu_out_h
FROM e JOIN mimiciv_icu.icustays i USING(stay_id)
ORDER BY e.stay_id;
