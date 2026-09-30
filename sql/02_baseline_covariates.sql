WITH e AS (SELECT * FROM json_to_recordset($1::json) AS q(stay_id integer,subject_id integer,hadm_id integer,time_zero timestamp)), vitals_raw AS (SELECT e.stay_id,ce.itemid,ce.valuenum,ROW_NUMBER() OVER(PARTITION BY e.stay_id,ce.itemid ORDER BY ce.charttime DESC) rn FROM e JOIN mimiciv_icu.chartevents ce ON ce.stay_id=e.stay_id AND ce.charttime<=e.time_zero AND ce.charttime>=e.time_zero-INTERVAL '24 hours' AND ce.valuenum IS NOT NULL WHERE ce.itemid IN (220045,220052,220181,225312,220277,223835,227009,220339,224700)), labs_raw AS (SELECT e.stay_id,le.itemid,le.valuenum,ROW_NUMBER() OVER(PARTITION BY e.stay_id,le.itemid ORDER BY le.charttime DESC) rn FROM e JOIN mimiciv_hosp.labevents le ON le.subject_id=e.subject_id AND le.hadm_id=e.hadm_id AND le.charttime<=e.time_zero AND le.charttime>=e.time_zero-INTERVAL '48 hours' AND le.valuenum>0 WHERE le.itemid IN (50813,50912,50885,51265,50821)), vitals_pivot AS (
    SELECT 
        stay_id,
        MAX(CASE WHEN itemid = 220045 THEN valuenum END) AS heart_rate,
        MAX(CASE WHEN itemid IN (220052, 220181, 225312) THEN valuenum END) AS map_value,
        MAX(CASE WHEN itemid = 220277 THEN valuenum END) AS spo2_value,
        MAX(CASE WHEN itemid IN (223835, 227009) THEN valuenum END) AS fio2_value,
        MAX(CASE WHEN itemid IN (220339, 224700) THEN valuenum END) AS peep_value
    FROM vitals_raw
    WHERE rn = 1
    GROUP BY stay_id
),

labs_pivot AS (
    SELECT 
        stay_id,
        MAX(CASE WHEN itemid = 50813 THEN valuenum END) AS lactate,
        MAX(CASE WHEN itemid = 50912 THEN valuenum END) AS creatinine,
        MAX(CASE WHEN itemid = 50885 THEN valuenum END) AS bilirubin,
        MAX(CASE WHEN itemid = 51265 THEN valuenum END) AS platelets,
        MAX(CASE WHEN itemid = 50821 THEN valuenum END) AS pao2
    FROM labs_raw
    WHERE rn = 1
    GROUP BY stay_id
),

vaso_at_t0 AS (
    SELECT 
        e.stay_id,
        TRUE AS vasopressor_use
    FROM (SELECT q.*,c.age,c.gender,c.race,c.admission_type,c.hospital_expire_flag,c.hospital_deathtime,c.icu_intime,c.icu_outtime,c.total_mv_hours,c.first_mv_start FROM e q JOIN analysis._tmp_exclusion_chain c USING(stay_id)) e
    INNER JOIN mimiciv_icu.inputevents ie ON e.stay_id = ie.stay_id
    WHERE ie.itemid IN (221906, 221289, 222315, 221749, 221662)
      AND ie.starttime <= e.time_zero AND ie.endtime >= e.time_zero
      AND ie.rate > 0
    GROUP BY e.stay_id
),

initial_sed AS (
    SELECT DISTINCT ON (e.stay_id) 
        e.stay_id, s.drug_class AS initial_sedative
    FROM (SELECT q.*,c.age,c.gender,c.race,c.admission_type,c.hospital_expire_flag,c.hospital_deathtime,c.icu_intime,c.icu_outtime,c.total_mv_hours,c.first_mv_start FROM e q JOIN analysis._tmp_exclusion_chain c USING(stay_id)) e
    INNER JOIN analysis._tmp_sedative_raw s ON e.stay_id = s.stay_id
    WHERE s.infusion_start <= e.time_zero
    ORDER BY e.stay_id, s.infusion_start DESC
),

comorbid AS (
    SELECT 
        hadm_id,
        MAX(CASE WHEN icd_code LIKE '250%' OR icd_code LIKE 'E11%' OR icd_code LIKE 'E10%' THEN 1 ELSE 0 END) AS diabetes,
        MAX(CASE WHEN icd_code LIKE '585%' OR icd_code LIKE 'N18%' THEN 1 ELSE 0 END) AS ckd,
        MAX(CASE WHEN icd_code LIKE '571%' OR icd_code LIKE 'K74%' OR icd_code LIKE 'K70%' THEN 1 ELSE 0 END) AS liver_disease,
        MAX(CASE WHEN icd_code LIKE '428%' OR icd_code LIKE 'I50%' THEN 1 ELSE 0 END) AS chf,
        MAX(CASE WHEN icd_code LIKE '490%' OR icd_code LIKE '491%' OR icd_code LIKE '492%'
            OR icd_code LIKE 'J44%' OR icd_code LIKE 'J43%' THEN 1 ELSE 0 END) AS copd
    FROM mimiciv_hosp.diagnoses_icd
    GROUP BY hadm_id
)

SELECT 
    e.stay_id, e.subject_id, e.hadm_id, e.time_zero,
    e.age, e.gender, e.race, e.admission_type,
    e.hospital_expire_flag, e.hospital_deathtime,
    e.icu_intime, e.icu_outtime, e.total_mv_hours, e.first_mv_start,
    vp.heart_rate, vp.map_value, vp.spo2_value,
    vp.fio2_value, vp.peep_value,
    lp.lactate, lp.creatinine, lp.bilirubin, lp.platelets, lp.pao2,
    CASE 
        WHEN lp.pao2 IS NOT NULL AND vp.fio2_value IS NOT NULL AND vp.fio2_value > 0 
        THEN ROUND((lp.pao2 / (vp.fio2_value / 100.0))::numeric, 1) 
    END AS pf_ratio,
    COALESCE(va.vasopressor_use, FALSE) AS vasopressor_at_t0,
    ised.initial_sedative,
    COALESCE(cm.diabetes, 0) AS diabetes,
    COALESCE(cm.ckd, 0) AS ckd,
    COALESCE(cm.liver_disease, 0) AS liver_disease,
    COALESCE(cm.chf, 0) AS chf,
    COALESCE(cm.copd, 0) AS copd
FROM (SELECT q.*,c.age,c.gender,c.race,c.admission_type,c.hospital_expire_flag,c.hospital_deathtime,c.icu_intime,c.icu_outtime,c.total_mv_hours,c.first_mv_start FROM e q JOIN analysis._tmp_exclusion_chain c USING(stay_id)) e
LEFT JOIN vitals_pivot vp ON e.stay_id = vp.stay_id
LEFT JOIN labs_pivot lp ON e.stay_id = lp.stay_id
LEFT JOIN vaso_at_t0 va ON e.stay_id = va.stay_id
LEFT JOIN initial_sed ised ON e.stay_id = ised.stay_id
LEFT JOIN comorbid cm ON e.hadm_id = cm.hadm_id;
