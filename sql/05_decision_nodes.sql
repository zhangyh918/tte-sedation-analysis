WITH e AS (SELECT * FROM json_to_recordset($1::json) AS q(stay_id integer,subject_id integer,hadm_id integer,time_zero timestamp)), patient_nodes AS (
    SELECT
        e.stay_id,
        e.subject_id,
        gs.node_id,
        e.time_zero + gs.node_id * INTERVAL '12 hours' AS node_start,
        e.time_zero + (gs.node_id + 1) * INTERVAL '12 hours' AS node_end
    FROM e
    CROSS JOIN LATERAL generate_series(0, 13) AS gs(node_id)
),

node_rass AS (
    SELECT
        pn.stay_id,
        pn.node_id,
        MODE() WITHIN GROUP (ORDER BY rr.rass_score) AS modal_rass,
        AVG(rr.rass_score) AS mean_rass,
        MIN(rr.rass_score) AS min_rass,
        COUNT(rr.rass_score) AS rass_count
    FROM patient_nodes pn
    LEFT JOIN analysis._tmp_rass_raw rr
        ON rr.stay_id = pn.stay_id
       AND rr.rass_time >= pn.node_start
       AND rr.rass_time < pn.node_end
    GROUP BY pn.stay_id, pn.node_id
),

-- Assign each observed RASS node and its subsequent empty nodes to the same
-- carry-forward group. This implements true LOCF across any number of
-- consecutive empty nodes without a repeated correlated lookup.
node_rass_groups AS (
    SELECT
        nr.*,
        COUNT(*) FILTER (WHERE nr.rass_count > 0) OVER (
            PARTITION BY nr.stay_id
            ORDER BY nr.node_id
            ROWS UNBOUNDED PRECEDING
        ) AS locf_group
    FROM node_rass nr
),

node_rass_locf AS (
    SELECT
        rg.stay_id,
        rg.node_id,
        rg.modal_rass,
        rg.mean_rass,
        rg.min_rass,
        rg.rass_count,
        MAX(rg.modal_rass) OVER (
            PARTITION BY rg.stay_id, rg.locf_group
        ) AS modal_rass_locf,
        MAX(rg.mean_rass) OVER (
            PARTITION BY rg.stay_id, rg.locf_group
        ) AS mean_rass_locf,
        MAX(rg.min_rass) OVER (
            PARTITION BY rg.stay_id, rg.locf_group
        ) AS min_rass_locf,
        (rg.rass_count = 0) AS rass_imputed
    FROM node_rass_groups rg
),

qualifying_sedatives AS (
    SELECT
        sr.stay_id,
        sr.infusion_start,
        sr.infusion_end,
        sr.drug_class,
        sr.rate,
        sr.rateuom
    FROM analysis._tmp_sedative_raw sr
    WHERE sr.infusion_end > sr.infusion_start
      AND (
          (
              sr.drug_class = 'propofol'
              AND sr.rateuom = 'mcg/kg/min'
              AND sr.rate >= 5
          )
          OR (
              sr.drug_class = 'midazolam'
              AND sr.rateuom = 'mg/hour'
              AND sr.rate >= 0.5
          )
          OR (
              sr.drug_class = 'dexmedetomidine'
              AND sr.rateuom = 'mcg/kg/hour'
              AND sr.rate >= 0.1
          )
      )
),

-- Clip every qualifying infusion to its overlapping 12-hour node.
node_intervals AS (
    SELECT
        pn.stay_id,
        pn.node_id,
        GREATEST(qs.infusion_start, pn.node_start) AS interval_start,
        LEAST(qs.infusion_end, pn.node_end) AS interval_end,
        qs.drug_class,
        qs.rate
    FROM patient_nodes pn
    JOIN qualifying_sedatives qs
      ON qs.stay_id = pn.stay_id
     AND qs.infusion_start < pn.node_end
     AND qs.infusion_end > pn.node_start
),

-- Running maximum is required for correct gaps-and-islands merging when an
-- early long interval contains several shorter intervals.
ordered_intervals AS (
    SELECT
        ni.*,
        MAX(ni.interval_end) OVER (
            PARTITION BY ni.stay_id, ni.node_id
            ORDER BY ni.interval_start, ni.interval_end
            ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING
        ) AS prior_max_end
    FROM node_intervals ni
    WHERE ni.interval_end > ni.interval_start
),

island_markers AS (
    SELECT
        oi.*,
        CASE
            WHEN oi.prior_max_end IS NULL
              OR oi.interval_start > oi.prior_max_end
            THEN 1 ELSE 0
        END AS starts_new_island
    FROM ordered_intervals oi
),

island_ids AS (
    SELECT
        im.*,
        SUM(im.starts_new_island) OVER (
            PARTITION BY im.stay_id, im.node_id
            ORDER BY im.interval_start, im.interval_end
            ROWS UNBOUNDED PRECEDING
        ) AS island_id
    FROM island_markers im
),

merged_intervals AS (
    SELECT
        stay_id,
        node_id,
        island_id,
        MIN(interval_start) AS merged_start,
        MAX(interval_end) AS merged_end
    FROM island_ids
    GROUP BY stay_id, node_id, island_id
),

node_coverage AS (
    SELECT
        stay_id,
        node_id,
        SUM(EXTRACT(EPOCH FROM (merged_end - merged_start)) / 3600.0)
            AS qualifying_sed_hours
    FROM merged_intervals
    GROUP BY stay_id, node_id
),

node_drug_details AS (
    SELECT
        stay_id,
        node_id,
        STRING_AGG(DISTINCT drug_class, ', ' ORDER BY drug_class)
            AS qualifying_drugs,
        MAX(rate) FILTER (WHERE drug_class = 'propofol')
            AS max_propofol_rate,
        MAX(rate) FILTER (WHERE drug_class = 'midazolam')
            AS max_midazolam_rate,
        MAX(rate) FILTER (WHERE drug_class = 'dexmedetomidine')
            AS max_dexmedetomidine_rate
    FROM node_intervals
    GROUP BY stay_id, node_id
)

SELECT
    pn.stay_id,
    pn.subject_id,
    pn.node_id,
    pn.node_start,
    pn.node_end,
    rl.modal_rass,
    rl.mean_rass,
    rl.min_rass,
    rl.rass_count,
    rl.modal_rass_locf,
    rl.mean_rass_locf,
    rl.min_rass_locf,
    rl.rass_imputed,
    COALESCE(nc.qualifying_sed_hours, 0.0) AS qualifying_sed_hours,
    (COALESCE(nc.qualifying_sed_hours, 0.0) >= 5.999999)
        AS qualifying_sed_active,
    nd.qualifying_drugs,
    nd.max_propofol_rate,
    nd.max_midazolam_rate,
    nd.max_dexmedetomidine_rate,
    CASE
        WHEN rl.modal_rass_locf IS NULL THEN 'missing'
        WHEN rl.modal_rass_locf >= -2
         AND COALESCE(nc.qualifying_sed_hours, 0.0) >= 5.999999
            THEN 'light_with_qualifying_infusion'
        WHEN rl.modal_rass_locf >= -2
         AND COALESCE(nc.qualifying_sed_hours, 0.0) < 5.999999
            THEN 'light_without_qualifying_infusion'
        WHEN rl.modal_rass_locf <= -3
         AND COALESCE(nc.qualifying_sed_hours, 0.0) >= 5.999999
            THEN 'deep_with_qualifying_infusion'
        WHEN rl.modal_rass_locf <= -3
         AND COALESCE(nc.qualifying_sed_hours, 0.0) < 5.999999
            THEN 'deep_without_qualifying_infusion'
        ELSE 'unclassified'
    END AS joint_class_primary
FROM patient_nodes pn
JOIN node_rass_locf rl
  ON rl.stay_id = pn.stay_id
 AND rl.node_id = pn.node_id
LEFT JOIN node_coverage nc
  ON nc.stay_id = pn.stay_id
 AND nc.node_id = pn.node_id
LEFT JOIN node_drug_details nd
  ON nd.stay_id = pn.stay_id
 AND nd.node_id = pn.node_id
ORDER BY pn.stay_id, pn.node_id;
