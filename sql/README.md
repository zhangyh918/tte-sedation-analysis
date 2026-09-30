# SQL notes

The SQL files are intended for a local PostgreSQL installation of MIMIC-IV version 3.1.

The standard source schemas are `mimiciv_hosp` and `mimiciv_icu`. Project working tables are referenced through an `analysis` schema. The working schema contains locally generated tables for the cohort-exclusion chain, mechanical-ventilation intervals, RASS records and sedative-infusion intervals. These tables remain within the credentialed database environment and are not exported to this repository.

Several queries receive cohort identifiers and time zero through a parameterised JSON recordset (`$1::json`). The JSON input is generated locally and must not be committed.

Run order:

1. `01_candidate_cohorts.sql`
2. `02_baseline_covariates.sql` to `06_vasopressor_sofa.sql`
3. `07_primary_icu_scope.sql` to `09_primary_sedatives.sql`
4. `10_sensitivity_icu_scope.sql` to `12_sensitivity_sedatives.sql`

Change the `analysis` schema name in a local copy if a different working schema is used.

