# Target trial emulation of sedation strategies

This repository contains the analysis code supporting the manuscript **“Early Deep-to-Light Sedation versus a Sustained Light-or-Awake Strategy in Mechanically Ventilated Adults: A Target Trial Emulation.”**

The study used MIMIC-IV version 3.1. The repository contains code only. It does not contain MIMIC-IV records, patient-level analytic data, database credentials, or record-level derived datasets.

## Repository contents

- `sql/`: cohort, baseline, decision-node and analysis-scope queries for a local PostgreSQL installation of MIMIC-IV.
- `R/`: event ordering, ventilation-episode handling, strategy-history construction, inverse-probability-of-censoring weighting, outcome estimation and prespecified sensitivity analyses.
- `data/README.md`: data-access and data-sharing restrictions.
- `environment/install_packages.R`: R package installation commands.

## Data access

MIMIC-IV version 3.1 is available through PhysioNet to credentialed users who complete the required training and data use agreement:

- https://physionet.org/content/mimiciv/3.1/
- https://doi.org/10.13026/kpb9-mt58

MIMIC-IV data and patient-level derivatives must not be committed to this repository. See `DATA_ACCESS.md` and `data/README.md`.

## Software

- PostgreSQL with MIMIC-IV version 3.1 loaded in the standard `mimiciv_hosp` and `mimiciv_icu` schemas.
- R 4.x.
- R packages: `data.table` and `survival`.

Install the R packages with:

```r
source("environment/install_packages.R")
```

## Running the code

1. Clone or download this repository.
2. Start R from the repository root.
3. Set the repository root and load the analysis modules:

```r
Sys.setenv(TTE_REPO_ROOT = normalizePath("."))
source("R/load_modules.R")
```

4. Run the SQL queries against a credentialed local MIMIC-IV installation. The queries use a working schema called `analysis`; change that schema name locally if required. Parameterised queries using `$1::json` receive cohort identifiers and time zero from the local analysis session.
5. Construct the strategy histories and fit the models using the functions documented in the R source files. The primary analysis uses a 144-hour adherence window, a 48-hour transition window, modal RASS within each 12-hour node, and unstabilized inverse-probability-of-censoring weights.

The R modules are numbered in dependency order. `R/load_modules.R` loads them in that order.

## Reproducibility and confidentiality

All record-level inputs and outputs remain within the credentialed computing environment. Before running the code, users should confirm that their local MIMIC-IV schemas, item identifiers and derived tables match the manuscript definitions. Outputs should be reviewed for small cells before public release.

