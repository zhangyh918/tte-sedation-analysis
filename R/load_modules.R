repo_root <- Sys.getenv("TTE_REPO_ROOT", unset = normalizePath(".", mustWork = TRUE))

modules <- c(
  "01_event_order_core.R",
  "02_outcome_time_core.R",
  "03_event_history_core.R",
  "04_mv_pause_resume_core.R",
  "05_complete_history_core.R",
  "06_ipcw_estimation_core.R",
  "07_sofa_subgroup_core.R",
  "08_sensitivity_core.R",
  "09_coverage_sensitivity_core.R"
)

for (module in modules) {
  source(file.path(repo_root, "R", module))
}

