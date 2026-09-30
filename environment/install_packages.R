required_packages <- c("data.table", "survival")
missing_packages <- setdiff(required_packages, rownames(installed.packages()))

if (length(missing_packages)) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

