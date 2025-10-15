# 01_setup_packages.R
options(warn = 1)
options(readr.show_col_types = FALSE)

pkgs <- c(
  "tidyverse","data.table","janitor","lubridate","stringr","readr",
  "ggplot2","scales","patchwork",
  "pdftools",
  "heemod","hesim",
  "fitdistrplus","truncnorm",
  "xgboost","lightgbm","glmnet","mgcv","splines",
  "conflicted"   # <-- add this
)

to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

# Prefer dplyr when functions are masked
conflicted::conflict_prefer("select",  "dplyr", quiet = TRUE)
conflicted::conflict_prefer("filter",  "dplyr", quiet = TRUE)
conflicted::conflict_prefer("summarise","dplyr", quiet = TRUE)
conflicted::conflict_prefer("mutate",  "dplyr", quiet = TRUE)
conflicted::conflict_prefer("arrange", "dplyr", quiet = TRUE)
conflicted::conflict_prefer("left_join","dplyr", quiet = TRUE)