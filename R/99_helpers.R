# 99_helpers.R
# -------------------------------

# ---- Project root: set ONCE here ----
DIR.PROJ <- "/Users/eshtiaghi/Desktop/Khashi/Chapman/Research/My works/Cost Analysis/Analysis"
DIR.CLEAN <- file.path(DIR.PROJ, "data_clean")
DIR.OUT.FIG <- file.path(DIR.PROJ, "outputs", "figures")
DIR.OUT.TAB <- file.path(DIR.PROJ, "outputs", "tables")
dir.create(DIR.CLEAN, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR.OUT.FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR.OUT.TAB, recursive = TRUE, showWarnings = FALSE)

as_prob <- function(rate_per_100k) rate_per_100k / 1e5

norm_age <- function(x) {
  x <- tolower(x)
  x <- stringr::str_replace_all(x, "–", "-")
  x <- stringr::str_replace_all(x, "\\s+", "")
  x <- stringr::str_replace_all(x, "years|year|\\+", "")
  x
}

# piecewise-linear waning curve builder
make_waning_fn <- function(time_points_months, ve_points) {
  stopifnot(length(time_points_months) == length(ve_points))
  function(t_months) {
    approx(x = time_points_months, y = ve_points, xout = t_months,
           rule = 2, method = "linear")$y
  }
}

# assert columns exist
ensure_cols <- function(df, cols, df_name = "data") {
  missing <- setdiff(cols, names(df))
  if (length(missing)) {
    stop(sprintf("Missing columns in %s: %s", df_name, paste(missing, collapse=", ")))
  }
  invisible(TRUE)
}

# safe pull with default
pull_or <- function(df, expr, default) {
  v <- tryCatch(dplyr::pull(df, {{ expr }}), error = function(e) NA_real_)
  if (length(v) == 0 || all(is.na(v))) default else as.numeric(v[1])
}

# probability guards
cap01 <- function(x, eps=1e-6) pmax(pmin(as.numeric(x), 1 - eps), 0)
renorm_H <- function(p_h, p_o, p_d, maxsum=0.999) {
  s <- p_h + p_o + p_d
  if (s <= maxsum) return(c(p_h, p_o, p_d))
  scale <- maxsum / s
  c(p_h * scale, p_o * scale, p_d * scale)
}