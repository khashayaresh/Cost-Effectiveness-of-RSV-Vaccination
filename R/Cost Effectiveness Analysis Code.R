# -------------------------------
# 02_ingest_rsvnet.R
# -------------------------------
source(file.path(DIR.PROJ, "R", "01_setup_packages.R"))
source(file.path(DIR.PROJ, "R", "99_helpers.R"))

# ensure output dirs exist
dir.create(DIR.OUT.FIG, showWarnings = FALSE, recursive = TRUE)

csv_path <- file.path(
  DIR.PROJ,
  "Weekly_Rates_of_Laboratory-Confirmed_RSV_Hospitalizations_from_the_RSV-NET_Surveillance_System_20251002.csv"
)

rsvnet_raw <- readr::read_csv(csv_path, show_col_types = FALSE) %>% janitor::clean_names()
ensure_cols(rsvnet_raw, c("age_category","season","cumulative_rate"), "rsvnet_raw")

adult_labels <- c("50–64 years","50-64 years","65–74 years","65-74 years",
                  "75–84 years","75-84 years","≥85 years","85+ years",
                  "≥18 years (adults)","18+ years (adults)")

rsvnet_adult <- rsvnet_raw %>% dplyr::filter(age_category %in% adult_labels)

# Make 'season' an ordered factor (nice plotting)
rsvnet_season_age <- rsvnet_adult %>%
  dplyr::mutate(
    age_norm = norm_age(age_category),
    season = factor(season, levels = sort(unique(season)))
  ) %>%
  dplyr::group_by(season, age_category, age_norm) %>%
  dplyr::summarise(
    season_cum_rate_per100k = max(cumulative_rate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(season, age_category)

readr::write_csv(rsvnet_season_age, file.path(DIR.PROJ, "rsvnet_season_age_unadjusted.csv"))

# sanity plot
p <- ggplot(rsvnet_season_age,
            aes(x = season, y = season_cum_rate_per100k,
                color = age_category, group = age_category)) +
  geom_line() + geom_point() + coord_flip() +
  labs(title="RSV-NET seasonal cumulative hospitalization rate (unadjusted)",
       y="per 100,000", x="Season") + theme_minimal()

ggsave(file.path(DIR.OUT.FIG, "rsvnet_unadjusted_by_season_age.png"),
       p, width=9, height=6, dpi=300)

# -------------------------------
# 03_create_curated_inputs.R
# -------------------------------
source(file.path(DIR.PROJ, "R", "01_setup_packages.R"))
source(file.path(DIR.PROJ, "R", "99_helpers.R"))

# 1) Havers — keep labels aligned with norm_age()
havers_structured <- tibble::tribble(
  ~age_band,     ~adjusted_rate_min_per100k, ~adjusted_rate_max_per100k, ~icu_pct, ~in_hosp_death_pct, ~undertesting_multiplier,
  "50-64 years", NA_real_,                   NA_real_,                   0.191,    0.043,              1.5,
  "65-74 years", NA_real_,                   NA_real_,                   0.191,    0.043,              1.5,
  "75-84 years", 244.7,                      411.4,                      0.191,    0.058,              1.5,
  "85+ years",   NA_real_,                   NA_real_,                   0.191,    0.058,              1.5
) %>% janitor::clean_names()

readr::write_csv(havers_structured, file.path(DIR.CLEAN, "havers_structured.csv"))

# 2) ACIP (incidence & costs; VE goes to its own file)
acip_structured <- tibble::tribble(
  ~age_band, ~chronic_status,         ~incidence_range_per100k_low, ~incidence_range_per100k_high, ~hosp_cost_mean, ~hosp_cost_low, ~hosp_cost_high, ~ve_months, ~ve_values,
  ">=60",    "with_chronic",          198,                           527,                           22000,           NA,             NA,              NA,         NA,
  ">=60",    "without_chronic",       32,                            121,                           11876,           8407,           47512,           NA,         NA,
  "60-64",   "general_outpatient",    1833,                          1833,                          NA,              NA,             NA,              NA,         NA,
  ">=65",    "general_outpatient",    2478,                          2478,                          NA,              NA,             NA,              NA,         NA,
  "50-59",   "COPD_outpatient",       2925,                          2925,                          NA,              NA,             NA,              NA,         NA,
  "50-59",   "with_chronic_hosp",     106,                           106,                           20330,           NA,             NA,              NA,         NA,
  "50-59",   "COPD_hosp",             169,                           169,                           35308,           NA,             NA,              NA,         NA
) %>% janitor::clean_names()

readr::write_csv(acip_structured, file.path(DIR.CLEAN, "acip_structured.csv"))

# 3) Landi 28-day O→I
landi_structured <- tibble::tribble(
  ~age_band,   ~risk_group, ~outpatient_to_hosp_28d_pct,
  "adult_all", "overall",   6.2,
  "adult_all", "overall",   6.0,
  "adult_all", "overall",   4.5,
  "adult_all", "high_risk", 7.6,
  "adult_all", "high_risk", 8.5,
  "adult_all", "high_risk", 6.5
) %>% janitor::clean_names()

readr::write_csv(landi_structured, file.path(DIR.CLEAN, "landi_structured.csv"))

# 4) EQ-5D age norms (keep in sync with norm_age())
eq5d_age_util_structured <- tibble::tribble(
  ~age_norm, ~utility_base,
  "55-64",   0.815,
  "65-74",   0.824,
  "75-84",   0.811,
  "85plus",  0.800     # placeholder until you add exact value
)
readr::write_csv(eq5d_age_util_structured, file.path(DIR.CLEAN, "eq5d_age_util_structured.csv"))

# 5) Costs & disutilities
costs_utilities <- tibble::tribble(
  ~item,                         ~value,
  "outpatient_cost",              200,
  "hosp_cost_mean",               22000,
  "hosp_cost_low",                11876,
  "hosp_cost_high",               47512,
  "disutility_outpt",             0.02,
  "disutility_hosp_nonicu",       0.10,
  "disutility_hosp_icu",          0.20,
  "vaccine_price",                270,
  "admin_cost",                   25
)
readr::write_csv(costs_utilities, file.path(DIR.CLEAN, "costs_utilities.csv"))

message("Curated input CSVs written to: ", DIR.CLEAN)

# -------------------------------
# 04_build_inputs.R
# -------------------------------
source(file.path(DIR.PROJ, "R", "01_setup_packages.R"))
source(file.path(DIR.PROJ, "R", "99_helpers.R"))

rsvnet <- readr::read_csv(file.path(DIR.PROJ, "rsvnet_season_age_unadjusted.csv"), show_col_types = FALSE) %>%
  dplyr::mutate(age_norm = norm_age(age_category))
havers_struct <- readr::read_csv(file.path(DIR.CLEAN, "havers_structured.csv"), show_col_types = FALSE) %>%
  dplyr::mutate(age_norm = norm_age(age_band))

ensure_cols(rsvnet,       c("season","age_category","age_norm","season_cum_rate_per100k"), "rsvnet")
ensure_cols(havers_struct,c("age_norm","icu_pct","in_hosp_death_pct","undertesting_multiplier"), "havers_struct")

df <- rsvnet %>%
  dplyr::left_join(
    dplyr::select(havers_struct, age_norm, icu_pct, in_hosp_death_pct, undertesting_multiplier),
    by = "age_norm"
  ) %>%
  dplyr::mutate(
    undertest_mult = dplyr::if_else(is.na(undertesting_multiplier), 1.0, undertesting_multiplier),
    pcr_mult       = 1.5,
    adj_per100k    = season_cum_rate_per100k * undertest_mult * pcr_mult
  )

readr::write_csv(df, file.path(DIR.CLEAN, "adjusted_incidence_age.csv"))

base_incidence_age <- df %>%
  dplyr::group_by(age_norm, age_category) %>%
  dplyr::summarise(
    adj_rate_per100k_median = median(adj_per100k, na.rm = TRUE),
    icu_pct = dplyr::first(icu_pct),
    in_hosp_death_pct = dplyr::first(in_hosp_death_pct),
    .groups = "drop"
  ) %>%
  dplyr::mutate(annual_hosp_prob = adj_rate_per100k_median / 1e5)

readr::write_csv(base_incidence_age, file.path(DIR.CLEAN, "base_incidence_age.csv"))
message("Adjusted incidence + base_incidence_age written to: ", DIR.CLEAN)

# -------------------------------
# 05_ml_risk_model.R
# -------------------------------
source(file.path(DIR.PROJ, "R", "01_setup_packages.R"))
source(file.path(DIR.PROJ, "R", "99_helpers.R"))

base_age    <- readr::read_csv(file.path(DIR.CLEAN, "base_incidence_age.csv"), show_col_types = FALSE)
acip_struct <- readr::read_csv(file.path(DIR.CLEAN, "acip_structured.csv"), show_col_types = FALSE)

ensure_cols(base_age, c("age_norm","annual_hosp_prob"), "base_incidence_age")

set.seed(1)
N <- 100000
age_levels <- unique(base_age$age_norm)
age_probs  <- rep(1 / length(age_levels), length(age_levels))

pop <- tibble::tibble(
  age_norm = sample(age_levels, N, replace = TRUE, prob = age_probs),
  COPD   = rbinom(N, 1, 0.15),
  CHF    = rbinom(N, 1, 0.12),
  ASTHMA = rbinom(N, 1, 0.12),
  DM     = rbinom(N, 1, 0.25),
  CKD    = rbinom(N, 1, 0.10),
  IMMUNO = rbinom(N, 1, 0.06),
  OBESE  = rbinom(N, 1, 0.35)
) %>%
  dplyr::mutate(CHRONIC = as.integer(COPD | CHF | ASTHMA | DM | CKD | IMMUNO | OBESE)) %>%
  dplyr::left_join(dplyr::select(base_age, age_norm, annual_hosp_prob), by = "age_norm")

design  <- model.matrix(~ 0 + age_norm + CHRONIC + COPD + CHF + ASTHMA + DM + CKD + IMMUNO + OBESE, data = pop)
beta    <- rep(0, ncol(design)); names(beta) <- colnames(design)
sigmoid <- function(z) 1/(1+exp(-z))

objective <- function(b) {
  lin <- as.vector(design %*% b); p <- sigmoid(lin)
  # 1) match age-marginal risks
  pen_age <- pop %>% dplyr::mutate(p=p) %>%
    dplyr::group_by(age_norm) %>%
    dplyr::summarise(diff2 = (mean(p) - mean(annual_hosp_prob))^2, .groups = "drop") %>%
    dplyr::summarise(sum(diff2)) %>% dplyr::pull()
  # 2) chronic vs non-chronic (soft target from ACIP)
  pen_chronic <- 0
  if (all(c("incidence_range_per100k_low","incidence_range_per100k_high","chronic_status","age_band") %in% names(acip_struct))) {
    targ <- acip_struct %>%
      dplyr::mutate(age_norm = norm_age(age_band)) %>%
      dplyr::group_by(chronic_status) %>%
      dplyr::summarise(lo = mean(incidence_range_per100k_low, na.rm = TRUE)/1e5,
                       hi = mean(incidence_range_per100k_high, na.rm = TRUE)/1e5, .groups="drop")
    p_chronic <- mean(p[pop$CHRONIC==1]); p_non <- mean(p[pop$CHRONIC==0])
    lo_ch <- targ$lo[grepl("with_chronic|chronic", tolower(targ$chronic_status))]
    hi_ch <- targ$hi[grepl("with_chronic|chronic", tolower(targ$chronic_status))]
    pen_chronic <- (max(0, lo_ch - p_chronic) + max(0, p_chronic - hi_ch))^2 +
      0.1*max(0, (p_non - p_chronic + 0.002))^2
  }
  # 3) optional ≥75 share target ~45.6% (Havers)
  target_share_75p <- 0.456
  idx75 <- pop$age_norm %in% c("75-84years","85years","85+","75yearsplus","75-84","85plus")
  pen_age_share <- if (any(idx75)) {
    w_tot <- sum(p); w_75 <- sum(p[idx75]); (w_75/w_tot - target_share_75p)^2
  } else 0
  pen_age + pen_chronic + pen_age_share
}

opt <- optim(beta, objective, method="BFGS", control=list(maxit=1500))
beta_hat <- opt$par
pop$annual_hosp_prob_ml <- sigmoid(as.vector(design %*% beta_hat))

risk_quintiles <- pop %>%
  dplyr::group_by(age_norm) %>%
  dplyr::mutate(quint = dplyr::ntile(annual_hosp_prob_ml, 5)) %>%
  dplyr::group_by(age_norm, quint) %>%
  dplyr::summarise(mean_prob = mean(annual_hosp_prob_ml), n = dplyr::n(), .groups = "drop")

readr::write_csv(risk_quintiles, file.path(DIR.CLEAN, "risk_quintiles_by_age.csv"))
message("Risk quintiles written to: ", DIR.CLEAN)

# --- VE waning points (Arexvy-like) ---
ve_points <- tibble::tibble(
  outcome = c(rep("hosp",4), rep("outpt",4)),
  months  = c(6, 12, 24, 36,   6, 12, 24, 36),
  ve      = c(0.82, 0.75, 0.62, 0.50, 0.50, 0.45, 0.35, 0.30)
)
readr::write_csv(ve_points, file.path(DIR.CLEAN, "ve_waning_points.csv"))
message("Wrote: ", file.path(DIR.CLEAN, "ve_waning_points.csv"))

# -------------------------------
# 06_markov_cea.R  (revised)
# -------------------------------
source(file.path(DIR.PROJ, "R", "01_setup_packages.R"))
source(file.path(DIR.PROJ, "R", "99_helpers.R"))

# --- Load inputs ---
inc_by_age  <- readr::read_csv(file.path(DIR.CLEAN, "base_incidence_age.csv"), show_col_types = FALSE)
eq_util     <- readr::read_csv(file.path(DIR.CLEAN, "eq5d_age_util_structured.csv"), show_col_types = FALSE)
costs_utils <- readr::read_csv(file.path(DIR.CLEAN, "costs_utilities.csv"), show_col_types = FALSE)

ensure_cols(inc_by_age, c("age_norm","annual_hosp_prob","icu_pct","in_hosp_death_pct"), "base_incidence_age")
ensure_cols(eq_util,    c("age_norm","utility_base"),                                   "eq5d_age_util_structured")
ensure_cols(costs_utils,c("item","value"),                                             "costs_utilities")

# --- Normalize + choose cohort ---
inc_by_age <- inc_by_age %>% dplyr::mutate(age_norm = norm_age(age_norm))
eq_util    <- eq_util    %>% dplyr::mutate(age_norm = norm_age(age_norm))

age_sel_raw <- "75-84"
age_sel     <- norm_age(age_sel_raw)
if (!(age_sel %in% inc_by_age$age_norm)) {
  message("age_sel not found; using first available.")
  age_sel <- inc_by_age$age_norm[1]
}

# --- Scalars (gentle floor only if missing) ---
p_hosp_base_raw <- pull_or(inc_by_age %>% dplyr::filter(age_norm == age_sel), annual_hosp_prob, NA_real_)
p_hosp_base <- if (!is.na(p_hosp_base_raw)) p_hosp_base_raw else dplyr::case_when(
  grepl("^85", age_sel) ~ 0.0080,
  grepl("^75", age_sel) ~ 0.0035,
  grepl("^65", age_sel) ~ 0.0025,
  TRUE                  ~ 0.0015
)

icu_pct       <- pull_or(inc_by_age %>% dplyr::filter(age_norm == age_sel), icu_pct,        0.191)
in_hosp_death <- pull_or(inc_by_age %>% dplyr::filter(age_norm == age_sel), in_hosp_death_pct, 0.058)
util_base     <- pull_or(eq_util    %>% dplyr::filter(age_norm == age_sel), utility_base,   0.811)

# Costs
outpatient_cost <- pull_or(costs_utils %>% dplyr::filter(item=="outpatient_cost"), value, 150)
hosp_cost_mean  <- pull_or(costs_utils %>% dplyr::filter(item=="hosp_cost_mean"),  value, 22000)
vax_price       <- pull_or(costs_utils %>% dplyr::filter(item=="vaccine_price"),   value, 270)
admin_cost      <- pull_or(costs_utils %>% dplyr::filter(item=="admin_cost"),      value, 25)
vax_total_cost  <- vax_price + admin_cost

# Disutilities
du_O        <- pull_or(costs_utils %>% dplyr::filter(item=="disutility_outpt"),       value, 0.012)
du_I_nonicu <- pull_or(costs_utils %>% dplyr::filter(item=="disutility_hosp_nonicu"), value, 0.13)
du_I_icu    <- pull_or(costs_utils %>% dplyr::filter(item=="disutility_hosp_icu"),    value, 0.25)

# Outpatient & background mortality (proxy)
p_outpt <- dplyr::case_when(
  grepl("^5", age_sel)          ~ 1833/1e5,
  grepl("^6", age_sel)          ~ 1833/1e5,
  grepl("^7|85|\\+|≥", age_sel) ~ 2478/1e5,
  TRUE                          ~ 2000/1e5
)
p_bg_mort <- dplyr::case_when(
  grepl("^5", age_sel)          ~ 0.005,
  grepl("^6", age_sel)          ~ 0.007,
  grepl("^7", age_sel)          ~ 0.020,
  grepl("85|\\+|≥", age_sel)    ~ 0.050,
  TRUE                          ~ 0.010
)

# VE waning (dual)
ve_file <- file.path(DIR.CLEAN, "ve_waning_points.csv")
if (file.exists(ve_file)) {
  ve_pts <- readr::read_csv(ve_file, show_col_types = FALSE) %>% dplyr::arrange(outcome, months)
  waning_hosp  <- make_waning_fn(ve_pts %>% dplyr::filter(outcome=="hosp")  %>% dplyr::pull(months),
                                 ve_pts %>% dplyr::filter(outcome=="hosp")  %>% dplyr::pull(ve))
  waning_outpt <- make_waning_fn(ve_pts %>% dplyr::filter(outcome=="outpt") %>% dplyr::pull(months),
                                 ve_pts %>% dplyr::filter(outcome=="outpt") %>% dplyr::pull(ve))
} else {
  waning_hosp  <- function(m) rep(0.82, length(m))
  waning_outpt <- function(m) rep(0.50, length(m))
}

disc_rate <- 0.03
cycles    <- 15
ve_months <- seq(6, cycles*12, by=12)
ve_hosp   <- waning_hosp(ve_months)
ve_outpt <- pmax(0, pmin(1, 0.6 * ve_hosp))

# NO-DOUBLE-COUNTING (conservative)
share_hosp_via_outpt <- if (grepl("^75|85|\\+|≥", age_sel)) 0.25 else 0.20
share_hosp_via_outpt <- max(0, min(share_hosp_via_outpt, 0.90))

cap01   <- function(x, eps=1e-6) pmax(pmin(as.numeric(x), 1 - eps), 0)
renorm_H<- function(p_i, p_o, p_d, maxsum=0.999) {
  s <- p_i + p_o + p_d
  if (s <= maxsum) return(c(p_i, p_o, p_d))
  scale <- maxsum / s; c(p_i*scale, p_o*scale, p_d*scale)
}

make_tp_consistent <- function(p_hosp, p_outpt, p_bg_mort, p_case_fatal,
                               ve_h = 0, ve_o = 0, age_sel = "75-84",
                               share_via_O = 0.30) {
  p_hosp_eff  <- cap01(p_hosp  * (1 - ve_h))
  p_outpt_eff <- cap01(p_outpt * (1 - ve_o))
  p_bg_mort   <- cap01(p_bg_mort)
  p_case_fatal<- cap01(p_case_fatal)
  
  # choose O->I so H->I + (H->O)*O->I = p_hosp_eff
  p_OI <- if (p_outpt_eff > 0) (p_hosp_eff * share_via_O) / p_outpt_eff else 0
  p_OI <- cap01(min(p_OI, 0.25))  # annual cap
  
  p_HI_direct <- p_hosp_eff - p_outpt_eff * p_OI
  p_HI_direct <- cap01(p_HI_direct)
  
  H_vals <- renorm_H(p_HI_direct, p_outpt_eff, p_bg_mort, maxsum = 0.999)
  p_HI_direct <- H_vals[1]; p_outpt_eff <- H_vals[2]; p_bg_mort <- H_vals[3]
  
  m <- matrix(0, nrow=4, ncol=4, dimnames=list(c("H","O","I","D"), c("H","O","I","D")))
  m["H","I"] <- p_HI_direct
  m["H","O"] <- p_outpt_eff
  m["H","D"] <- p_bg_mort
  m["H","H"] <- 1 - (m["H","I"] + m["H","O"] + m["H","D"])
  m["O","I"] <- p_OI
  m["O","D"] <- cap01(0.001)
  m["O","H"] <- 1 - (m["O","I"] + m["O","D"])
  m["I","D"] <- p_case_fatal
  m["I","H"] <- 1 - m["I","D"]
  m["D","D"] <- 1
  m
}

tp_no_vax <- lapply(seq_len(cycles), function(t) {
  make_tp_consistent(p_hosp_base, p_outpt, p_bg_mort, in_hosp_death,
                     ve_h = 0, ve_o = 0, age_sel = age_sel, share_via_O = share_hosp_via_outpt)
})
tp_vax <- lapply(seq_len(cycles), function(t) {
  make_tp_consistent(p_hosp_base, p_outpt, p_bg_mort, in_hosp_death,
                     ve_h = ve_hosp[t], ve_o = ve_outpt[t], age_sel = age_sel, share_via_O = share_hosp_via_outpt)
})

# --- simulate occupancy (used for incident counting) ---
cohort_size <- 100000
sim_arm <- function(tp_list, cohort_size=100000) {
  dist <- matrix(0, nrow=cycles+1, ncol=4, dimnames=list(0:cycles, c("H","O","I","D")))
  dist[1,"H"] <- cohort_size
  for (t in seq_len(cycles)) dist[t+1,] <- as.numeric(dist[t,] %*% tp_list[[t]])
  dist
}
dist_no <- sim_arm(tp_no_vax, cohort_size)
dist_vx <- sim_arm(tp_vax, cohort_size)

# ---------- EVENT-BASED COSTING & QALYS ----------
# durations (conservative)
dur_outpt_days       <- 5
dur_hosp_nonicu_days <- 6
dur_hosp_icu_days    <- 12

qaly_loss_outpt <- du_O * (dur_outpt_days / 365)
qaly_loss_hosp  <- (icu_pct * du_I_icu * (dur_hosp_icu_days / 365)) +
  ((1 - icu_pct) * du_I_nonicu * (dur_hosp_nonicu_days / 365))

count_incidents <- function(dist_t, m) {
  inc_O <- as.numeric(dist_t["H"] * m["H","O"])
  inc_I <- as.numeric(dist_t["H"] * m["H","I"] + dist_t["O"] * m["O","I"])
  list(inc_O = inc_O, inc_I = inc_I)
}

get_incident_series <- function(tp_list, cohort_size = 100000) {
  dist <- matrix(0, nrow = cycles + 1, ncol = 4, dimnames = list(0:cycles, c("H","O","I","D")))
  dist[1,"H"] <- cohort_size
  inc_O <- numeric(cycles); inc_I <- numeric(cycles)
  for (t in seq_len(cycles)) {
    inc <- count_incidents(dist[t,], tp_list[[t]])
    inc_O[t] <- inc$inc_O; inc_I[t] <- inc$inc_I
    dist[t+1,] <- as.numeric(dist[t,] %*% tp_list[[t]])
  }
  list(dist = dist, inc_O = inc_O, inc_I = inc_I)
}

res_no <- get_incident_series(tp_no_vax, cohort_size)
res_vx <- get_incident_series(tp_vax, cohort_size)

# diagnostics: incidence per 1,000 person-years & prevented events
total_py <- cohort_size * cycles
hosp_py_no <- sum(res_no$inc_I)  / total_py * 1000
hosp_py_vx <- sum(res_vx$inc_I)  / total_py * 1000
outp_py_no <- sum(res_no$inc_O)  / total_py * 1000
outp_py_vx <- sum(res_vx$inc_O)  / total_py * 1000
prevented_hosp_per1000 <- (sum(res_no$inc_I) - sum(res_vx$inc_I)) / cohort_size * 1000
prevented_outp_per1000 <- (sum(res_no$inc_O) - sum(res_vx$inc_O)) / cohort_size * 1000

cat("\n--- INCIDENCE DIAGNOSTICS ---\n")
print(list(
  hosp_per_1000py_no   = hosp_py_no,
  hosp_per_1000py_vax  = hosp_py_vx,
  outpt_per_1000py_no  = outp_py_no,
  outpt_per_1000py_vax = outp_py_vx,
  prevented_hosp_per1000_pp  = prevented_hosp_per1000,
  prevented_outpt_per1000_pp = prevented_outp_per1000
))

# discount factors
disc <- (1 / (1 + disc_rate))^(0:cycles)

# costs — use conservative ICU multiplier
icu_multiplier_cost <- 1.2
cost_I <- icu_pct*hosp_cost_mean*icu_multiplier_cost + (1-icu_pct)*hosp_cost_mean

costs_no <- c(0, res_no$inc_O * outpatient_cost + res_no$inc_I * cost_I)
costs_vx <- c(0, res_vx$inc_O * outpatient_cost + res_vx$inc_I * cost_I)
costs_vx[1] <- costs_vx[1] + vax_total_cost * cohort_size  # vaccine once at start

# baseline QALYs equal in both arms; event losses differ
base_qalys <- rep(util_base * cohort_size, cycles + 1)
event_qaly_loss_no <- c(0, res_no$inc_O * qaly_loss_outpt + res_no$inc_I * qaly_loss_hosp)
event_qaly_loss_vx <- c(0, res_vx$inc_O * qaly_loss_outpt + res_vx$inc_I * qaly_loss_hosp)

qalys_no <- base_qalys - event_qaly_loss_no
qalys_vx <- base_qalys - event_qaly_loss_vx

# discounted totals
dcost_no <- sum(costs_no * disc)
dcost_vx <- sum(costs_vx * disc)
dqaly_no <- sum(qalys_no * disc)
dqaly_vx <- sum(qalys_vx * disc)

# per-person outcomes & ICER
pp_cost_no <- dcost_no / cohort_size
pp_cost_vx <- dcost_vx / cohort_size
pp_qaly_no <- dqaly_no / cohort_size
pp_qaly_vx <- dqaly_vx / cohort_size

incr_cost <- pp_cost_vx - pp_cost_no
incr_qaly <- pp_qaly_vx - pp_qaly_no
icer <- incr_cost / incr_qaly

cat("\n--- EVENT-BASED CEA SUMMARY ---\n")
print(list(
  age_sel = age_sel,
  p_hosp_base = p_hosp_base,
  p_outpt = p_outpt,
  p_bg_mort = p_bg_mort,
  util_base = util_base,
  outpatient_cost = outpatient_cost,
  hosp_cost_mean = hosp_cost_mean,
  vax_total_cost = vax_total_cost,
  qaly_loss_outpt = qaly_loss_outpt,
  qaly_loss_hosp  = qaly_loss_hosp,
  incremental_cost_pp = incr_cost,
  incremental_qalys_pp = incr_qaly,
  icer = icer
))

# ================================================
# 07_psa_ceac_evppi.R  (updated to latest model)
# ================================================

source(file.path(DIR.PROJ,"R","01_setup_packages.R"))
source(file.path(DIR.PROJ,"R","99_helpers.R"))
make_dirs()

# --- Define output dir if not already ---
DIR.OUT <- file.path(DIR.PROJ, "outputs")
if (!dir.exists(DIR.OUT)) dir.create(DIR.OUT, recursive = TRUE)

# --- Function to run the Markov model with parameter overrides ----
run_markov <- function(
    p_hosp_base = 0.002583,
    p_outpt = 0.02478,
    p_bg_mort = 0.02,
    icu_pct = 0.191,
    in_hosp_death = 0.058,
    util_base = 0.811,
    outpatient_cost = 200,
    hosp_cost_mean = 22000,
    vax_total_cost = 295,
    du_O = 0.012,
    du_I_nonicu = 0.13,
    du_I_icu = 0.25,
    ve_shift = 0,      # shifts VE curve up/down
    inc_mult = 1       # multiplies p_hosp_base and p_outpt
) {
  
  # --- Waning VE arrays (apply shift multiplicatively) ---
  ve_file <- file.path(DIR.CLEAN, "ve_waning_points.csv")
  if (file.exists(ve_file)) {
    ve_pts <- readr::read_csv(ve_file)
    waning_hosp <- make_waning_fn(ve_pts$months[ve_pts$outcome=="hosp"],
                                  ve_pts$ve[ve_pts$outcome=="hosp"] + ve_shift)
    waning_outpt <- make_waning_fn(ve_pts$months[ve_pts$outcome=="outpt"],
                                   ve_pts$ve[ve_pts$outcome=="outpt"] + ve_shift)
  } else {
    waning_hosp  <- function(m) rep(0.82+ve_shift, length(m))
    waning_outpt <- function(m) rep(0.50+ve_shift, length(m))
  }
  
  cycles    <- 15
  disc_rate <- 0.03
  ve_months <- seq(6, cycles*12, by=12)
  ve_hosp   <- pmax(0, waning_hosp(ve_months))
  ve_outpt  <- pmax(0, waning_outpt(ve_months))
  
  # Adjust incidence by multiplier
  p_hosp    <- p_hosp_base * inc_mult
  p_out     <- p_outpt * inc_mult
  
  tp_no <- lapply(seq_len(cycles), function(t)
    make_tp_consistent(p_hosp = p_hosp,
                       p_outpt = p_out,
                       p_bg_mort = p_bg_mort,
                       p_case_fatal = in_hosp_death,
                       ve_h = 0, ve_o = 0,
                       share_via_O = 0.40))
  tp_vx <- lapply(seq_len(cycles), function(t)
    make_tp_consistent(p_hosp = p_hosp,
                       p_outpt = p_out,
                       p_bg_mort = p_bg_mort,
                       p_case_fatal = in_hosp_death,
                       ve_h = ve_hosp[t],
                       ve_o = ve_outpt[t],
                       share_via_O = 0.40))
  
  # --- Simulate and count incidents ---
  get_incident_series <- function(tp_list, cohort = 1e5){
    dist <- matrix(0, nrow=cycles+1, ncol=4, dimnames=list(0:cycles,c("H","O","I","D")))
    dist[1,"H"] <- cohort
    inc_O <- inc_I <- numeric(cycles)
    for (t in seq_len(cycles)) {
      inc_O[t] <- dist[t,"H"]*tp_list[[t]]["H","O"]
      inc_I[t] <- dist[t,"H"]*tp_list[[t]]["H","I"] + dist[t,"O"]*tp_list[[t]]["O","I"]
      dist[t+1,] <- as.numeric(dist[t,] %*% tp_list[[t]])
    }
    list(dist=dist, inc_O=inc_O, inc_I=inc_I)
  }
  
  res_no <- get_incident_series(tp_no)
  res_vx <- get_incident_series(tp_vx)
  
  disc <- (1/(1+disc_rate))^(0:cycles)
  
  qaly_loss_outpt <- du_O * (7/365)
  qaly_loss_hosp  <- (icu_pct*du_I_icu*(14/365)) + ((1-icu_pct)*du_I_nonicu*(7/365))
  
  # Costs
  costs_no <- c(0, res_no$inc_O*outpatient_cost +
                  res_no$inc_I*((icu_pct*1.5 + (1-icu_pct))*hosp_cost_mean))
  costs_vx <- c(0, res_vx$inc_O*outpatient_cost +
                  res_vx$inc_I*((icu_pct*1.5 + (1-icu_pct))*hosp_cost_mean))
  costs_vx[1] <- costs_vx[1] + vax_total_cost*1e5
  
  # QALYs
  base_q <- rep(util_base*1e5, cycles+1)
  qloss_no <- c(0, res_no$inc_O*qaly_loss_outpt + res_no$inc_I*qaly_loss_hosp)
  qloss_vx <- c(0, res_vx$inc_O*qaly_loss_outpt + res_vx$inc_I*qaly_loss_hosp)
  
  dcost_no <- sum(costs_no*disc)
  dcost_vx <- sum(costs_vx*disc)
  dqaly_no <- sum((base_q-qloss_no)*disc)
  dqaly_vx <- sum((base_q-qloss_vx)*disc)
  
  data.frame(
    incr_cost = (dcost_vx - dcost_no)/1e5,
    incr_qaly = (dqaly_vx - dqaly_no)/1e5
  )
}

# --- PSA draws ---
set.seed(123)
n_sim <- 2000
WTP   <- c(50000,100000,150000)

hosp_cost_mean <- 22000; hosp_cost_sd <- 8000
shape <- (hosp_cost_mean/hosp_cost_sd)^2
rate  <- hosp_cost_mean/hosp_cost_sd^2
hosp_cost_draws <- rgamma(n_sim, shape=shape, rate=rate)

ve_shift_draws  <- rnorm(n_sim, mean=0, sd=0.05)
inc_mult_draws  <- rlnorm(n_sim, meanlog=0, sdlog=0.2)

res <- purrr::map2_df(hosp_cost_draws, seq_len(n_sim), function(hc, i){
  ve_s <- ve_shift_draws[i]; inc_m <- inc_mult_draws[i]
  out <- run_markov(hosp_cost_mean = hc,
                    ve_shift = ve_s,
                    inc_mult = inc_m)
  out
})
# --- CEAC ---
ceac_probs <- purrr::map_dbl(WTP, function(lambda) {
  mean(res$incr_cost - lambda*res$incr_qaly < 0)
})

ceac <- data.frame(
  WTP = WTP,
  Prob_CE = ceac_probs
)

readr::write_csv(ceac, file.path(DIR.OUT,"ceac.csv"))

# Plot
p <- ggplot(ceac, aes(x = WTP, y = Prob_CE)) +
  geom_line() + geom_point() +
  labs(x = "Willingness-to-pay ($/QALY)",
       y = "Probability Cost-Effective") +
  theme_minimal()

print(p)
ggsave(file.path(DIR.OUT, "ceac_plot.png"), p, width = 8, height = 6, dpi = 300)
