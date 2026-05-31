# =============================================================================
# Predicting Commercial-Carrier Safety Violations at National Scale
# -----------------------------------------------------------------------------
# A machine-learning analysis of US FMCSA Safety Measurement System data:
#   Q1 — Predict which carriers record Out-of-Service (OOS) violations
#   Q2 — Segment the carrier population into distinct risk profiles (K-Means)
#
# DATA (public): FMCSA SMS — Census | Inspection | Violations | Crash
#   Join key: dot_number  (Violations -> Inspection via inspection_id first)
#   Scale: 90M+ raw rows across four files; each raw file is ~1+ GB.
#   NOTE: This pipeline requires a high-memory machine (32GB+ RAM recommended)
#         or a cloud environment. Raw CSVs are not included in this repo due to
#         size — download them from the FMCSA SMS public data portal.
#
# PIPELINE
#   1. Load raw CSVs
#   2. Clean each table (Census, Inspection, Violations, Crash)
#   3. Aggregate to one row per carrier (dot_number)
#   4. Master merge onto Census (left joins, integrity-checked)
#   5. Feature engineering (OOS rate, violation-intensity score, ratios)
#   6. Target definition (binary OOS flag)
#   7. EDA (data quality, distributions, per-BASIC behaviour, geography)
#   8. Q1 — Random Forest: Strict (8 feat) vs Behavioral Fingerprint (17 feat)
#   9. Q2 — K-Means clustering with elbow + silhouette selection
#
# METHODOLOGICAL NOTE ON LEAKAGE (stated transparently):
#   Model B's rate features derive from the same inspection records that produce
#   the OOS target, so they carry weak target echoes. This is mitigated by
#   (a) using rates not raw counts, (b) framing Model B as descriptive/profiling,
#   and (c) always reporting it alongside the leakage-free Model A baseline.
# =============================================================================


# ── LIBRARIES ────────────────────────────────────────────────────────────────
library(tidyverse)      # dplyr, ggplot2, tidyr, readr, purrr
library(lubridate)      # date parsing
library(janitor)        # clean_names()
library(naniar)         # missingness visualisation
library(scales)         # axis formatting
library(ggcorrplot)     # correlation heatmaps
library(patchwork)      # multi-panel composites
library(ggridges)       # ridge density plots
library(hexbin)         # hex binning for dense scatterplots
library(ranger)         # fast random forest
library(pROC)           # ROC + threshold optimisation
library(caret)          # confusion matrix, fold creation
library(cluster)        # silhouette
library(factoextra)     # cluster visualisation


# ── GLOBAL HELPERS & THEME ───────────────────────────────────────────────────
N_CORES <- max(1L, parallel::detectCores() - 1L)

# Safe division: returns 0 when denominator is missing or zero
safe_div <- function(x, y) ifelse(is.na(y) | y == 0, 0, x / y)

# Robust Y/N -> integer parser used across all four raw tables
yn_to_int <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "NULL", "N/A")] <- NA_character_
  as.integer(x %in% c("Y", "YES", "1", "T", "TRUE"))
}

to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

# Fast Mann-Whitney AUC — ~10-50x faster than pROC::auc() on large test sets.
# as.numeric() casts both operands to double so n_pos * n_neg does not overflow
# R's 2.1-billion integer cap on large test sets.
fast_auc <- function(labels, scores, positive = "OOS") {
  pos <- scores[labels == positive]
  neg <- scores[labels != positive]
  n_pos <- as.numeric(length(pos))
  n_neg <- as.numeric(length(neg))
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  r <- rank(c(pos, neg))
  (sum(r[seq_len(n_pos)]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

theme_eda <- function() {
  theme_bw(base_size = 11) +
    theme(plot.title       = element_text(face = "bold"),
          plot.subtitle    = element_text(color = "grey30"),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92"),
          strip.text       = element_text(face = "bold"))
}

PAL_OOS    <- c("Clean" = "#16A34A", "OOS" = "#DC2626")
PAL_MODELS <- c("Strict (8 feat)" = "#F97316", "Fingerprint (17 feat)" = "#1F4E79")
PAL_OPS    <- c("A" = "#1f77b4", "B" = "#d62728", "C" = "#7f7f7f")
PAL_HM     <- c("Hazmat" = "#d73027", "Non-Hazmat" = "#bababa")
OP_LABELS  <- c("A" = "Interstate", "B" = "Intrastate Hazmat",
                "C" = "Intrastate Non-Hazmat")

FIG_DIR <- "figs"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, name, w = 8, h = 5) {
  ggsave(file.path(FIG_DIR, paste0(name, ".png")),
         plot = p, width = w, height = h, dpi = 150, bg = "white")
}


# STEP 1 — LOAD RAW DATA

cat("\n── Loading raw CSVs ──\n")

census_raw <- read_csv("data/SMS_Input_-_Motor_Carrier_Census_Information_20260315.csv") |> clean_names()

insp_raw  <- read.csv('data/SMS_Input_-_Inspection_20260315.csv') |> clean_names()

viol_raw  <- read.csv('data/Vehicle_Inspections_and_Violations_20260316.csv') |> clean_names()

crash_raw  <- read_csv("data/SMS_Input_-_Crash_20260315.csv") |> clean_names()

# Quick shape report
cat(sprintf("Census     : %d rows × %d cols\n", nrow(census_raw), ncol(census_raw)))
cat(sprintf("Inspection : %d rows × %d cols\n", nrow(insp_raw),   ncol(insp_raw)))
cat(sprintf("Violations : %d rows × %d cols\n", nrow(viol_raw),   ncol(viol_raw)))
cat(sprintf("Crash      : %d rows × %d cols\n", nrow(crash_raw),  ncol(crash_raw)))


#  2. CLEANED EACH TABLE INDIVIDUALLY 

# 2a) Census: master table, one row per carrier
census <- census_raw %>%
  mutate(
    dot_number          = trimws(dot_number),
    nbr_power_unit      = to_num(nbr_power_unit),
    driver_total        = to_num(driver_total),
    mcs150_mileage      = to_num(mcs150_mileage),
    recent_mileage      = to_num(recent_mileage),
    hm_flag             = yn_to_int(hm_flag),
    pc_flag             = yn_to_int(pc_flag),
    private_only        = yn_to_int(private_only),
    authorized_for_hire = yn_to_int(authorized_for_hire),
    carrier_operation   = toupper(trimws(carrier_operation))
  ) %>%
  filter(!is.na(dot_number), dot_number != "") %>%
  distinct(dot_number, .keep_all = TRUE)   # safety: master must be unique

# 2b) Inspections
inspection <- insp_raw %>%
  mutate(
    dot_number = trimws(dot_number),
    unique_id  = trimws(unique_id),
    insp_date  = suppressWarnings(ymd(insp_date)),  # change to mdy() if needed
    across(c(driver_oos_total, vehicle_oos_total, oos_total, hazmat_oos_total,
             total_hazmat_sent, basic_viol,
             unsafe_viol, fatigued_viol, dr_fitness_viol,
             subt_alcohol_viol, vh_maint_viol, hm_viol,
             time_weight),
           to_num),
    across(c(unsafe_insp, fatigued_insp, dr_fitness_insp,
             subt_alcohol_insp, vh_maint_insp, hm_insp,
             hazmat_placard_req),
           yn_to_int)
  ) %>%
  filter(!is.na(dot_number), dot_number != "")

# 2c) Violations - line-level, FK = inspection_id -> Inspection.unique_id
violation <- viol_raw %>%
  mutate(
    inspection_id      = trimws(inspection_id),
    out_of_service_flag = yn_to_int(out_of_service_indicator)
  ) %>%
  filter(!is.na(inspection_id), inspection_id != "")

# 2d) Crashes
crash <- crash_raw %>%
  mutate(
    dot_number      = trimws(dot_number),
    fatalities      = to_num(fatalities),
    injuries        = to_num(injuries),
    severity_weight = to_num(severity_weight),
    tow_away        = yn_to_int(tow_away),
    hazmat_released = yn_to_int(hazmat_released),
    not_preventable = yn_to_int(not_preventable),
    report_date     = suppressWarnings(ymd(report_date))
  ) %>%
  filter(!is.na(dot_number), dot_number != "") %>%
  filter(is.na(not_preventable) | not_preventable == 0)   # SMS exclusion


#  3. AGGREGATED TO CARRIER LEVEL 
# 3a) Inspection aggregates per dot_number
insp_agg <- inspection %>%
  group_by(dot_number) %>%
  summarise(
    total_inspections = n(),
    sum_oos_total     = sum(oos_total,         na.rm = TRUE),
    sum_driver_oos    = sum(driver_oos_total,  na.rm = TRUE),
    sum_vehicle_oos   = sum(vehicle_oos_total, na.rm = TRUE),
    sum_hazmat_oos    = sum(hazmat_oos_total,  na.rm = TRUE),
    sum_basic_viol    = sum(basic_viol,        na.rm = TRUE),
    viol_unsafe       = sum(unsafe_viol,       na.rm = TRUE),
    viol_fatigued     = sum(fatigued_viol,     na.rm = TRUE),
    viol_dr_fitness   = sum(dr_fitness_viol,   na.rm = TRUE),
    viol_subt_alcohol = sum(subt_alcohol_viol, na.rm = TRUE),
    viol_vh_maint     = sum(vh_maint_viol,     na.rm = TRUE),
    viol_hm           = sum(hm_viol,           na.rm = TRUE),
    insp_unsafe       = sum(unsafe_insp,       na.rm = TRUE),
    insp_fatigued     = sum(fatigued_insp,     na.rm = TRUE),
    insp_dr_fitness   = sum(dr_fitness_insp,   na.rm = TRUE),
    insp_subt_alcohol = sum(subt_alcohol_insp, na.rm = TRUE),
    insp_vh_maint     = sum(vh_maint_insp,     na.rm = TRUE),
    insp_hm           = sum(hm_insp,           na.rm = TRUE),
    last_inspection   = suppressWarnings(max(insp_date, na.rm = TRUE)),
    .groups = "drop"
  )

# 3b) Violation aggregates - bridge through Inspection to get dot_number
viol_with_dot <- violation %>%
  inner_join(inspection %>% select(unique_id, dot_number),
             by = c("inspection_id" = "unique_id"))

viol_agg <- viol_with_dot %>%
  group_by(dot_number) %>%
  summarise(
    total_violations_recorded = n(),
    total_oos_violations      = sum(out_of_service_flag, na.rm = TRUE),
    .groups = "drop"
  )

# 3c) Crash aggregates per dot_number
crash_agg <- crash %>%
  group_by(dot_number) %>%
  summarise(
    total_crashes        = n(),
    total_fatalities     = sum(fatalities,      na.rm = TRUE),
    total_injuries       = sum(injuries,        na.rm = TRUE),
    total_tow_away       = sum(tow_away,        na.rm = TRUE),
    total_hazmat_release = sum(hazmat_released, na.rm = TRUE),
    sum_severity_weight  = sum(severity_weight, na.rm = TRUE),
    .groups = "drop"
  )


# 4. MASTER MERGE 
# Census is the LEFT base. Each child table is already aggregated to one row
# per dot_number, so the result has exactly nrow(census) rows.

master <- census %>%
  left_join(insp_agg,  by = "dot_number") %>%
  left_join(viol_agg,  by = "dot_number") %>%
  left_join(crash_agg, by = "dot_number") %>%
  mutate(across(c(total_inspections, sum_oos_total, sum_driver_oos,
                  sum_vehicle_oos, sum_hazmat_oos, sum_basic_viol,
                  viol_unsafe, viol_fatigued, viol_dr_fitness,
                  viol_subt_alcohol, viol_vh_maint, viol_hm,
                  insp_unsafe, insp_fatigued, insp_dr_fitness,
                  insp_subt_alcohol, insp_vh_maint, insp_hm,
                  total_violations_recorded, total_oos_violations,
                  total_crashes, total_fatalities, total_injuries,
                  total_tow_away, total_hazmat_release,
                  sum_severity_weight),
                ~ replace_na(., 0)))

# Hard sanity check - if this fails, the merge fan-out happened.
stopifnot(nrow(master) == nrow(census))
cat(sprintf("Master rows = %d  |  Census rows = %d  ->  MERGE OK\n",
            nrow(master), nrow(census)))


# =============================================================================
# MISSING-VALUE TREATMENT (documented, four-group strategy)
# =============================================================================
# NAs in this dataset do not all mean the same thing. We treat them in four
# groups, each with a different justification. (Event-count columns from the
# child tables were already set to 0 in the merge above — that is GROUP 1,
# made explicit here.)

master <- master %>%

  # GROUP 1 — Structural missingness ("safe" carriers).  ACTION: impute 0.
  # In a federal relational database, a carrier with no crashes / no violations /
  # never inspected returns NULL on join, not 0. Those NAs literally mean
  # "zero events", so 0 is the correct value, not a guess. (Applied in the merge
  # above; re-stated here for the engineered/ratio columns.)
  mutate(across(c(sum_driver_oos, sum_vehicle_oos, sum_hazmat_oos,
                  total_oos_violations, sum_basic_viol,
                  insp_unsafe, insp_fatigued, insp_dr_fitness,
                  insp_subt_alcohol, insp_vh_maint, insp_hm,
                  total_fatalities, total_injuries, sum_severity_weight,
                  total_tow_away, total_hazmat_release),
                ~ replace_na(., 0))) %>%

  # GROUP 2 — Core operational metrics.  ACTION: impute median.
  # An active commercial carrier physically cannot have 0 trucks or 0 drivers,
  # so 0 would be wrong here. The median (not mean) is used so mega-fleets like
  # FedEx/UPS with 100,000+ units do not skew the imputed value.
  mutate(across(c(nbr_power_unit, driver_total),
                ~ replace_na(., median(., na.rm = TRUE)))) %>%

  # GROUP 3 — Categorical / high-volume text.  ACTION: impute "UNKNOWN".
  # Fields like dba_name are missing for ~1.5M carriers simply because most
  # firms trade under their legal name. We keep the row and flag the absence so
  # downstream factor-based steps do not choke on NA. (any_of() guards against
  # columns that may be absent in a given data vintage.)
  mutate(across(any_of(c("dba_name", "legal_name", "vmt_source_id",
                         "op_other", "telephone")),
                ~ replace_na(as.character(.), "UNKNOWN"))) %>%

  # GROUP 4 — Minor geographic fields.  ACTION: drop the rows.
  # phy_state / phy_city etc. are missing for ~1,500 of ~2M carriers (<0.1%).
  # Dropping is cleaner than inventing fake locations, especially since
  # geography is used as a feature later. (any_of via where() keeps this safe
  # if a column is absent.)
  filter(if_all(any_of(c("phy_state", "phy_city", "phy_zip")),
                ~ !is.na(.)))

cat(sprintf("After missing-value treatment: %d carriers\n", nrow(master)))


# 5. FEATURE ENGINEERING 
# BASIC weights for Violation Intensity Score - tune to your evidence.
W_UNSAFE     <- 1.0
W_FATIGUE    <- 1.0
W_DR_FITNESS <- 0.7
W_ALCOHOL    <- 1.5   # high crash-risk per FMCSA SMS
W_VH_MAINT   <- 0.8
W_HM         <- 1.2   # higher consequence for hazmat releases

master <- master %>%
  mutate(
    # Out-of-Service Rate
    oos_rate                  = safe_div(sum_oos_total, total_inspections),
    # Violation Intensity Score (weighted)
    violation_intensity_score = W_UNSAFE     * viol_unsafe       +
      W_FATIGUE    * viol_fatigued     +
      W_DR_FITNESS * viol_dr_fitness   +
      W_ALCOHOL    * viol_subt_alcohol +
      W_VH_MAINT   * viol_vh_maint     +
      W_HM         * viol_hm,
    # Crash-to-Inspection Ratio
    crash_to_inspection_ratio = safe_div(total_crashes, total_inspections),
    # Convenience
    fatal_crash_rate      = safe_div(total_fatalities, total_crashes),
    injury_per_inspection = safe_div(total_injuries,  total_inspections),
    drivers_per_unit      = safe_div(driver_total,    nbr_power_unit),
    fleet_size_log        = log1p(pmax(nbr_power_unit, 0, na.rm = TRUE)),
    # carrier_operation dummies (A=Interstate, B=Intra-Hazmat, C=Intra-NonHM)
    is_interstate   = as.integer(carrier_operation == "A"),
    is_intra_hazmat = as.integer(carrier_operation == "B"),
    is_intra_nonhaz = as.integer(carrier_operation == "C")
  )


# 6. TARGET DEFINITION (For Research Question 1) 
# Binary: did this carrier record any OOS violation across its inspections?
master <- master %>%
  mutate(oos_flag = factor(if_else(sum_oos_total > 0, "OOS", "Clean"),
                           levels = c("Clean", "OOS")))

cat("\nTarget distribution (full master):\n"); print(table(master$oos_flag))

# Restrict modelling to carriers with >=1 inspection (otherwise target is
# undefined - "no inspection" is not the same as "passed every inspection").
model_df <- master %>% filter(total_inspections >= 1)


# 7. EDA / DATA QUALITY CHECKS 
cat("\nMissing values per modelling column:\n")
print(sapply(model_df, function(v) sum(is.na(v))))

num_features <- model_df %>%
  select(nbr_power_unit, driver_total, total_inspections,
         oos_rate, violation_intensity_score,
         crash_to_inspection_ratio, total_crashes,
         total_fatalities, sum_severity_weight) %>%
  mutate(across(everything(), ~ replace_na(., 0)))

png(file.path(FIG_DIR, "corr_plot.png"), 900, 800)
corrplot::corrplot(cor(num_features), method = "color", type = "upper",
                   tl.cex = 0.8)
dev.off()

# ---------------------------------------------------------------------------
# CHECKPOINT: persist the carrier-level master table to disk.
# The 90M-row merge above is expensive; writing it here lets the EDA and
# modelling sections below resume from this point without recomputing it.
# ---------------------------------------------------------------------------
write_csv(master, "data/Final_file.csv")


# =============================================================================
# EXPLORATORY DATA ANALYSIS
# =============================================================================

# Reading the Data-set
df <- read.csv("data/Final_file.csv")
summary(df)
colnames(df)

#Fixing the missing values

missing_summary <- df %>%
  summarise_all(~ mean(is.na(.)) * 100) %>%
  pivot_longer(everything(), names_to = "column", values_to = "missing_pct") %>%
  filter(missing_pct > 0) %>%
  arrange(desc(missing_pct))

# View the percentages
print(missing_summary, n = Inf)


# 1. NORMALISE TYPES (df is assumed to already exist) 
to_lgl <- function(x) {
  if (is.logical(x)) return(x)
  s <- toupper(trimws(as.character(x)))
  out <- rep(NA, length(s))
  out[s %in% c("Y", "YES", "TRUE", "T", "1")] <- TRUE
  out[s %in% c("N", "NO",  "FALSE", "F", "0")] <- FALSE
  as.logical(out)
}

df <- df %>%
  mutate(
    hm_flag_lgl    = to_lgl(hm_flag),
    pc_flag_lgl    = to_lgl(pc_flag),
    HM_Status      = if_else(hm_flag_lgl, "Hazmat", "Non-Hazmat"),
    oos_flag       = factor(oos_flag, levels = c("Clean", "OOS")),
    has_inspection = total_inspections >= 1,
    # convenience derived feature for B2 - raw violations per inspection,
    # used to contrast against the engineered violation_intensity_score
    raw_viol_per_insp = if_else(total_inspections > 0,
                                sum_basic_viol / total_inspections, 0)
  )

# Sanity print: confirm key columns exist
needed <- c("nbr_power_unit", "driver_total", "carrier_operation", "hm_flag_lgl",
            "phy_state", "total_inspections", "sum_oos_total", "sum_basic_viol",
            "sum_driver_oos", "sum_vehicle_oos", "sum_hazmat_oos",
            "viol_unsafe", "viol_fatigued", "viol_dr_fitness",
            "viol_subt_alcohol", "viol_vh_maint", "viol_hm",
            "insp_unsafe", "insp_fatigued", "insp_dr_fitness",
            "insp_subt_alcohol", "insp_vh_maint", "insp_hm",
            "total_crashes", "total_fatalities", "total_injuries",
            "sum_severity_weight", "total_hazmat_release",
            "oos_rate", "violation_intensity_score",
            "crash_to_inspection_ratio", "oos_flag")
cat("Columns missing for EDA:",
    paste(setdiff(needed, names(df)), collapse = ", "), "\n")


# A. UNIVERSE & DATA QUALITY

# A1. Missing data percent (top 25)
a1 <- df %>%
  summarise(across(everything(), ~ mean(is.na(.)) * 100)) %>%
  pivot_longer(everything(), names_to = "column", values_to = "missing_pct") %>%
  filter(missing_pct > 0) %>%
  arrange(desc(missing_pct)) %>%
  slice_head(n = 25) %>%
  ggplot(aes(missing_pct, reorder(column, missing_pct))) +
  geom_col(fill = "#1F4E79") +
  geom_text(aes(label = sprintf("%.1f%%", missing_pct)),
            hjust = -0.1, size = 3) +
  scale_x_continuous(expand = expansion(c(0, 0.15))) +
  labs(title = "Top 25 columns by missing-value share",
       subtitle = "Confirms the imputation map applied earlier",
       x = "% missing", y = NULL) +
  theme_eda()
save_plot(a1, "A1_missing_data", h = 7)

# A2. Carrier operation distribution
a2 <- df %>%
  count(carrier_operation) %>%
  ggplot(aes(reorder(carrier_operation, -n), n, fill = carrier_operation)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = comma(n)), vjust = -0.4) +
  scale_x_discrete(labels = OP_LABELS) +
  scale_y_continuous(labels = comma, expand = expansion(c(0, 0.1))) +
  scale_fill_manual(values = PAL_OPS) +
  labs(title = "Active carriers by operation type",
       subtitle = "Interstate dominates the universe - any model must handle this imbalance",
       x = NULL, y = "Number of carriers") +
  theme_eda()
save_plot(a2, "A2_carrier_operation_distribution")

# A3. HM and Passenger flags by operation
a3_data <- bind_rows(
  df %>% filter(hm_flag_lgl)  %>% count(carrier_operation) %>% mutate(Type = "Hazmat"),
  df %>% filter(pc_flag_lgl)  %>% count(carrier_operation) %>% mutate(Type = "Passenger")
)
a3 <- ggplot(a3_data,
             aes(reorder(carrier_operation, n), n, fill = Type)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = comma(n)),
            position = position_dodge(width = 0.9),
            hjust = -0.1, size = 3) +
  coord_flip() +
  scale_x_discrete(labels = OP_LABELS) +
  scale_y_continuous(labels = comma, expand = expansion(c(0, 0.18))) +
  scale_fill_manual(values = c("Hazmat" = "#d73027", "Passenger" = "#4575b4")) +
  labs(title = "Hazmat and Passenger carriers across operation types",
       subtitle = "Hazmat concentrates in interstate; passenger is rarer overall",
       x = NULL, y = "Carrier count") +
  theme_eda()
save_plot(a3, "A3_hm_pax_by_operation")

# A4. Fleet size distribution
a4 <- df %>%
  filter(nbr_power_unit > 0) %>%
  ggplot(aes(nbr_power_unit)) +
  geom_histogram(bins = 60, fill = "#1F4E79", color = "white") +
  scale_x_log10(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(title = "Fleet size distribution (log10 scale)",
       subtitle = "Long tail - tiny owner-operators dominate; mega fleets are rare but matter",
       x = "Power units (log10)", y = "Carriers") +
  theme_eda()
save_plot(a4, "A4_fleet_size_distribution")

# A5. Drivers per power unit (uses pre-computed drivers_per_unit)
a5 <- df %>%
  filter(drivers_per_unit > 0, drivers_per_unit < 10) %>%
  ggplot(aes(drivers_per_unit)) +
  geom_histogram(bins = 60, fill = "#7F1D1D", color = "white") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(title = "Drivers per power unit",
       subtitle = "Bunched around 1; values >1 imply multi-shift / team operations",
       x = "Drivers / power units", y = "Carriers") +
  theme_eda()
save_plot(a5, "A5_drivers_per_unit")


# B. QUESTION 1 - OOS VIOLATION PREDICTION

# B1. TARGET VARIABLE EXPLORATION 
# B1a. OOS rate buckets
b1a <- df %>%
  filter(has_inspection) %>%
  mutate(oos_bucket = cut(oos_rate,
                          breaks = c(-0.01, 0, 0.1, 0.25, 0.5, 1.01),
                          labels = c("0%", "0-10%", "10-25%", "25-50%", ">50%"))) %>%
  count(oos_bucket) %>% drop_na() %>%
  mutate(pct = n / sum(n)) %>%
  ggplot(aes(oos_bucket, pct, fill = oos_bucket)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = percent(pct, 0.1)), vjust = -0.4) +
  scale_y_continuous(labels = percent, expand = expansion(c(0, 0.15))) +
  scale_fill_manual(values = c("#16A34A","#84CC16","#F59E0B","#F97316","#DC2626")) +
  labs(title = "Distribution of carrier Out-of-Service rates",
       subtitle = "Most inspected carriers sit at 0% - but the right tail is the policy target",
       x = "OOS rate bucket", y = "Share of inspected carriers") +
  theme_eda()
save_plot(b1a, "B1a_oos_rate_distribution")

# B1b. Class balance of binary target
b1b <- df %>%
  filter(has_inspection) %>%
  count(oos_flag) %>%
  mutate(pct = n / sum(n)) %>%
  ggplot(aes(oos_flag, pct, fill = oos_flag)) +
  geom_col(show.legend = FALSE, width = 0.6) +
  geom_text(aes(label = sprintf("%s (%s)", comma(n), percent(pct, 0.1))),
            vjust = -0.4) +
  scale_y_continuous(labels = percent, expand = expansion(c(0, 0.15))) +
  scale_fill_manual(values = PAL_OOS) +
  labs(title = "Binary target: did the carrier ever record an OOS violation?",
       subtitle = "Class imbalance drives the threshold-tuning question - keep this in mind",
       x = NULL, y = "Share of inspected carriers") +
  theme_eda()
save_plot(b1b, "B1b_oos_class_balance")

# B1c. OOS rate by operation (positive-rate boxplot)
b1c <- df %>%
  filter(has_inspection, oos_rate > 0) %>%
  ggplot(aes(carrier_operation, oos_rate, fill = carrier_operation)) +
  geom_boxplot(outlier.alpha = 0.2, show.legend = FALSE) +
  scale_x_discrete(labels = OP_LABELS) +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(values = PAL_OPS) +
  labs(title = "OOS rate by operation type (carriers with at least one OOS)",
       subtitle = "Spread differs across operation types - feature is informative",
       x = NULL, y = "OOS rate") +
  theme_eda()
save_plot(b1c, "B1c_oos_rate_by_operation")

# B1d. OOS prevalence by HM status
b1d <- df %>%
  filter(has_inspection) %>%
  group_by(HM_Status) %>%
  summarise(oos_share = mean(oos_flag == "OOS"), n = n(), .groups = "drop") %>%
  ggplot(aes(HM_Status, oos_share, fill = HM_Status)) +
  geom_col(show.legend = FALSE, width = 0.55) +
  geom_text(aes(label = sprintf("%s\n(n=%s)", percent(oos_share, 0.1), comma(n))),
            vjust = -0.2) +
  scale_y_continuous(labels = percent, expand = expansion(c(0, 0.18))) +
  scale_fill_manual(values = PAL_HM) +
  labs(title = "OOS prevalence: Hazmat vs Non-Hazmat",
       subtitle = "Difference here motivates the Hazmat sub-model in Q1's third support question",
       x = NULL, y = "Share with at least one OOS") +
  theme_eda()
save_plot(b1d, "B1d_oos_by_hm")


# B2. ENGINEERED-FEATURE IMPACT 
# This is the heart of "do engineered features improve the model?".

# B2a. oos_rate density by class (sanity check)
b2a <- df %>%
  filter(has_inspection) %>%
  ggplot(aes(oos_rate, fill = oos_flag, color = oos_flag)) +
  geom_density(alpha = 0.4) +
  scale_x_continuous(labels = percent, limits = c(0, 1)) +
  scale_fill_manual(values = PAL_OOS) +
  scale_color_manual(values = PAL_OOS) +
  labs(title = "OOS rate density by OOS class",
       subtitle = "Sanity check: the engineered rate cleanly separates the two classes",
       x = "OOS rate", y = "Density", fill = NULL, color = NULL) +
  theme_eda() + theme(legend.position = "top")
save_plot(b2a, "B2a_oos_rate_density")

# B2b. Crash-to-Inspection ratio density
b2b <- df %>%
  filter(has_inspection, crash_to_inspection_ratio < 1) %>%
  ggplot(aes(crash_to_inspection_ratio, fill = oos_flag, color = oos_flag)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = PAL_OOS) +
  scale_color_manual(values = PAL_OOS) +
  labs(title = "Crash-to-Inspection ratio by OOS class",
       subtitle = "If OOS carriers have a heavier right tail, the feature carries signal",
       x = "Crashes / Inspections", y = "Density", fill = NULL, color = NULL) +
  theme_eda() + theme(legend.position = "top")
save_plot(b2b, "B2b_crash_insp_ratio_density")

# B2c. Violation Intensity Score density - THE engineered feature for Q1
b2c <- df %>%
  filter(has_inspection, violation_intensity_score > 0) %>%
  ggplot(aes(violation_intensity_score, fill = oos_flag, color = oos_flag)) +
  geom_density(alpha = 0.4) +
  scale_x_log10(labels = comma) +
  scale_fill_manual(values = PAL_OOS) +
  scale_color_manual(values = PAL_OOS) +
  labs(title = "Violation Intensity Score by OOS class (log10 scale)",
       subtitle = "Right-shift for OOS class is the visual proof the score adds predictive signal",
       x = "Violation Intensity Score (log10)", y = "Density",
       fill = NULL, color = NULL) +
  theme_eda() + theme(legend.position = "top")
save_plot(b2c, "B2c_violation_intensity_density")

# B2d. Raw vs engineered comparison - side by side
b2d_raw <- df %>%
  filter(has_inspection, raw_viol_per_insp > 0, raw_viol_per_insp < 30) %>%
  ggplot(aes(raw_viol_per_insp, fill = oos_flag, color = oos_flag)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = PAL_OOS) +
  scale_color_manual(values = PAL_OOS) +
  labs(title = "Raw: violations / inspection",
       x = "Violations / Inspection", y = "Density", fill = NULL, color = NULL) +
  theme_eda() + theme(legend.position = "none")

b2d_eng <- df %>%
  filter(has_inspection, violation_intensity_score > 0) %>%
  ggplot(aes(violation_intensity_score, fill = oos_flag, color = oos_flag)) +
  geom_density(alpha = 0.4) +
  scale_x_log10(labels = comma) +
  scale_fill_manual(values = PAL_OOS) +
  scale_color_manual(values = PAL_OOS) +
  labs(title = "Engineered: weighted intensity score (log10)",
       x = "Violation Intensity Score (log10)", y = NULL, fill = NULL, color = NULL) +
  theme_eda() + theme(legend.position = "top")

b2d <- (b2d_raw | b2d_eng) +
  plot_annotation(
    title = "Raw vs engineered features: which separates the classes better?",
    subtitle = "Wider gap between green and red densities = stronger predictor"
  )
save_plot(b2d, "B2d_raw_vs_engineered", w = 11, h = 5)

# B2e. Joint hex - violation_intensity_score x crash_to_inspection_ratio
b2e <- df %>%
  filter(has_inspection, violation_intensity_score > 0,
         crash_to_inspection_ratio > 0, crash_to_inspection_ratio < 1) %>%
  ggplot(aes(violation_intensity_score, crash_to_inspection_ratio)) +
  geom_hex(bins = 40) +
  scale_fill_viridis_c(trans = "log10", labels = comma) +
  scale_x_log10(labels = comma) +
  labs(title = "Joint distribution: Intensity Score vs Crash-Inspection Ratio",
       subtitle = "Carriers in the upper-right corner are highest-priority enforcement targets",
       x = "Violation Intensity Score (log10)",
       y = "Crashes / Inspections", fill = "Carriers (log)") +
  theme_eda()
save_plot(b2e, "B2e_joint_intensity_crash")


# B3. PER-BASIC BEHAVIOUR 
# B3a. Heatmap of mean OOS volumes by Driver/Vehicle/Hazmat split
b3a_data <- df %>%
  filter(has_inspection) %>%
  group_by(carrier_operation, HM_Status) %>%
  summarise(`Driver OOS`  = mean(sum_driver_oos,  na.rm = TRUE),
            `Vehicle OOS` = mean(sum_vehicle_oos, na.rm = TRUE),
            `Hazmat OOS`  = mean(sum_hazmat_oos,  na.rm = TRUE),
            .groups = "drop") %>%
  pivot_longer(c(`Driver OOS`, `Vehicle OOS`, `Hazmat OOS`),
               names_to = "OOS_Category", values_to = "Average_Violations")
b3a <- ggplot(b3a_data,
              aes(OOS_Category, carrier_operation, fill = Average_Violations)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", Average_Violations)),
            color = "grey20", size = 3) +
  scale_fill_gradient(low = "white", high = "#d73027", name = "Avg per carrier") +
  scale_y_discrete(labels = OP_LABELS) +
  facet_wrap(~ HM_Status, ncol = 1) +
  labs(title = "Where do OOS violations happen? (Driver / Vehicle / Hazmat split)",
       subtitle = "Hazmat-flagged carriers and interstate operations carry the bulk",
       x = "OOS category", y = NULL) +
  theme_eda() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"))
save_plot(b3a, "B3a_oos_split_heatmap", w = 8, h = 6)

# B3b. Stacked composition of OOS violations by operation
b3b <- df %>%
  group_by(carrier_operation) %>%
  summarise(Driver  = sum(sum_driver_oos,  na.rm = TRUE),
            Vehicle = sum(sum_vehicle_oos, na.rm = TRUE),
            Hazmat  = sum(sum_hazmat_oos,  na.rm = TRUE),
            .groups = "drop") %>%
  pivot_longer(-carrier_operation, names_to = "Category", values_to = "n") %>%
  group_by(carrier_operation) %>%
  mutate(share = n / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(carrier_operation, share, fill = Category)) +
  geom_col(position = "fill") +
  geom_text(aes(label = percent(share, 0.1)),
            position = position_fill(vjust = 0.5), color = "white", size = 3.2) +
  scale_x_discrete(labels = OP_LABELS) +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(values = c("Driver" = "#1f77b4",
                               "Vehicle" = "#ff7f0e",
                               "Hazmat" = "#d62728")) +
  labs(title = "Composition of OOS violations by operation",
       subtitle = "Driver-related vs vehicle-related mix shifts across operation types",
       x = NULL, y = "Share of OOS violations", fill = "OOS type") +
  theme_eda()
save_plot(b3b, "B3b_oos_composition")

# B3c. BASIC inspection coverage - share of carriers ever inspected per BASIC
b3c <- df %>%
  filter(has_inspection) %>%
  transmute(carrier_operation,
            `Unsafe Driving`    = insp_unsafe       > 0,
            `Hours-of-Service`  = insp_fatigued     > 0,
            `Driver Fitness`    = insp_dr_fitness   > 0,
            `Substance/Alcohol` = insp_subt_alcohol > 0,
            `Vehicle Maint.`    = insp_vh_maint     > 0,
            `HM Compliance`     = insp_hm           > 0) %>%
  pivot_longer(-carrier_operation, names_to = "BASIC", values_to = "covered") %>%
  group_by(carrier_operation, BASIC) %>%
  summarise(coverage = mean(covered, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(BASIC, carrier_operation, fill = coverage)) +
  geom_tile(color = "white") +
  geom_text(aes(label = percent(coverage, 1)), size = 3, color = "grey15") +
  scale_y_discrete(labels = OP_LABELS) +
  scale_fill_gradient(low = "white", high = "#1F4E79", labels = percent) +
  labs(title = "BASIC inspection coverage by operation type",
       subtitle = "Share of carriers ever inspected against each BASIC - feature-availability check",
       x = NULL, y = NULL, fill = "Coverage") +
  theme_eda() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(b3c, "B3c_basic_coverage", h = 5)

# B3d. Per-BASIC violation distribution by OOS class (uses viol_<basic>)
b3d_data <- df %>%
  filter(has_inspection) %>%
  select(oos_flag,
         `Unsafe Driving`    = viol_unsafe,
         `Hours-of-Service`  = viol_fatigued,
         `Driver Fitness`    = viol_dr_fitness,
         `Substance/Alcohol` = viol_subt_alcohol,
         `Vehicle Maint.`    = viol_vh_maint,
         `HM Compliance`     = viol_hm) %>%
  pivot_longer(-oos_flag, names_to = "BASIC", values_to = "viol") %>%
  filter(viol > 0)
b3d <- ggplot(b3d_data, aes(oos_flag, viol, fill = oos_flag)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.6) +
  scale_y_log10(labels = comma) +
  scale_fill_manual(values = PAL_OOS) +
  facet_wrap(~ BASIC, scales = "free_y") +
  labs(title = "Per-BASIC violation counts: OOS vs Clean carriers",
       subtitle = "Wider separation between green and red = more predictive BASIC",
       x = NULL, y = "Violations (log10)") +
  theme_eda() + theme(legend.position = "none")
save_plot(b3d, "B3d_viol_per_basic_by_class", w = 10, h = 6)


# B4. FEATURE-GROUP PREVIEW 

# B4a. Correlation matrix grouped
num_cols_for_cor <- df %>%
  select(any_of(c("nbr_power_unit", "driver_total", "recent_mileage",
                  "total_inspections", "sum_oos_total", "sum_basic_viol",
                  "violation_intensity_score",
                  "total_crashes", "total_fatalities",
                  "sum_severity_weight",
                  "oos_rate", "crash_to_inspection_ratio"))) %>%
  mutate(across(everything(), as.numeric))
cor_mat <- cor(num_cols_for_cor, use = "pairwise.complete.obs")
b4a <- ggcorrplot(cor_mat, type = "lower", lab = TRUE, lab_size = 3,
                  colors = c("#2563EB","white","#DC2626"),
                  outline.color = "white", tl.srt = 45) +
  labs(title = "Correlation between fleet, inspection, and crash metrics",
       subtitle = "Block structure: exposure (size/inspections), violations, and crashes cluster") +
  theme_eda() + theme(panel.grid = element_blank())
save_plot(b4a, "B4a_correlation_matrix", w = 8, h = 7)

# B4b. Top numeric features split by oos_flag - violins
b4b_data <- df %>%
  filter(has_inspection) %>%
  select(oos_flag, nbr_power_unit, driver_total, total_inspections,
         sum_basic_viol, total_crashes, violation_intensity_score) %>%
  pivot_longer(-oos_flag, names_to = "feature", values_to = "value") %>%
  filter(value > 0)
b4b <- ggplot(b4b_data, aes(oos_flag, value, fill = oos_flag)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.6) +
  scale_y_log10(labels = comma) +
  scale_fill_manual(values = PAL_OOS) +
  facet_wrap(~ feature, scales = "free_y") +
  labs(title = "How key features differ between OOS and Clean carriers",
       subtitle = "Log-y violins; broader violin = wider distribution at that magnitude",
       x = NULL, y = "Value (log10)") +
  theme_eda() + theme(legend.position = "none")
save_plot(b4b, "B4b_features_by_class", w = 9, h = 6)


# B5. HAZMAT vs NON-HAZMAT 

# B5a. OOS rate ridges - HM vs Non-HM by operation
b5a <- df %>%
  filter(has_inspection, oos_rate > 0) %>%
  ggplot(aes(oos_rate, carrier_operation, fill = HM_Status)) +
  geom_density_ridges(alpha = 0.6, scale = 1.05) +
  scale_x_continuous(labels = percent, limits = c(0, 1)) +
  scale_y_discrete(labels = OP_LABELS) +
  scale_fill_manual(values = PAL_HM) +
  labs(title = "OOS rate distribution by operation, split by Hazmat status",
       subtitle = "Where Hazmat ridges sit further right, Hazmat-specific risk is higher",
       x = "OOS rate", y = NULL, fill = NULL) +
  theme_eda() + theme(legend.position = "top")
save_plot(b5a, "B5a_oos_rate_ridges_by_hm", h = 6)

# B5b. Crash-to-Inspection ratio by HM
b5b <- df %>%
  filter(has_inspection, crash_to_inspection_ratio > 0,
         crash_to_inspection_ratio < 1) %>%
  ggplot(aes(crash_to_inspection_ratio, fill = HM_Status, color = HM_Status)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = PAL_HM) +
  scale_color_manual(values = PAL_HM) +
  labs(title = "Crash-to-Inspection ratio: Hazmat vs Non-Hazmat",
       subtitle = "Right-shift for Hazmat would mean each inspection covers a riskier carrier",
       x = "Crashes / Inspections", y = "Density",
       fill = NULL, color = NULL) +
  theme_eda() + theme(legend.position = "top")
save_plot(b5b, "B5b_crash_rate_by_hm")

# B5c. Severity weight by HM (log boxplot)
b5c <- df %>%
  filter(sum_severity_weight > 0) %>%
  ggplot(aes(HM_Status, sum_severity_weight, fill = HM_Status)) +
  geom_boxplot(outlier.alpha = 0.15, show.legend = FALSE) +
  scale_y_log10(labels = comma) +
  scale_fill_manual(values = PAL_HM) +
  labs(title = "Aggregate crash severity weight by Hazmat status",
       subtitle = "Log scale - tail mass is what to watch, not the median",
       x = NULL, y = "Total severity weight (log10)") +
  theme_eda()
save_plot(b5c, "B5c_severity_by_hm")


# C. QUESTION 2 - K-MEANS CLUSTERING SETUP - This is answering why clustering is necessary


# C1. MULTIVARIATE SIGNAL 

# C1a. Crash severity vs fleet size
c1a <- df %>%
  filter(nbr_power_unit > 0, total_crashes > 0) %>%
  ggplot(aes(nbr_power_unit, sum_severity_weight, color = HM_Status)) +
  geom_point(alpha = 0.4, size = 1.2) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  scale_x_log10(labels = comma) +
  scale_y_log10(labels = comma) +
  scale_color_manual(values = PAL_HM) +
  facet_wrap(~ carrier_operation, labeller = as_labeller(OP_LABELS)) +
  labs(title = "Crash severity vs fleet size across operation types",
       subtitle = "If size alone explained severity the slope would be ~1; deviation justifies clustering",
       x = "Power units (log10)", y = "Total severity weight (log10)",
       color = NULL) +
  theme_eda() + theme(legend.position = "top")
save_plot(c1a, "C1a_severity_vs_fleet_size", w = 10, h = 5)

# C1b. OOS rate vs Crash-to-Inspection ratio
c1b <- df %>%
  filter(has_inspection, oos_rate > 0,
         crash_to_inspection_ratio > 0, crash_to_inspection_ratio < 1) %>%
  ggplot(aes(oos_rate, crash_to_inspection_ratio, color = HM_Status)) +
  geom_point(alpha = 0.3, size = 1) +
  scale_x_continuous(labels = percent) +
  scale_color_manual(values = PAL_HM) +
  labs(title = "Two axes of risk: OOS rate vs Crash-to-Inspection ratio",
       subtitle = "Carriers in the upper-right corner are the cluster K-Means should isolate",
       x = "OOS rate", y = "Crashes / Inspections", color = NULL) +
  theme_eda() + theme(legend.position = "top")
save_plot(c1b, "C1b_oos_vs_crash_ratio")

# C1c. Drivers vs power units
c1c <- df %>%
  filter(driver_total > 0, nbr_power_unit > 0) %>%
  ggplot(aes(nbr_power_unit, driver_total)) +
  geom_hex(bins = 50) +
  scale_fill_viridis_c(trans = "log10", labels = comma) +
  scale_x_log10(labels = comma) + scale_y_log10(labels = comma) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Drivers vs power units (log-log)",
       subtitle = "Red line = 1 driver per unit; carriers above are multi-shift / team operations",
       x = "Power units (log10)", y = "Drivers (log10)", fill = "Carriers") +
  theme_eda()
save_plot(c1c, "C1c_drivers_vs_units")


#C2. CRASH & SEVERITY SHAPE 

# C2a. Severity weight log distribution
c2a <- df %>%
  filter(sum_severity_weight > 0) %>%
  ggplot(aes(sum_severity_weight)) +
  geom_histogram(bins = 60, fill = "#7F1D1D", color = "white") +
  scale_x_log10(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(title = "Aggregate crash severity weight per carrier",
       subtitle = "Heavy right tail: a small minority of carriers carry most severity",
       x = "Total severity weight (log10)", y = "Carriers") +
  theme_eda()
save_plot(c2a, "C2a_severity_distribution")

# C2b. Fatalities vs Injuries
c2b <- df %>%
  filter(total_crashes > 0) %>%
  ggplot(aes(total_injuries, total_fatalities)) +
  geom_jitter(width = 0.2, height = 0.2, alpha = 0.3, size = 1,
              color = "#7F1D1D") +
  scale_x_continuous(limits = c(0, quantile(df$total_injuries, 0.995, na.rm = TRUE))) +
  scale_y_continuous(limits = c(0, quantile(df$total_fatalities, 0.999, na.rm = TRUE))) +
  labs(title = "Carrier fatality count vs injury count",
       subtitle = "Points off the X-axis (fatalities > 0) are the highest-impact cluster",
       x = "Total injuries", y = "Total fatalities") +
  theme_eda()
save_plot(c2b, "C2b_fatalities_vs_injuries")

# C2c. Hazmat release rate by operation
c2c <- df %>%
  filter(total_crashes > 0) %>%
  group_by(carrier_operation, HM_Status) %>%
  summarise(release_rate = sum(total_hazmat_release, na.rm = TRUE) /
              sum(total_crashes, na.rm = TRUE),
            n_carriers = n(), .groups = "drop") %>%
  ggplot(aes(carrier_operation, release_rate, fill = HM_Status)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = percent(release_rate, 0.01)),
            position = position_dodge(width = 0.9),
            vjust = -0.3, size = 3) +
  scale_x_discrete(labels = OP_LABELS) +
  scale_y_continuous(labels = percent, expand = expansion(c(0, 0.15))) +
  scale_fill_manual(values = PAL_HM) +
  labs(title = "Hazmat release rate among crashes by operation type",
       subtitle = "Critical for the 'which cluster causes chemical spills' policy question",
       x = NULL, y = "Hazmat releases / crashes", fill = NULL) +
  theme_eda() + theme(legend.position = "top")
save_plot(c2c, "C2c_hazmat_release_rate")

# C2d. Per-BASIC OOS rate heatmap by operation - foreshadows clustering features
c2d <- df %>%
  filter(has_inspection) %>%
  group_by(carrier_operation) %>%
  summarise(`Unsafe Driving`    = sum(viol_unsafe)       / sum(total_inspections),
            `Hours-of-Service`  = sum(viol_fatigued)     / sum(total_inspections),
            `Driver Fitness`    = sum(viol_dr_fitness)   / sum(total_inspections),
            `Substance/Alcohol` = sum(viol_subt_alcohol) / sum(total_inspections),
            `Vehicle Maint.`    = sum(viol_vh_maint)     / sum(total_inspections),
            `HM Compliance`     = sum(viol_hm)           / sum(total_inspections),
            .groups = "drop") %>%
  pivot_longer(-carrier_operation, names_to = "BASIC", values_to = "rate") %>%
  ggplot(aes(BASIC, carrier_operation, fill = rate)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f", rate)), size = 3, color = "grey15") +
  scale_y_discrete(labels = OP_LABELS) +
  scale_fill_gradient(low = "white", high = "#7F1D1D") +
  labs(title = "Per-BASIC violation rate (violations / inspection) by operation type",
       subtitle = "These rates are exactly the features K-Means uses - patterns here preview the clusters",
       x = NULL, y = NULL, fill = "Rate") +
  theme_eda() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(c2d, "C2d_per_basic_rate_heatmap", h = 5)


# C3. GEOGRAPHY 
# C3a. Top 15 carrier states
c3a <- df %>% count(phy_state, sort = TRUE) %>% slice_head(n = 15) %>%
  ggplot(aes(reorder(phy_state, n), n)) +
  geom_col(fill = "#1F4E79") +
  geom_text(aes(label = comma(n)), hjust = -0.1, size = 3.3) +
  coord_flip() +
  scale_y_continuous(labels = comma, expand = expansion(c(0, 0.15))) +
  labs(title = "Top 15 states by carrier headquarters count",
       x = NULL, y = "Number of carriers") +
  theme_eda()
save_plot(c3a, "C3a_top_states", h = 6)

# C3b. State-level avg OOS rate (top 20 by carrier count)
c3b <- df %>%
  filter(has_inspection) %>%
  group_by(phy_state) %>%
  summarise(n = n(), avg_oos = mean(oos_rate, na.rm = TRUE),
            .groups = "drop") %>%
  filter(n >= 500) %>%
  slice_max(n, n = 20) %>%
  ggplot(aes(reorder(phy_state, avg_oos), avg_oos)) +
  geom_col(fill = "#DC2626") +
  geom_text(aes(label = percent(avg_oos, 0.1)),
            hjust = -0.1, size = 3) +
  coord_flip() +
  scale_y_continuous(labels = percent, expand = expansion(c(0, 0.15))) +
  labs(title = "Average carrier OOS rate by HQ state (top 20 by population)",
       subtitle = "Variance here is what the Location feature group captures in Q1",
       x = NULL, y = "Average OOS rate") +
  theme_eda()
save_plot(c3b, "C3b_state_oos_rate", h = 6)

# C3c. State-level total severity weight
c3c <- df %>%
  group_by(phy_state) %>%
  summarise(total_sev = sum(sum_severity_weight, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  slice_max(total_sev, n = 15) %>%
  ggplot(aes(reorder(phy_state, total_sev), total_sev)) +
  geom_col(fill = "#7F1D1D") +
  geom_text(aes(label = comma(round(total_sev))), hjust = -0.1, size = 3) +
  coord_flip() +
  scale_y_continuous(labels = comma, expand = expansion(c(0, 0.18))) +
  labs(title = "Top 15 HQ states by aggregate crash severity weight",
       x = NULL, y = "Sum of severity weight") +
  theme_eda()
save_plot(c3c, "C3c_state_severity", h = 6)


# C4. CLUSTER PREVIEW 

# C4a. Top 10% severity carriers vs the rest
c4a_data <- df %>%
  filter(has_inspection) %>%
  mutate(risk_decile = ntile(sum_severity_weight, 10),
         tier = if_else(risk_decile == 10, "Top 10% severity", "Bottom 90%")) %>%
  group_by(tier) %>%
  summarise(`Avg power units`        = mean(nbr_power_unit, na.rm = TRUE),
            `Avg drivers`            = mean(driver_total,   na.rm = TRUE),
            `Avg OOS rate`           = mean(oos_rate,       na.rm = TRUE),
            `Avg crashes`            = mean(total_crashes,  na.rm = TRUE),
            `Avg intensity score`    = mean(violation_intensity_score, na.rm = TRUE),
            `% Hazmat`               = mean(hm_flag_lgl)    * 100,
            `% Interstate`           = mean(carrier_operation == "A") * 100,
            .groups = "drop") %>%
  pivot_longer(-tier, names_to = "metric", values_to = "value")
c4a <- ggplot(c4a_data, aes(metric, value, fill = tier)) +
  geom_col(position = "dodge") +
  coord_flip() +
  scale_fill_manual(values = c("Bottom 90%" = "#9CA3AF",
                               "Top 10% severity" = "#DC2626")) +
  facet_wrap(~ metric, scales = "free", ncol = 2) +
  labs(title = "Top 10% severity carriers vs the rest",
       subtitle = "Differences across these axes are the seams K-Means will cut along",
       x = NULL, y = NULL, fill = NULL) +
  theme_eda() +
  theme(strip.text = element_text(size = 9),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top")
save_plot(c4a, "C4a_top_decile_profile", w = 10, h = 6)

# C4b. PCA preview of clustering features
clust_features <- df %>%
  filter(has_inspection) %>%
  transmute(nbr_power_unit_log = log1p(nbr_power_unit),
            driver_total_log   = log1p(driver_total),
            oos_rate           = oos_rate,
            crash_ratio        = crash_to_inspection_ratio,
            severity_log       = log1p(sum_severity_weight),
            fatalities_log     = log1p(total_fatalities),
            intensity_log      = log1p(violation_intensity_score),
            HM_Status          = HM_Status) %>%
  drop_na()
samp <- if (nrow(clust_features) > 30000)
  sample(seq_len(nrow(clust_features)), 30000) else seq_len(nrow(clust_features))
X      <- scale(clust_features[samp, 1:7])
pca    <- prcomp(X, center = FALSE, scale. = FALSE)
pca_df <- tibble(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                 HM_Status = clust_features$HM_Status[samp])
c4b <- ggplot(pca_df, aes(PC1, PC2, color = HM_Status)) +
  geom_point(alpha = 0.3, size = 0.8) +
  scale_color_manual(values = PAL_HM) +
  labs(title = "PCA projection of clustering features",
       subtitle = "If Hazmat carriers form a visible band/region, K-Means will likely separate them",
       x = sprintf("PC1 (%.1f%% var)", summary(pca)$importance[2, 1] * 100),
       y = sprintf("PC2 (%.1f%% var)", summary(pca)$importance[2, 2] * 100),
       color = NULL) +
  theme_eda() + theme(legend.position = "top")
save_plot(c4b, "C4b_pca_preview")


# =============================================================================
# Q1 — RANDOM FOREST: OUT-OF-SERVICE VIOLATION PREDICTION
# =============================================================================

# Dual-model comparison: Strict (leakage-free baseline) vs Behavioral Fingerprint.
# An early version reached AUC = 1.0 — a clear sign of target leakage / overfitting,
# which was diagnosed and removed. The two models below isolate genuine signal.

# WHY TWO MODELS?
# This script fits and compares two specifications of the same RF classifier:

#   MODEL A — STRICT (8 features):
#     Carrier profile + crash history only. NO inspection-derived features.
#     Tests the question: "Can carrier identity ALONE predict OOS?"
#     Expected AUC: ~0.62. Honest baseline.

#   MODEL B — BEHAVIORAL FINGERPRINT (17 features):
#     Strict features + violation rates + drivers/unit + carrier age + state tier.
#     Tests the question: "How much does historical behavior pattern improve
#     identification of OOS-flagged carriers?"
#     Expected AUC: ~0.80-0.88.

#   The GAP between A and B quantifies how much information the rate features
#   add. This comparison is the analytical contribution of the report.

# TRANSPARENCY CAVEAT (must appear in report):
#   Model B's rate features are computed from the same inspection records that
#   produce oos_flag (OOS violations are a subset of BASIC violations). Rates
#   therefore carry weak target echoes. We mitigate via (a) rates not raw
#   counts, (b) framing as descriptive/profiling not predictive, (c) reporting
#   alongside Model A as the leakage-free baseline.


N_CORES <- max(1L, parallel::detectCores() - 1L)


#Helpers & theme 


#AGAIN WRitten by AI since even M2 compute was taking time.. so told it to make it fase.
# Mann-Whitney AUC. as.numeric() casts BOTH operands to double so n_pos*n_neg
# (which can exceed 2.1 billion = R's integer cap on big test sets) does NOT
# overflow. Returns the same value as pROC::auc() but ~10-50x faster.

# ── STEP 1: BUILD ALL ENGINEERED FEATURES (used by both models) 
# We build the SUPERSET of features once, then each model picks its own
# subset. Single dataframe, single train/test split — fair comparison.


cat("Building engineered features...\n")

model_df <- df %>%
  filter(total_inspections >= 1) %>%
  mutate(
    # Carrier profile transforms (used by BOTH models)
    fleet_size_log    = log1p(pmax(nbr_power_unit, 0, na.rm = TRUE)),
    log_driver_total  = log1p(pmax(driver_total,   0, na.rm = TRUE)),
    
    # Drivers per power unit (Model B only) — staffing intensity proxy
    drivers_per_unit  = ifelse(nbr_power_unit > 0,
                               driver_total / nbr_power_unit, 0),
    
    # Behavioral rates (Model B only)
    rate_unsafe       = safe_div(viol_unsafe,       total_inspections),
    rate_fatigued     = safe_div(viol_fatigued,     total_inspections),
    rate_dr_fitness   = safe_div(viol_dr_fitness,   total_inspections),
    rate_subt_alc     = safe_div(viol_subt_alcohol, total_inspections),
    rate_vh_maint     = safe_div(viol_vh_maint,     total_inspections),
    rate_hm           = safe_div(viol_hm,           total_inspections),
    
    # Crash history (used by BOTH models)
    log_total_crashes   = log1p(pmax(total_crashes,       0, na.rm = TRUE)),
    log_severity_weight = log1p(pmax(sum_severity_weight, 0, na.rm = TRUE)),
    
    # Target
    oos_flag = factor(oos_flag, levels = c("Clean", "OOS"))
  ) %>%
  # Outlier handling (separate mutate so quantile sees populated column)
  mutate(
    drivers_per_unit  = pmin(drivers_per_unit,
                             quantile(drivers_per_unit, 0.99, na.rm = TRUE))) %>%
  filter(!is.na(phy_state), phy_state != "") %>%
  drop_na(oos_flag, fleet_size_log, log_driver_total, log_total_crashes,
          log_severity_weight, drivers_per_unit,
          hm_flag, is_interstate, is_intra_hazmat, is_intra_nonhaz,
          rate_unsafe, rate_fatigued, rate_dr_fitness, rate_subt_alc,
          rate_vh_maint, rate_hm, phy_state)

cat(sprintf("Modelling subset : %d carriers\n", nrow(model_df)))
cat(sprintf("  OOS   : %d  (%.1f%%)\n",
            sum(model_df$oos_flag == "OOS"),
            mean(model_df$oos_flag == "OOS") * 100))
cat(sprintf("  Clean : %d  (%.1f%%)\n",
            sum(model_df$oos_flag == "Clean"),
            mean(model_df$oos_flag == "Clean") * 100))


# ── STEP 2: SHARED 70/30 STRATIFIED SPLIT 
set.seed(42)
oos_idx   <- which(model_df$oos_flag == "OOS")
clean_idx <- which(model_df$oos_flag == "Clean")
train_idx <- c(sample(oos_idx,   floor(0.70 * length(oos_idx))),
               sample(clean_idx, floor(0.70 * length(clean_idx))))
test_idx  <- setdiff(seq_len(nrow(model_df)), train_idx)

train_data <- model_df[train_idx, ]
test_data  <- model_df[test_idx,  ]

cat(sprintf("\nTrain: %d  |  Test: %d\n", nrow(train_data), nrow(test_data)))


# ── STEP 3: STATE RISK TIER (Model B feature, train-encoded) 
state_oos_rate <- train_data %>%
  group_by(phy_state) %>%
  summarise(n = n(),
            state_oos_rate = mean(oos_flag == "OOS"),
            .groups = "drop") %>%
  filter(n >= 30)

tier_cuts <- quantile(state_oos_rate$state_oos_rate,
                      probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)

state_oos_rate <- state_oos_rate %>%
  mutate(state_risk_tier = cut(state_oos_rate,
                               breaks         = tier_cuts,
                               include.lowest = TRUE,
                               labels         = c("LOW", "MED", "HIGH"))) %>%
  select(phy_state, state_risk_tier)

attach_tier <- function(d) {
  d %>%
    left_join(state_oos_rate, by = "phy_state") %>%
    mutate(state_risk_tier = factor(
      ifelse(is.na(state_risk_tier), "MED", as.character(state_risk_tier)),
      levels = c("LOW", "MED", "HIGH")
    ))
}

train_data <- attach_tier(train_data)
test_data  <- attach_tier(test_data)

cat("\nState risk tier (train):\n")
print(table(train_data$state_risk_tier))


# ── STEP 4: CASE WEIGHTS (shared) 
n_train    <- nrow(train_data)
n_oos_tr   <- sum(train_data$oos_flag == "OOS")
n_clean_tr <- sum(train_data$oos_flag == "Clean")

w_oos   <- n_train / (2 * n_oos_tr)
w_clean <- n_train / (2 * n_clean_tr)
case_weights <- ifelse(train_data$oos_flag == "OOS", w_oos, w_clean)


# ── STEP 5: FEATURE SETS
# Model A — STRICT (8 features). Fully leakage-free.
FEATURES_A <- c(
  "fleet_size_log", "log_driver_total",
  "hm_flag", "is_interstate", "is_intra_hazmat", "is_intra_nonhaz",
  "log_total_crashes", "log_severity_weight"
)

# Model B — FINGERPRINT (17 features). Strict + 9 additional.
FEATURES_B <- c(
  FEATURES_A,
  "drivers_per_unit", "state_risk_tier",
  "rate_unsafe", "rate_fatigued", "rate_dr_fitness",
  "rate_subt_alc", "rate_vh_maint", "rate_hm"
)

cat(sprintf("\nModel A (Strict)      : %d features\n", length(FEATURES_A)))
cat(sprintf("Model B (Fingerprint) : %d features\n", length(FEATURES_B)))


# ── STEP 6: REUSABLE TRAIN-AND-EVALUATE FUNCTION 
# Takes a feature vector, returns a list with model object, predictions on
# train/test/CV, AUCs, and the optimal-threshold ROC. Used twice: once for
# Model A, once for Model B.

train_and_eval <- function(features, label, n_trees_final = 500) {
  cat(sprintf("\n─── %s (%d features) ───\n", label, length(features)))
  
  # CV tuning on 10% of train. Hyperparameter grid scaled to feature count.
  p <- length(features)
  tune_grid <- expand.grid(
    mtry          = unique(c(floor(sqrt(p)), floor(sqrt(p)) + 2,
                             max(2, floor(p / 3)))),
    min.node.size = c(20, 100, 300)
  ) %>% distinct()
  
  set.seed(42)
  tune_train <- train_data %>%
    group_by(oos_flag) %>%
    slice_sample(prop = 0.10) %>%
    ungroup()
  tune_cw <- ifelse(tune_train$oos_flag == "OOS", w_oos, w_clean)
  folds   <- caret::createFolds(tune_train$oos_flag, k = 5,
                                returnTrain = FALSE)
  
  cat(sprintf("  CV tuning on %d rows...\n", nrow(tune_train)))
  t0 <- Sys.time()
  
  cv_results <- map_dfr(seq_len(nrow(tune_grid)), function(i) {
    mt  <- tune_grid$mtry[i]
    mns <- tune_grid$min.node.size[i]
    fold_aucs <- map_dbl(folds, function(val_idx) {
      tr  <- tune_train[-val_idx, ]
      val <- tune_train[ val_idx, ]
      cw  <- tune_cw[-val_idx]
      fit <- ranger(
        formula       = oos_flag ~ .,
        data          = tr[, c(features, "oos_flag")],
        num.trees     = 300, mtry = mt, min.node.size = mns,
        case.weights  = cw, probability = TRUE,
        num.threads   = N_CORES, seed = 42
      )
      prob <- predict(fit, data = val[, features])$predictions[, "OOS"]
      fast_auc(val$oos_flag, prob, "OOS")
    })
    tibble(mtry = mt, min.node.size = mns,
           mean_auc = mean(fold_aucs), sd_auc = sd(fold_aucs))
  })
  
  cat(sprintf("  CV done in %.1fs\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  
  best <- cv_results %>% slice_max(mean_auc, n = 1)
  cat(sprintf("  Best: mtry=%d  min.node.size=%d  CV AUC=%.4f\n",
              best$mtry, best$min.node.size, best$mean_auc))
  
  # Final model on full training set
  cat("  Fitting final forest...\n")
  t0 <- Sys.time()
  set.seed(42)
  rf <- ranger(
    formula       = oos_flag ~ .,
    data          = train_data[, c(features, "oos_flag")],
    num.trees     = n_trees_final,
    mtry          = best$mtry,
    min.node.size = best$min.node.size,
    case.weights  = case_weights,
    probability   = TRUE,
    importance    = "permutation",
    num.threads   = N_CORES,
    seed          = 42
  )
  cat(sprintf("  Fit done in %.1fs   OOB error = %.4f\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs")),
              rf$prediction.error))
  
  # Predictions
  test_probs <- predict(rf, data = test_data[, features],
                        num.threads = N_CORES)$predictions[, "OOS"]
  
  # Train AUC on a 50k sample
  set.seed(42)
  ts_idx <- c(
    sample(which(train_data$oos_flag == "OOS"),
           min(25000, sum(train_data$oos_flag == "OOS"))),
    sample(which(train_data$oos_flag == "Clean"),
           min(25000, sum(train_data$oos_flag == "Clean")))
  )
  train_sample <- train_data[ts_idx, ]
  train_probs  <- predict(rf, data = train_sample[, features],
                          num.threads = N_CORES)$predictions[, "OOS"]
  
  auc_train <- fast_auc(train_sample$oos_flag, train_probs, "OOS")
  auc_cv    <- best$mean_auc
  auc_test  <- fast_auc(test_data$oos_flag,    test_probs,  "OOS")
  
  cat(sprintf("  AUC: train=%.4f  cv=%.4f  test=%.4f  (gap=%.4f)\n",
              auc_train, auc_cv, auc_test, auc_train - auc_test))
  
  # Test ROC and Youden threshold
  roc_test <- roc(test_data$oos_flag, test_probs,
                  levels = c("Clean", "OOS"), direction = "<", quiet = TRUE)
  yi <- which.max(roc_test$sensitivities + roc_test$specificities - 1)
  opt_thresh <- roc_test$thresholds[yi]
  opt_sens   <- roc_test$sensitivities[yi]
  opt_spec   <- roc_test$specificities[yi]
  
  pred_factor <- factor(ifelse(test_probs >= opt_thresh, "OOS", "Clean"),
                        levels = c("Clean", "OOS"))
  cm <- confusionMatrix(pred_factor, test_data$oos_flag, positive = "OOS")
  
  list(
    label       = label,
    features    = features,
    model       = rf,
    cv_results  = cv_results,
    best_params = best,
    test_probs  = test_probs,
    auc_train   = auc_train,
    auc_cv      = auc_cv,
    auc_test    = auc_test,
    roc_test    = roc_test,
    opt_thresh  = opt_thresh,
    opt_sens    = opt_sens,
    opt_spec    = opt_spec,
    cm          = cm
  )
}


# STEP 7: RUN BOTH MODELS 
res_A <- train_and_eval(FEATURES_A, "Model A — Strict")
res_B <- train_and_eval(FEATURES_B, "Model B — Fingerprint")


# STEP 8: SIDE-BY-SIDE COMPARISON TABLE
cat("\n\n══════ MODEL COMPARISON ══════\n")

comparison_tbl <- tibble(
  Metric = c("Features", "mtry", "min.node.size", "OOB Error",
             "Train AUC", "CV AUC", "Test AUC",
             "Train-Test Gap", "Sensitivity", "Specificity",
             "Balanced Accuracy", "Kappa"),
  `Strict (8)` = c(
    length(res_A$features),
    res_A$best_params$mtry,
    res_A$best_params$min.node.size,
    sprintf("%.4f", res_A$model$prediction.error),
    sprintf("%.4f", res_A$auc_train),
    sprintf("%.4f", res_A$auc_cv),
    sprintf("%.4f", res_A$auc_test),
    sprintf("%.4f", res_A$auc_train - res_A$auc_test),
    sprintf("%.4f", res_A$opt_sens),
    sprintf("%.4f", res_A$opt_spec),
    sprintf("%.4f", res_A$cm$byClass["Balanced Accuracy"]),
    sprintf("%.4f", res_A$cm$overall["Kappa"])
  ),
  `Fingerprint (17)` = c(
    length(res_B$features),
    res_B$best_params$mtry,
    res_B$best_params$min.node.size,
    sprintf("%.4f", res_B$model$prediction.error),
    sprintf("%.4f", res_B$auc_train),
    sprintf("%.4f", res_B$auc_cv),
    sprintf("%.4f", res_B$auc_test),
    sprintf("%.4f", res_B$auc_train - res_B$auc_test),
    sprintf("%.4f", res_B$opt_sens),
    sprintf("%.4f", res_B$opt_spec),
    sprintf("%.4f", res_B$cm$byClass["Balanced Accuracy"]),
    sprintf("%.4f", res_B$cm$overall["Kappa"])
  )
)
print(comparison_tbl, n = Inf)

cat(sprintf("\n→ AUC IMPROVEMENT (test) : %.4f → %.4f  (delta = +%.4f)\n",
            res_A$auc_test, res_B$auc_test,
            res_B$auc_test - res_A$auc_test))


# ── STEP 9: OVERLAID ROC CURVES 
roc_df <- bind_rows(
  tibble(FPR  = 1 - res_A$roc_test$specificities,
         TPR  = res_A$roc_test$sensitivities,
         Model = "Strict (8 feat)"),
  tibble(FPR  = 1 - res_B$roc_test$specificities,
         TPR  = res_B$roc_test$sensitivities,
         Model = "Fingerprint (17 feat)")
)

p_roc_compare <- ggplot(roc_df, aes(FPR, TPR, color = Model)) +
  geom_line(linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey50") +
  annotate("point",
           x = 1 - res_A$opt_spec, y = res_A$opt_sens,
           color = "#F97316", size = 3.5, shape = 18) +
  annotate("point",
           x = 1 - res_B$opt_spec, y = res_B$opt_sens,
           color = "#1F4E79", size = 3.5, shape = 18) +
  annotate("text", x = 0.6, y = 0.18,
           label = sprintf("Strict      AUC = %.4f", res_A$auc_test),
           size = 4, fontface = "bold", color = "#F97316", hjust = 0) +
  annotate("text", x = 0.6, y = 0.10,
           label = sprintf("Fingerprint AUC = %.4f", res_B$auc_test),
           size = 4, fontface = "bold", color = "#1F4E79", hjust = 0) +
  scale_color_manual(values = PAL_MODELS) +
  labs(
    title    = "Q1 — ROC Comparison: Strict vs Fingerprint",
    subtitle = "Diamond = Youden's J optimal threshold per model",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = NULL
  ) +
  theme_eda() +
  theme(legend.position = "top")

save_plot(p_roc_compare, "Q1_RF_roc_comparison", w = 9, h = 6)
print(p_roc_compare)


#STEP 10: AUC BAR CHART (Train / CV / Test × 2 models) 
auc_long <- tibble(
  Model = factor(rep(c("Strict (8 feat)", "Fingerprint (17 feat)"), each = 3),
                 levels = c("Strict (8 feat)", "Fingerprint (17 feat)")),
  Split = factor(rep(c("Train (50k)", "5-fold CV", "Test"), 2),
                 levels = c("Train (50k)", "5-fold CV", "Test")),
  AUC = c(res_A$auc_train, res_A$auc_cv, res_A$auc_test,
          res_B$auc_train, res_B$auc_cv, res_B$auc_test)
)

p_auc_bar <- ggplot(auc_long, aes(Split, AUC, fill = Model)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", AUC)),
            position = position_dodge(width = 0.8),
            vjust = -0.4, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = PAL_MODELS) +
  scale_y_continuous(limits = c(0, 1.0),
                     labels = number_format(accuracy = 0.01)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  annotate("text", x = 3, y = 0.53,
           label = "Random baseline (AUC = 0.5)",
           size = 3, color = "grey40", hjust = 1) +
  labs(
    title    = "Q1 — AUC Comparison Across Splits and Models",
    subtitle = sprintf("Test AUC gain: %.4f → %.4f  (Δ = +%.4f)",
                       res_A$auc_test, res_B$auc_test,
                       res_B$auc_test - res_A$auc_test),
    x = NULL, y = "AUC", fill = NULL
  ) +
  theme_eda() +
  theme(legend.position = "top")

save_plot(p_auc_bar, "Q1_RF_auc_comparison_bar", w = 10, h = 6)
print(p_auc_bar)


# STEP 11: SIDE-BY-SIDE IMPORTANCE PLOTS 
recoder <- c(
  fleet_size_log      = "Fleet Size (log)",
  log_driver_total    = "Driver Total (log)",
  drivers_per_unit    = "Drivers per Power Unit",
  hm_flag             = "Hazmat Flag",
  is_interstate       = "Interstate Op.",
  is_intra_hazmat     = "Intrastate Hazmat Op.",
  is_intra_nonhaz     = "Intrastate Non-Haz Op.",
  state_risk_tier     = "State Risk Tier",
  rate_unsafe         = "Unsafe Driving Rate",
  rate_fatigued       = "Fatigued Driving Rate",
  rate_dr_fitness     = "Driver Fitness Rate",
  rate_subt_alc       = "Substance/Alcohol Rate",
  rate_vh_maint       = "Vehicle Maintenance Rate",
  rate_hm             = "Hazmat Compliance Rate",
  log_total_crashes   = "Total Crashes (log)",
  log_severity_weight = "Crash Severity Weight (log)"
)

build_imp <- function(res) {
  tibble(feature    = names(res$model$variable.importance),
         importance = res$model$variable.importance) %>%
    mutate(feature = recode(feature, !!!recoder)) %>%
    arrange(desc(importance))
}

imp_A <- build_imp(res_A)
imp_B <- build_imp(res_B)

p_imp_A <- ggplot(imp_A, aes(importance, reorder(feature, importance))) +
  geom_col(fill = "#F97316") +
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.4) +
  labs(title = "Strict (8 features)",
       x = "Permutation Importance", y = NULL) +
  theme_eda()

p_imp_B <- ggplot(imp_B, aes(importance, reorder(feature, importance))) +
  geom_col(fill = "#1F4E79") +
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.4) +
  labs(title = "Fingerprint (17 features)",
       x = "Permutation Importance", y = NULL) +
  theme_eda()

p_imp_compare <- p_imp_A + p_imp_B +
  plot_annotation(
    title    = "Q1 — Variable Importance: Strict vs Fingerprint",
    subtitle = "Note how rate features dominate Model B's importance ranking",
    theme    = theme(plot.title = element_text(face = "bold", size = 13))
  )

save_plot(p_imp_compare, "Q1_RF_importance_comparison", w = 13, h = 6)
print(p_imp_compare)


# STEP 12: PROBABILITY DENSITY OVERLAY 
prob_df <- bind_rows(
  tibble(prob = res_A$test_probs, oos_flag = test_data$oos_flag,
         Model = "Strict (8 feat)"),
  tibble(prob = res_B$test_probs, oos_flag = test_data$oos_flag,
         Model = "Fingerprint (17 feat)")
) %>%
  mutate(Model = factor(Model,
                        levels = c("Strict (8 feat)", "Fingerprint (17 feat)")))

p_prob <- ggplot(prob_df, aes(prob, fill = oos_flag, color = oos_flag)) +
  geom_density(alpha = 0.4, linewidth = 0.7) +
  facet_wrap(~ Model, ncol = 2) +
  scale_fill_manual(values  = PAL_OOS) +
  scale_color_manual(values = PAL_OOS) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title    = "Q1 — Predicted P(OOS) Density by Class",
    subtitle = "Wider separation (right panel) = stronger discrimination",
    x = "Predicted P(OOS)", y = "Density",
    fill = NULL, color = NULL
  ) +
  theme_eda() +
  theme(legend.position = "top")

save_plot(p_prob, "Q1_RF_probability_comparison", w = 12, h = 5)
print(p_prob)


# STEP 13: OOB CONVERGENCE FOR MODEL B (the headline model) 
cat("\nComputing OOB convergence for Model B...\n")
t0 <- Sys.time()
set.seed(42)
rf_conv <- ranger(
  formula       = oos_flag ~ .,
  data          = train_data[, c(FEATURES_B, "oos_flag")],
  num.trees     = 500,
  mtry          = res_B$best_params$mtry,
  min.node.size = res_B$best_params$min.node.size,
  case.weights  = case_weights,
  probability   = FALSE,
  keep.inbag    = TRUE,
  num.threads   = N_CORES,
  seed          = 42
)
tree_preds <- predict(rf_conv, data = train_data[, FEATURES_B],
                      predict.all = TRUE,
                      num.threads = N_CORES)$predictions
inbag    <- do.call(cbind, rf_conv$inbag.counts)
oob_mask <- inbag == 0L
y_int    <- as.integer(train_data$oos_flag)

tree_seq <- c(10, 25, 50, 100, 150, 200, 300, 400, 500)
oob_seq <- map_dbl(tree_seq, function(nt) {
  preds_sub <- tree_preds[, seq_len(nt), drop = FALSE]
  mask_sub  <- oob_mask[,  seq_len(nt), drop = FALSE]
  oos_votes <- rowSums((preds_sub == 2L) & mask_sub)
  total     <- rowSums(mask_sub)
  has_vote  <- total > 0
  pred      <- ifelse(oos_votes / pmax(total, 1) >= 0.5, 2L, 1L)
  mean(pred[has_vote] != y_int[has_vote])
})
cat(sprintf("OOB convergence done in %.1fs\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

oob_df <- tibble(num_trees = tree_seq, oob_error = oob_seq)
p_oob <- ggplot(oob_df, aes(num_trees, oob_error)) +
  geom_line(color = "#1F4E79", linewidth = 1.2) +
  geom_point(size = 3, color = "#1F4E79") +
  geom_hline(yintercept = res_B$model$prediction.error,
             linetype = "dashed", color = "#DC2626") +
  scale_x_continuous(breaks = tree_seq) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(title    = "Q1 (Model B) — OOB Error Convergence",
       subtitle = "Flat plateau = stable; no overfitting as trees grow",
       x = "Number of Trees", y = "OOB Classification Error") +
  theme_eda()
save_plot(p_oob, "Q1_RF_oob_convergence", w = 9, h = 5)
print(p_oob)


# Q1 SUMMARY 
cat("\nMODEL A — STRICT (8 features, leakage-free)\n")
cat(sprintf("  Test AUC = %.4f  |  Bal.Acc = %.4f  |  Kappa = %.4f\n",
            res_A$auc_test,
            res_A$cm$byClass["Balanced Accuracy"],
            res_A$cm$overall["Kappa"]))

cat("\nMODEL B — FINGERPRINT (17 features, includes rates)\n")
cat(sprintf("  Test AUC = %.4f  |  Bal.Acc = %.4f  |  Kappa = %.4f\n",
            res_B$auc_test,
            res_B$cm$byClass["Balanced Accuracy"],
            res_B$cm$overall["Kappa"]))

cat(sprintf("\nINFORMATION GAIN FROM RATE FEATURES: AUC +%.4f\n",
            res_B$auc_test - res_A$auc_test))

cat("\n✔  Q1 complete.\n")


# =============================================================================
# QUESTION 2 — K-MEANS UNSUPERVISED CLUSTERING
# Segment carriers by fleet attributes + violation behaviour + crash history
# =============================================================================


cat("\n── Q2: K-MEANS CLUSTERING ──\n")


# ── STEP 1: Build clustering feature matrix 
# Features exactly match Q2's research question:
#   Fleet/Ops    : fleet_size_log, is_interstate, is_intra_hazmat,
#                  is_intra_nonhaz, hm_flag
#   Violation    : per-BASIC OOS rates (violations / total_inspections)
#   Crash hist.  : total_crashes (log), sum_severity_weight (log),
#                  total_fatalities (log), total_injuries (log)

cluster_base <- df %>%
  filter(total_inspections >= 1) %>%     # need inspections for BASIC rates
  transmute(
    dot_number,
    # Fleet & Operations
    fleet_size_log   = log1p(pmax(nbr_power_unit, 0, na.rm = TRUE)),
    is_interstate    = as.numeric(is_interstate),
    is_intra_hazmat  = as.numeric(is_intra_hazmat),
    is_intra_nonhaz  = as.numeric(is_intra_nonhaz),
    hm_flag          = as.numeric(hm_flag),
    # Per-BASIC OOS rates (violations ÷ inspections — bounded, comparable)
    rate_unsafe      = safe_div(viol_unsafe,       total_inspections),
    rate_fatigued    = safe_div(viol_fatigued,      total_inspections),
    rate_dr_fitness  = safe_div(viol_dr_fitness,    total_inspections),
    rate_subt_alc    = safe_div(viol_subt_alcohol,  total_inspections),
    rate_vh_maint    = safe_div(viol_vh_maint,      total_inspections),
    rate_hm          = safe_div(viol_hm,            total_inspections),
    # Crash history (log-transform to reduce extreme skew)
    crashes_log      = log1p(total_crashes),
    severity_log     = log1p(sum_severity_weight),
    fatalities_log   = log1p(total_fatalities),
    injuries_log     = log1p(total_injuries),
    # Raw labels — NOT fed to algorithm, only for post-cluster profiling
    oos_rate_raw         = oos_rate,
    total_crashes_raw    = total_crashes,
    severity_weight_raw  = sum_severity_weight,
    total_fatalities_raw = total_fatalities,
    total_injuries_raw   = total_injuries,
    carrier_operation,
    hm_flag_raw          = hm_flag
  ) %>%
  drop_na()

cat(sprintf("Clustering subset: %d carriers\n", nrow(cluster_base)))

# Separate feature matrix from label columns
FEATURE_COLS <- c("fleet_size_log", "is_interstate", "is_intra_hazmat",
                  "is_intra_nonhaz", "hm_flag",
                  "rate_unsafe", "rate_fatigued", "rate_dr_fitness",
                  "rate_subt_alc", "rate_vh_maint", "rate_hm",
                  "crashes_log", "severity_log", "fatalities_log", "injuries_log")

X_raw    <- cluster_base[, FEATURE_COLS]
X_scaled <- scale(X_raw)    # z-score normalisation (mandatory for K-Means)


# ── STEP 2: Elbow Method (Within-Cluster Sum of Squares) 
K_MAX <- 10
set.seed(42)

wss_vals <- map_dbl(1:K_MAX, function(k) {
  kmeans(X_scaled, centers = k, nstart = 25, iter.max = 300)$tot.withinss
})

elbow_df <- tibble(k = 1:K_MAX, wss = wss_vals)

p_elbow <- ggplot(elbow_df, aes(k, wss)) +
  geom_line(color = "#1F4E79", linewidth = 1.1) +
  geom_point(size = 3.5, color = "#1F4E79") +
  geom_point(data = elbow_df %>% filter(k == 4),   # annotate the elbow
             aes(k, wss), color = "#DC2626", size = 5, shape = 18) +
  scale_x_continuous(breaks = 1:K_MAX) +
  scale_y_continuous(labels = comma) +
  labs(
    title    = "Q2 K-Means — Elbow Method (Total Within-Cluster SS)",
    subtitle = "Look for the 'elbow' — the k where adding more clusters yields diminishing returns",
    x        = "Number of Clusters (k)",
    y        = "Total Within-Cluster Sum of Squares"
  ) +
  theme_eda()

save_plot(p_elbow, "Q2_elbow_method")
print(p_elbow)
cat(sprintf("\nWSS values:\n"))
print(elbow_df)


# ── STEP 3: Silhouette Score 
# Silhouette ∈ [-1, 1]; higher = better-defined, well-separated clusters.
# Computed on a sample for speed (silhouette is O(n²) in memory).

set.seed(42)
sil_n   <- min(nrow(X_scaled), 20000)   # cap at 20k for speed
sil_idx <- sample(seq_len(nrow(X_scaled)), sil_n)
X_sil   <- X_scaled[sil_idx, ]

sil_scores <- map_dbl(2:K_MAX, function(k) {
  km  <- kmeans(X_sil, centers = k, nstart = 25, iter.max = 300)
  sil <- silhouette(km$cluster, dist(X_sil))
  mean(sil[, "sil_width"])
})

sil_df <- tibble(k = 2:K_MAX, avg_silhouette = sil_scores)
best_k_sil   <- sil_df$k[which.max(sil_df$avg_silhouette)]
best_sil_val <- max(sil_df$avg_silhouette)

cat(sprintf("\nSilhouette scores by k:\n"))
print(sil_df)
cat(sprintf("\n▶  Best k by Silhouette = %d  (avg silhouette = %.4f)\n",
            best_k_sil, best_sil_val))

p_silhouette <- ggplot(sil_df, aes(k, avg_silhouette)) +
  geom_line(color = "#7F1D1D", linewidth = 1.1) +
  geom_point(size = 3.5, color = "#7F1D1D") +
  geom_point(data = sil_df %>% filter(k == best_k_sil),
             aes(k, avg_silhouette),
             color = "#16A34A", size = 5, shape = 18) +
  annotate("text",
           x     = best_k_sil + 0.15,
           y     = best_sil_val,
           label = sprintf("Best k = %d\n(%.3f)", best_k_sil, best_sil_val),
           hjust = 0, size = 3.5, color = "#16A34A") +
  scale_x_continuous(breaks = 2:K_MAX) +
  labs(
    title    = "Q2 K-Means — Average Silhouette Score by k",
    subtitle = "Higher = clusters are more cohesive and well-separated; green diamond = best k",
    x        = "Number of Clusters (k)",
    y        = "Average Silhouette Width"
  ) +
  theme_eda()

save_plot(p_silhouette, "Q2_silhouette_score")
print(p_silhouette)


# ── STEP 4: Side-by-side Elbow + Silhouette (combined figure)
p_combined_k <- p_elbow + p_silhouette +
  plot_annotation(
    title    = "Determining Optimal k — Elbow Method vs Silhouette Score",
    subtitle = "Both methods should agree, or pick the simpler/more interpretable k",
    theme    = theme(plot.title = element_text(face = "bold", size = 13))
  )

save_plot(p_combined_k, "Q2_k_selection_combined", w = 14, h = 5)
print(p_combined_k)


# ── STEP 5: Fit Final K-Means Model 
# Set K_FINAL based on the elbow / silhouette evidence above.

K_FINAL <- 4   # ← adjust after inspecting the elbow and silhouette plots

set.seed(42)
km_final <- kmeans(X_scaled, centers = K_FINAL,
                   nstart = 50, iter.max = 500)

# Attach cluster labels back to the carrier data
cluster_base <- cluster_base %>%
  mutate(cluster = factor(km_final$cluster,
                          labels = paste0("C", 1:K_FINAL)))

cat(sprintf("\nFinal model: K = %d\n", K_FINAL))
cat("\nCluster sizes:\n")
print(table(cluster_base$cluster))
cat(sprintf("Between-SS / Total-SS = %.2f%%\n",
            km_final$betweenss / km_final$totss * 100))


# ── STEP 6: Cluster Profiling Table 
profile_tbl <- cluster_base %>%
  group_by(cluster) %>%
  summarise(
    n_carriers       = n(),
    # Fleet
    med_fleet_size   = median(exp(fleet_size_log) - 1, na.rm = TRUE),
    pct_interstate   = mean(is_interstate) * 100,
    pct_hazmat       = mean(hm_flag_raw)   * 100,
    # Violation behaviour
    avg_oos_rate     = mean(oos_rate_raw,        na.rm = TRUE),
    avg_rate_unsafe  = mean(rate_unsafe,         na.rm = TRUE),
    avg_rate_fatigue = mean(rate_fatigued,        na.rm = TRUE),
    avg_rate_vhmaint = mean(rate_vh_maint,        na.rm = TRUE),
    avg_rate_hm      = mean(rate_hm,              na.rm = TRUE),
    # Crash history
    avg_crashes      = mean(total_crashes_raw,    na.rm = TRUE),
    avg_severity_wt  = mean(severity_weight_raw,  na.rm = TRUE),
    avg_fatalities   = mean(total_fatalities_raw, na.rm = TRUE),
    avg_injuries     = mean(total_injuries_raw,   na.rm = TRUE),
    .groups = "drop"
  )

cat("\n── Cluster Profile Table \n")
print(profile_tbl, n = K_FINAL, max.cols = 100)


# ── STEP 7: Statistical Tests — Do clusters differ significantly? ─────────────
# Kruskal-Wallis (non-parametric ANOVA) — appropriate given the skewed,
# count-heavy distributions.

kw_oos      <- kruskal.test(oos_rate_raw      ~ cluster, data = cluster_base)
kw_crashes  <- kruskal.test(total_crashes_raw ~ cluster, data = cluster_base)
kw_severity <- kruskal.test(severity_weight_raw ~ cluster, data = cluster_base)

cat("\n── Kruskal-Wallis Tests (H₀: all clusters have equal medians) ──────────\n")
cat(sprintf("OOS Rate       — chi² = %.2f, df = %d, p = %.2e  %s\n",
            kw_oos$statistic, kw_oos$parameter, kw_oos$p.value,
            ifelse(kw_oos$p.value < 0.05, "✔ SIGNIFICANT", "✘ not significant")))
cat(sprintf("Crash Frequency— chi² = %.2f, df = %d, p = %.2e  %s\n",
            kw_crashes$statistic, kw_crashes$parameter, kw_crashes$p.value,
            ifelse(kw_crashes$p.value < 0.05, "✔ SIGNIFICANT", "✘ not significant")))
cat(sprintf("Severity Weight— chi² = %.2f, df = %d, p = %.2e  %s\n",
            kw_severity$statistic, kw_severity$parameter, kw_severity$p.value,
            ifelse(kw_severity$p.value < 0.05, "✔ SIGNIFICANT", "✘ not significant")))


# ── STEP 8: PLOTS 
CLUSTER_COLORS <- c("C1" = "#1f77b4", "C2" = "#d62728",
                    "C3" = "#2ca02c", "C4" = "#ff7f0e",
                    "C5" = "#9467bd", "C6" = "#8c564b")[1:K_FINAL]


# ── 8a. PCA 2D Cluster Projection 
set.seed(42)
pca_n   <- min(nrow(X_scaled), 30000)
pca_idx <- sample(seq_len(nrow(X_scaled)), pca_n)
X_pca   <- X_scaled[pca_idx, ]
pca_fit <- prcomp(X_pca, center = FALSE, scale. = FALSE)

pca_plot_df <- tibble(
  PC1     = pca_fit$x[, 1],
  PC2     = pca_fit$x[, 2],
  cluster = cluster_base$cluster[pca_idx]
)

pct_var <- summary(pca_fit)$importance[2, 1:2] * 100

p_pca <- ggplot(pca_plot_df, aes(PC1, PC2, color = cluster)) +
  geom_point(alpha = 0.25, size = 0.7) +
  stat_ellipse(linewidth = 1.0, level = 0.90) +
  scale_color_manual(values = CLUSTER_COLORS) +
  labs(
    title    = "Q2 — K-Means Clusters in PCA Space (PC1 vs PC2)",
    subtitle = sprintf("Ellipses = 90%% confidence regions  |  PC1: %.1f%% var  |  PC2: %.1f%% var",
                       pct_var[1], pct_var[2]),
    x        = sprintf("PC1 (%.1f%% variance)", pct_var[1]),
    y        = sprintf("PC2 (%.1f%% variance)", pct_var[2]),
    color    = "Cluster"
  ) +
  theme_eda() +
  theme(legend.position = "right")

save_plot(p_pca, "Q2_pca_cluster_projection", w = 9, h = 6)
print(p_pca)


# ── 8b. OOS Rate by Cluster (Boxplot) 
p_oos_box <- ggplot(cluster_base,
                    aes(cluster, oos_rate_raw, fill = cluster)) +
  geom_boxplot(outlier.alpha = 0.15, show.legend = FALSE) +
  scale_fill_manual(values = CLUSTER_COLORS) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, quantile(cluster_base$oos_rate_raw,
                                            0.99, na.rm = TRUE))) +
  labs(
    title    = "Q2 — Out-of-Service Rate by Cluster",
    subtitle = sprintf("Kruskal-Wallis p = %.2e — clusters differ significantly in OOS rate",
                       kw_oos$p.value),
    x        = "Cluster",
    y        = "Out-of-Service Rate"
  ) +
  theme_eda()

save_plot(p_oos_box, "Q2_oos_rate_by_cluster")
print(p_oos_box)


# ── 8c. Total Crashes by Cluster (Boxplot) 
p_crash_box <- ggplot(
  cluster_base %>% filter(total_crashes_raw > 0),
  aes(cluster, total_crashes_raw, fill = cluster)) +
  geom_boxplot(outlier.alpha = 0.15, show.legend = FALSE) +
  scale_y_log10(labels = comma) +
  scale_fill_manual(values = CLUSTER_COLORS) +
  labs(
    title    = "Q2 — Total Crashes by Cluster (log10 scale, carriers with ≥1 crash)",
    subtitle = sprintf("Kruskal-Wallis p = %.2e", kw_crashes$p.value),
    x        = "Cluster",
    y        = "Total Crashes (log10)"
  ) +
  theme_eda()

save_plot(p_crash_box, "Q2_crashes_by_cluster")
print(p_crash_box)


# ── 8d. Aggregate Severity Weight by Cluster 
p_severity_box <- ggplot(
  cluster_base %>% filter(severity_weight_raw > 0),
  aes(cluster, severity_weight_raw, fill = cluster)) +
  geom_boxplot(outlier.alpha = 0.15, show.legend = FALSE) +
  scale_y_log10(labels = comma) +
  scale_fill_manual(values = CLUSTER_COLORS) +
  labs(
    title    = "Q2 — Aggregate Crash Severity Weight by Cluster (log10, non-zero only)",
    subtitle = sprintf("Kruskal-Wallis p = %.2e", kw_severity$p.value),
    x        = "Cluster",
    y        = "Sum of Severity Weight (log10)"
  ) +
  theme_eda()

save_plot(p_severity_box, "Q2_severity_by_cluster")
print(p_severity_box)


# ── 8e. Cluster Mean Profile Heatmap 
# Standardised cluster centroids → shows which dimension drives each cluster.

centroid_df <- as.data.frame(km_final$centers) %>%
  rownames_to_column("cluster_num") %>%
  mutate(cluster = paste0("C", cluster_num)) %>%
  select(-cluster_num) %>%
  pivot_longer(-cluster, names_to = "feature", values_to = "z_score") %>%
  mutate(feature = recode(feature,
                          fleet_size_log  = "Fleet Size (log)",
                          is_interstate   = "Interstate Op.",
                          is_intra_hazmat = "Intrastate Hazmat Op.",
                          is_intra_nonhaz = "Intrastate Non-Haz Op.",
                          hm_flag         = "Hazmat Flag",
                          rate_unsafe     = "Unsafe Driving Rate",
                          rate_fatigued   = "Fatigued Driving Rate",
                          rate_dr_fitness = "Driver Fitness Rate",
                          rate_subt_alc   = "Substance/Alcohol Rate",
                          rate_vh_maint   = "Vehicle Maint. Rate",
                          rate_hm         = "Hazmat Compliance Rate",
                          crashes_log     = "Total Crashes (log)",
                          severity_log    = "Severity Weight (log)",
                          fatalities_log  = "Total Fatalities (log)",
                          injuries_log    = "Total Injuries (log)"
  ))

p_heatmap <- ggplot(centroid_df, aes(cluster, feature, fill = z_score)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", z_score)),
            size = 3, color = "grey10") +
  scale_fill_gradient2(low = "#1F4E79", mid = "white", high = "#7F1D1D",
                       midpoint = 0, name = "Z-score") +
  labs(
    title    = "Q2 — Cluster Centroids Heatmap (Standardised Z-Scores)",
    subtitle = "Positive (red) = above average; Negative (blue) = below average for that feature",
    x        = "Cluster",
    y        = NULL
  ) +
  theme_eda() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 9))

save_plot(p_heatmap, "Q2_cluster_centroids_heatmap", w = 8, h = 7)
print(p_heatmap)


# ── 8f. Carrier Operation Mix by Cluster (Stacked Bar) 
p_op_mix <- cluster_base %>%
  count(cluster, carrier_operation) %>%
  group_by(cluster) %>%
  mutate(pct = n / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(cluster, pct, fill = carrier_operation)) +
  geom_col(position = "fill") +
  geom_text(aes(label = percent(pct, 1)),
            position = position_fill(vjust = 0.5), size = 3.5) +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(
    values = c("A" = "#1f77b4", "B" = "#d62728", "C" = "#7f7f7f"),
    labels = c("A" = "Interstate", "B" = "Intrastate Hazmat",
               "C" = "Intrastate Non-Hazmat")
  ) +
  labs(
    title    = "Q2 — Carrier Operation Mix within Each Cluster",
    subtitle = "If K-Means is capturing operation type, segments will have distinct operation compositions",
    x        = "Cluster",
    y        = "Share of carriers",
    fill     = "Operation Type"
  ) +
  theme_eda() +
  theme(legend.position = "top")

save_plot(p_op_mix, "Q2_operation_mix_by_cluster")
print(p_op_mix)


# ── 8g. Combined outcome comparison (3-panel)
g1 <- ggplot(cluster_base,
             aes(cluster, oos_rate_raw, fill = cluster)) +
  geom_violin(show.legend = FALSE, alpha = 0.7, trim = TRUE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA,
               show.legend = FALSE) +
  scale_fill_manual(values = CLUSTER_COLORS) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(title = "OOS Rate", x = NULL, y = NULL) +
  theme_eda()

g2 <- ggplot(
  cluster_base %>% filter(total_crashes_raw > 0),
  aes(cluster, total_crashes_raw, fill = cluster)) +
  geom_violin(show.legend = FALSE, alpha = 0.7, trim = TRUE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA,
               show.legend = FALSE) +
  scale_fill_manual(values = CLUSTER_COLORS) +
  scale_y_log10(labels = comma) +
  labs(title = "Total Crashes (log)", x = NULL, y = NULL) +
  theme_eda()

g3 <- ggplot(
  cluster_base %>% filter(severity_weight_raw > 0),
  aes(cluster, severity_weight_raw, fill = cluster)) +
  geom_violin(show.legend = FALSE, alpha = 0.7, trim = TRUE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA,
               show.legend = FALSE) +
  scale_fill_manual(values = CLUSTER_COLORS) +
  scale_y_log10(labels = comma) +
  labs(title = "Severity Weight (log)", x = NULL, y = NULL) +
  theme_eda()

p_combined_outcome <- g1 + g2 + g3 +
  plot_annotation(
    title    = "Q2 — Key Outcome Metrics Across Clusters (Violin + Boxplot)",
    subtitle = "Three axes Q2 specifically asks about: OOS rate | crash frequency | severity weight",
    theme    = theme(plot.title = element_text(face = "bold", size = 13))
  )

save_plot(p_combined_outcome, "Q2_combined_outcomes", w = 14, h = 5)
print(p_combined_outcome)


# ── STEP 9: Save cluster assignments 
clustered_output <- df %>%
  filter(total_inspections >= 1) %>%
  semi_join(cluster_base %>% select(dot_number), by = "dot_number") %>%
  left_join(cluster_base %>% select(dot_number, cluster), by = "dot_number")

write_csv(clustered_output, "Final_file_with_clusters.csv")
cat(sprintf("\nCluster assignments written to: Final_file_with_clusters.csv\n"))


# SUMMARY
cat("  EXECUTION COMPLETE — SUMMARY\n")

cat("\nQUESTION 2 — K-Means Clustering:\n")
cat(sprintf("  Final k                = %d\n", K_FINAL))
cat(sprintf("  Between-SS / Total-SS  = %.1f%%\n",
            km_final$betweenss / km_final$totss * 100))
cat(sprintf("  Best Silhouette k      = %d (score = %.4f)\n",
            best_k_sil, best_sil_val))
cat(sprintf("  Kruskal-Wallis (OOS)   p = %.2e\n", kw_oos$p.value))
cat(sprintf("  Kruskal-Wallis (Crash) p = %.2e\n", kw_crashes$p.value))
cat(sprintf("  Kruskal-Wallis (Sev)   p = %.2e\n", kw_severity$p.value))

cat(sprintf("\nAll figures saved to: %s/\n", FIG_DIR))
