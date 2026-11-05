#!/usr/bin/env Rscript




suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# Enrichment comparators
#
# Complementary landmark analyses.
#
# Analyses:
#   1) Maximum KDIGO stage 3 reached by each landmark (H24/H48),
#      using ALL available SCr measurements from ICU admission to landmark,
#      with RRT by landmark automatically counted as KDIGO stage 3.
#   2) M2 selection matched to the exact fraction selected by RRT
#      at each landmark.
#   3) M2 top 10% / top 20% enrichment among patients NOT on RRT
#      at each landmark.
#
# Outcome:
#   Subsequent ICU death after the landmark.
#
# Bootstrap:
#   2,000 ICU-stay-level bootstrap resamples.
# ============================================================

# ------------------------------------------------------------
# 0. SETTINGS
# ------------------------------------------------------------

BASE_DIR <- Sys.getenv("AKI_BASE_DIR")
if (!nzchar(BASE_DIR)) stop("Set AKI_BASE_DIR to the project data directory")

MIMIC_BF_RDS <- file.path(
  BASE_DIR,
  "DF/MIMIC/mimic31_aki_admission_base_v3.rds"
)

PRED_RDS <- file.path(
  BASE_DIR,
  "XGB_external_AUMC_train__MIMIC_G3_LANDMARK_RISKSET_v4",
  "LANDMARK_pred_dt.rds"
)

PRED_CSV <- file.path(
  BASE_DIR,
  "XGB_external_AUMC_train__MIMIC_G3_LANDMARK_RISKSET_v4",
  "LANDMARK_pred_dt.csv"
)

TS_CANDIDATES <- c(
  file.path(BASE_DIR, "DF/MIMIC/timeseries_0_14d_core_long.rds"),
  file.path(BASE_DIR, "DF/MIMIC/timeseries_mimic.rds"),
  file.path(BASE_DIR, "DF/MIMIC/MIMIC_datM_ts_AKI_SCr_h14.rds")
)

OUT_DIR <- file.path(BASE_DIR, "Enrichment_analysis", "Enrichment_comparators")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

MODEL_USE <- "M2"
TARGET_PROBABILITY <- "pC2"

LANDMARK_MAP <- data.table(
  manuscript_day = c("Day 1", "Day 2"),
  landmark_h = c(24L, 48L),
  pred_day_code = c(0L, 1L),
  observation_window = c("0–24 h", "0–48 h")
)

N_BOOT <- 2000L
BOOT_SEED <- 20260728L
SAVE_BOOT_REPLICATES <- FALSE

# ------------------------------------------------------------
# 1. HELPERS
# ------------------------------------------------------------

to_posix <- function(x) {
  if (inherits(x, "POSIXt")) return(as.POSIXct(x, tz = "UTC"))
  if (is.numeric(x)) {
    return(as.POSIXct(x, origin = "1970-01-01", tz = "UTC"))
  }
  as.POSIXct(x, tz = "UTC")
}

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

safe_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  as.numeric(
    quantile(
      x,
      probs = probability,
      na.rm = TRUE,
      names = FALSE,
      type = 7
    )
  )
}

fmt_pct_ci <- function(estimate, low, high, digits = 1) {
  if (any(is.na(c(estimate, low, high)))) return(NA_character_)
  sprintf(
    paste0("%.", digits, "f%% [%.", digits, "f–%.", digits, "f]"),
    100 * estimate,
    100 * low,
    100 * high
  )
}

fmt_num_ci <- function(estimate, low, high, digits = 2) {
  if (any(is.na(c(estimate, low, high)))) return(NA_character_)
  sprintf(
    paste0("%.", digits, "f [%.", digits, "f–%.", digits, "f]"),
    estimate,
    low,
    high
  )
}

kdigo_scr_stage <- function(scr, baseline) {

  scr <- safe_numeric(scr)
  baseline <- safe_numeric(baseline)

  out <- rep(NA_integer_, length(scr))
  ok <- is.finite(scr) & is.finite(baseline) & baseline > 0

  if (!any(ok)) return(out)

  ratio <- scr[ok] / baseline[ok]
  delta <- scr[ok] - baseline[ok]

  out[ok] <- ifelse(
    ratio >= 3 | scr[ok] >= 4,
    3L,
    ifelse(
      ratio >= 2,
      2L,
      ifelse(
        ratio >= 1.5 | delta >= 0.3,
        1L,
        0L
      )
    )
  )

  as.integer(out)
}

select_top_n <- function(dt, n_select, score_col = TARGET_PROBABILITY) {

  if (nrow(dt) == 0L || n_select <= 0L) {
    return(dt[0])
  }

  d <- copy(dt)

  setorderv(
    d,
    cols = c(score_col, "stay_id"),
    order = c(-1L, 1L),
    na.last = TRUE
  )

  n_select <- min(as.integer(n_select), nrow(d))
  d[seq_len(n_select)]
}

select_top_fraction <- function(dt, fraction, score_col = TARGET_PROBABILITY) {

  if (nrow(dt) == 0L) return(dt[0])

  n_select <- max(1L, ceiling(nrow(dt) * fraction))

  select_top_n(
    dt = dt,
    n_select = n_select,
    score_col = score_col
  )
}

make_metric_row <- function(
    landmark_name,
    landmark_h,
    strategy,
    reference_label,
    reference_dt,
    selected_dt
) {

  n_reference <- nrow(reference_dt)
  n_selected <- nrow(selected_dt)

  events_reference <- if (n_reference > 0L) {
    sum(reference_dt$outcome_after_landmark == 1L, na.rm = TRUE)
  } else {
    0L
  }

  rate_reference <- if (n_reference > 0L) {
    mean(reference_dt$outcome_after_landmark == 1L, na.rm = TRUE)
  } else {
    NA_real_
  }

  events_selected <- if (n_selected > 0L) {
    sum(selected_dt$outcome_after_landmark == 1L, na.rm = TRUE)
  } else {
    0L
  }

  rate_selected <- if (n_selected > 0L) {
    mean(selected_dt$outcome_after_landmark == 1L, na.rm = TRUE)
  } else {
    NA_real_
  }

  enrichment_factor <- if (
    is.finite(rate_reference) &&
    rate_reference > 0 &&
    is.finite(rate_selected)
  ) {
    rate_selected / rate_reference
  } else {
    NA_real_
  }

  data.table(
    manuscript_day = landmark_name,
    landmark_h = landmark_h,
    strategy = strategy,
    reference_population = reference_label,
    N_reference = n_reference,
    events_reference = events_reference,
    event_rate_reference = rate_reference,
    N_selected = n_selected,
    selected_fraction = if (n_reference > 0L) n_selected / n_reference else NA_real_,
    events_selected = events_selected,
    event_rate_selected = rate_selected,
    enrichment_factor = enrichment_factor
  )
}

compute_refinement_metrics <- function(dt) {

  dt <- as.data.table(copy(dt))

  landmark_name <- unique(dt$manuscript_day)
  landmark_h <- unique(dt$landmark_h)

  if (length(landmark_name) != 1L || length(landmark_h) != 1L) {
    stop("Each input dataset must correspond to exactly one landmark.")
  }

  full_dt <- dt
  rrt_dt <- full_dt[rrt_by_landmark == TRUE]
  kdigo3_max_dt <- full_dt[max_kdigo3_by_landmark == TRUE]

  m2_matched_rrt_dt <- select_top_n(
    full_dt,
    n_select = nrow(rrt_dt)
  )

  non_rrt_dt <- full_dt[rrt_by_landmark == FALSE]

  m2_nonrrt_top10_dt <- select_top_fraction(
    non_rrt_dt,
    fraction = 0.10
  )

  m2_nonrrt_top20_dt <- select_top_fraction(
    non_rrt_dt,
    fraction = 0.20
  )
  
  # ----------------------------------------------------------
  # Conventional renal comparators among NON-RRT patients
  # ----------------------------------------------------------
  
  kdigo3_onset_nonrrt_dt <- non_rrt_dt[
    kdigo3_onset_by_landmark == TRUE
  ]
  
  max_kdigo3_nonrrt_dt <- non_rrt_dt[
    max_kdigo3_by_landmark == TRUE
  ]

  rbindlist(
    list(
      make_metric_row(
        landmark_name,
        landmark_h,
        "RRT initiated by landmark",
        "All landmark-eligible patients",
        full_dt,
        rrt_dt
      ),
      make_metric_row(
        landmark_name,
        landmark_h,
        "Maximum KDIGO stage 3 reached by landmark",
        "All landmark-eligible patients",
        full_dt,
        kdigo3_max_dt
      ),
      make_metric_row(
        landmark_name,
        landmark_h,
        "M2 matched to RRT-selected fraction",
        "All landmark-eligible patients",
        full_dt,
        m2_matched_rrt_dt
      ),
      make_metric_row(
        landmark_name = landmark_name,
        landmark_h = landmark_h,
        strategy = "KDIGO stage 3 at AKI onset among non-RRT patients",
        reference_label = "Landmark-eligible patients not on RRT",
        reference_dt = non_rrt_dt,
        selected_dt = kdigo3_onset_nonrrt_dt
      ),
      
      make_metric_row(
        landmark_name = landmark_name,
        landmark_h = landmark_h,
        strategy = "Maximum SCr-based KDIGO stage 3 reached by landmark among non-RRT patients",
        reference_label = "Landmark-eligible patients not on RRT",
        reference_dt = non_rrt_dt,
        selected_dt = max_kdigo3_nonrrt_dt
      ),
      make_metric_row(
        landmark_name,
        landmark_h,
        "M2 top 10% among non-RRT patients",
        "Landmark-eligible patients not on RRT",
        non_rrt_dt,
        m2_nonrrt_top10_dt
      ),
      make_metric_row(
        landmark_name,
        landmark_h,
        "M2 top 20% among non-RRT patients",
        "Landmark-eligible patients not on RRT",
        non_rrt_dt,
        m2_nonrrt_top20_dt
      )
    ),
    use.names = TRUE,
    fill = TRUE
  )
}

# ------------------------------------------------------------
# 2. LOAD BASELINE + PREDICTIONS
# ------------------------------------------------------------

stopifnot(file.exists(MIMIC_BF_RDS))
stopifnot(file.exists(PRED_RDS) || file.exists(PRED_CSV))

cat("=== Loading MIMIC-IV baseline data ===\n")

bf_M <- as.data.table(readRDS(MIMIC_BF_RDS))

required_baseline_columns <- c(
  "stay_id",
  "intime",
  "outtime",
  "deathtime",
  "scr_baseline",
  "aki_onset_time",
  "aki_stage_at_onset",
  "rrt_start_time"
)

missing_baseline_columns <- setdiff(
  required_baseline_columns,
  names(bf_M)
)

if (length(missing_baseline_columns) > 0L) {
  stop(
    "Missing baseline columns: ",
    paste(missing_baseline_columns, collapse = ", ")
  )
}

bf_M[, stay_id := as.integer(stay_id)]
bf_M <- unique(bf_M, by = "stay_id")

bf_M[, intime := to_posix(intime)]
bf_M[, outtime := to_posix(outtime)]
bf_M[, deathtime := to_posix(deathtime)]
bf_M[, aki_onset_time := to_posix(aki_onset_time)]
bf_M[, rrt_start_time := to_posix(rrt_start_time)]
bf_M[, scr_baseline := safe_numeric(scr_baseline)]

bf_M[, death_icu := as.integer(
  !is.na(deathtime) &
    !is.na(intime) &
    !is.na(outtime) &
    deathtime >= intime &
    deathtime <= outtime
)]

cat("=== Loading saved M2 predictions ===\n")

if (file.exists(PRED_RDS)) {
  pred_dt <- as.data.table(readRDS(PRED_RDS))
} else {
  pred_dt <- fread(PRED_CSV)
}

required_prediction_columns <- c(
  "id",
  "model",
  "day",
  "true_class",
  TARGET_PROBABILITY
)

missing_prediction_columns <- setdiff(
  required_prediction_columns,
  names(pred_dt)
)

if (length(missing_prediction_columns) > 0L) {
  stop(
    "Missing prediction columns: ",
    paste(missing_prediction_columns, collapse = ", ")
  )
}

pred_dt[, id := as.integer(id)]
pred_dt[, day := as.integer(day)]

# ------------------------------------------------------------
# 3. LOAD LONGITUDINAL SCr
# ------------------------------------------------------------

existing_ts_paths <- TS_CANDIDATES[file.exists(TS_CANDIDATES)]

if (length(existing_ts_paths) == 0L) {
  stop(
    paste0(
      "No longitudinal MIMIC time-series file found.\n",
      "Checked:\n",
      paste(TS_CANDIDATES, collapse = "\n")
    )
  )
}

MIMIC_TS_RDS <- existing_ts_paths[1]

cat("=== Loading longitudinal MIMIC time series ===\n")
cat("Using: ", MIMIC_TS_RDS, "\n", sep = "")

ts_M <- as.data.table(readRDS(MIMIC_TS_RDS))

required_ts_columns <- c("stay_id", "t_hours", "var")

missing_ts_columns <- setdiff(required_ts_columns, names(ts_M))

if (length(missing_ts_columns) > 0L) {
  stop(
    "Time-series file is missing: ",
    paste(missing_ts_columns, collapse = ", ")
  )
}

value_col <- if ("value" %in% names(ts_M)) {
  "value"
} else if ("value_num" %in% names(ts_M)) {
  "value_num"
} else if ("med" %in% names(ts_M)) {
  "med"
} else {
  NA_character_
}

if (is.na(value_col)) {
  stop(
    "Could not identify time-series value column. ",
    "Expected one of: value, value_num, med."
  )
}

ts_M[, stay_id := as.integer(stay_id)]
ts_M[, t_hours := safe_numeric(t_hours)]
ts_M[, value_num_tmp := safe_numeric(get(value_col))]

scr_ts <- ts_M[
  var == "scr" &
    is.finite(t_hours) &
    t_hours >= 0 &
    t_hours <= 48 &
    is.finite(value_num_tmp),
  .(
    stay_id,
    t_hours,
    scr = value_num_tmp
  )
]

rm(ts_M)
gc()

cat(
  "SCr rows available from 0 to 48 h: ",
  format(nrow(scr_ts), big.mark = ","),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 4. MAXIMUM KDIGO STAGE 3 REACHED BY LANDMARK
# ------------------------------------------------------------

scr_ts <- merge(
  scr_ts,
  bf_M[, .(stay_id, scr_baseline)],
  by = "stay_id",
  all.x = TRUE
)

scr_ts[, kdigo_scr_stage := kdigo_scr_stage(
  scr = scr,
  baseline = scr_baseline
)]

kdigo_stage3_by_landmark <- rbindlist(
  lapply(
    seq_len(nrow(LANDMARK_MAP)),
    function(i) {

      map_i <- LANDMARK_MAP[i]

      tmp <- scr_ts[
        t_hours <= map_i$landmark_h,
        .(
          scr_stage3_by_landmark = any(
            kdigo_scr_stage == 3L,
            na.rm = TRUE
          )
        ),
        by = stay_id
      ]

      tmp[, `:=`(
        manuscript_day = map_i$manuscript_day,
        landmark_h = map_i$landmark_h
      )]

      tmp
    }
  ),
  use.names = TRUE,
  fill = TRUE
)

# ------------------------------------------------------------
# 5. BUILD LANDMARK RISK SETS
# ------------------------------------------------------------

# Primary landmark datasets:
# patients alive and still in ICU at the landmark, regardless of
# whether SCr-defined AKI had already been documented.
landmark_datasets <- list()

# Sensitivity landmark datasets requested by Reviewer 2:
# same landmark eligibility criteria, additionally restricted to
# patients who had already fulfilled SCr-based KDIGO AKI criteria
# by the landmark.
landmark_datasets_aki_by_landmark <- list()

audit_list <- list()

for (i in seq_len(nrow(LANDMARK_MAP))) {
  
  map_i <- LANDMARK_MAP[i]
  
  # ----------------------------------------------------------
  # Saved M2 predictions for the corresponding landmark
  # ----------------------------------------------------------
  
  pred_slice <- pred_dt[
    model == MODEL_USE &
      day == map_i$pred_day_code,
    .(
      stay_id = as.integer(id[1]),
      true_class = as.character(true_class[1]),
      pC2 = mean(as.numeric(pC2), na.rm = TRUE)
    ),
    by = id
  ][, id := NULL]
  
  pred_slice <- unique(pred_slice, by = "stay_id")
  
  # ----------------------------------------------------------
  # Merge predictions with baseline / outcome information
  # ----------------------------------------------------------
  
  dat_i <- merge(
    pred_slice,
    bf_M[
      ,
      .(
        stay_id,
        intime,
        outtime,
        deathtime,
        death_icu,
        scr_baseline,
        aki_onset_time,
        aki_stage_at_onset,
        rrt_start_time
      )
    ],
    by = "stay_id",
    all.x = TRUE
  )
  
  # ----------------------------------------------------------
  # Landmark time
  # ----------------------------------------------------------
  
  dat_i[, landmark_time :=
          intime + map_i$landmark_h * 3600
  ]
  
  # ----------------------------------------------------------
  # Primary landmark eligibility:
  # alive and still in ICU at the landmark
  # ----------------------------------------------------------
  
  dat_i[, eligible_at_landmark :=
          !is.na(intime) &
          !is.na(outtime) &
          outtime > landmark_time &
          (is.na(deathtime) | deathtime > landmark_time)
  ]
  
  
  # Audit des exclusions H24 : MIMIC-IV
  if (map_i$landmark_h == 24L) {
    
    audit_h24 <- copy(dat_i)
    
    
    # Vérifier les patients absents avant la construction du risk set
    cohort_ids <- unique(bf_M$stay_id)
    prediction_ids <- unique(pred_slice$stay_id)
    
    missing_from_predictions <- setdiff(
      cohort_ids,
      prediction_ids
    )
    
    cat("\n=== COHORT VS M2 PREDICTIONS ===\n")
    cat("Patients in baseline dataset:", length(cohort_ids), "\n")
    cat("Patients with M2 predictions:", length(prediction_ids), "\n")
    cat("Patients without M2 predictions:",
        length(missing_from_predictions), "\n")
    
    # Effectifs exclus après construction de dat_i
    cat("\n=== LANDMARK ELIGIBILITY ===\n")
    cat("N before landmark filtering:", nrow(dat_i), "\n")
    cat("N eligible by timestamps:",
        sum(dat_i$eligible_at_landmark, na.rm = TRUE), "\n")
    
    
    audit_h24[, `:=`(
      hours_to_death = as.numeric(
        difftime(deathtime, intime, units = "hours")
      ),
      hours_to_outtime = as.numeric(
        difftime(outtime, intime, units = "hours")
      ),
      
      death_by_h24 = !is.na(deathtime) &
        deathtime <= landmark_time,
      
      discharge_by_h24 = !is.na(outtime) &
        outtime <= landmark_time,
      
      missing_prediction = !is.finite(pC2),
      
      missing_class = is.na(true_class)
    )]
    
    audit_h24[, exclusion_reason := fifelse(
      death_by_h24,
      "Death at or before H24",
      fifelse(
        discharge_by_h24,
        "ICU discharge at or before H24",
        fifelse(
          missing_prediction | missing_class,
          "Missing prediction or class",
          fifelse(
            is.na(intime) | is.na(outtime),
            "Missing timestamps",
            "Eligible"
          )
        )
      )
    )]
    
    cat("\n=== H24 EXCLUSION REASONS ===\n")
    
    print(
      audit_h24[, .N, by = exclusion_reason][
        order(exclusion_reason)
      ]
    )
    
    cat("\n=== DEATHS AT OR BEFORE H24 ===\n")
    
    print(
      audit_h24[
        death_by_h24 == TRUE,
        .(
          stay_id,
          intime,
          deathtime,
          outtime,
          hours_to_death,
          hours_to_outtime,
          eligible_at_landmark
        )
      ][order(hours_to_death)]
    )
    
    cat("\n=== DEATHS H24 WITH OUTTIME > H24 ===\n")
    
    print(
      audit_h24[
        death_by_h24 == TRUE &
          hours_to_outtime > 24,
        .(
          N = .N,
          min_death_h = min(hours_to_death, na.rm = TRUE),
          median_death_h = median(hours_to_death, na.rm = TRUE),
          max_death_h = max(hours_to_death, na.rm = TRUE)
        )
      ]
    )
    cat("\n=== RECONCILIATION COHORT -> H24 ===\n")
    
    cat("Initial cohort:", length(cohort_ids), "\n")
    
    cat("Absent from M2 predictions:",
        length(missing_from_predictions), "\n")
    
    cat("Excluded by timestamps:",
        sum(!dat_i$eligible_at_landmark, na.rm = TRUE), "\n")
    
    cat("Excluded by prediction or class:",
        sum(
          dat_i$eligible_at_landmark &
            (!is.finite(dat_i$pC2) | is.na(dat_i$true_class)),
          na.rm = TRUE
        ), "\n")
    
    cat("Final H24 population:",
        sum(
          dat_i$eligible_at_landmark &
            is.finite(dat_i$pC2) &
            !is.na(dat_i$true_class),
          na.rm = TRUE
        ), "\n")
    
    fwrite(
      audit_h24,
      file.path(OUT_DIR, "audit_exclusions_H24_MIMIC.csv")
    )
  }
  
  # ----------------------------------------------------------
  # NEW:
  # SCr-defined AKI already documented by the landmark
  #
  # This defines the prospectively screenable AKI population
  # requested by the Reviewer.
  # ----------------------------------------------------------
  
  dat_i[, aki_by_landmark :=
          !is.na(aki_onset_time) &
          !is.na(intime) &
          aki_onset_time >= intime &
          aki_onset_time <= landmark_time
  ]
  
  # ----------------------------------------------------------
  # Subsequent ICU mortality
  # ----------------------------------------------------------
  
  dat_i[, outcome_after_landmark := as.integer(
    death_icu == 1L &
      !is.na(deathtime) &
      deathtime > landmark_time
  )]
  
  # ----------------------------------------------------------
  # RRT initiated by landmark
  # ----------------------------------------------------------
  
  dat_i[, rrt_by_landmark :=
          !is.na(rrt_start_time) &
          rrt_start_time >= intime &
          rrt_start_time <= landmark_time
  ]
  
  # ----------------------------------------------------------
  # KDIGO stage 3 at first AKI onset
  # ----------------------------------------------------------
  
  dat_i[, kdigo3_onset_by_landmark :=
          !is.na(aki_onset_time) &
          aki_onset_time >= intime &
          aki_onset_time <= landmark_time &
          aki_stage_at_onset == 3
  ]
  
  # ----------------------------------------------------------
  # Maximum SCr-based KDIGO stage reached by landmark
  # ----------------------------------------------------------
  
  max_kdigo_i <- kdigo_stage3_by_landmark[
    manuscript_day == map_i$manuscript_day,
    .(
      stay_id,
      scr_stage3_by_landmark
    )
  ]
  
  dat_i <- merge(
    dat_i,
    max_kdigo_i,
    by = "stay_id",
    all.x = TRUE
  )
  
  dat_i[
    is.na(scr_stage3_by_landmark),
    scr_stage3_by_landmark := FALSE
  ]
  
  dat_i[, max_kdigo3_by_landmark :=
          scr_stage3_by_landmark |
          rrt_by_landmark
  ]
  
  # ----------------------------------------------------------
  # PRIMARY LANDMARK RISK SET
  #
  # This reproduces the original analyses:
  # patients alive and in ICU at landmark, with available
  # trajectory prediction, irrespective of whether AKI had
  # already been documented by that landmark.
  # ----------------------------------------------------------
  
  eligible_i <- dat_i[
    eligible_at_landmark == TRUE &
      !is.na(true_class) &
      is.finite(pC2) &
      !is.na(outcome_after_landmark)
  ]
  
  setorder(eligible_i, stay_id)
  
  if (anyDuplicated(eligible_i$stay_id) > 0L) {
    stop(
      "Duplicate stay_id values remain in primary risk set at ",
      map_i$manuscript_day
    )
  }
  
  eligible_i[, `:=`(
    manuscript_day = map_i$manuscript_day,
    landmark_h = map_i$landmark_h,
    pred_day_code = map_i$pred_day_code,
    observation_window = map_i$observation_window
  )]
  
  landmark_datasets[[map_i$manuscript_day]] <- eligible_i
  
  # ----------------------------------------------------------
  # NEW SENSITIVITY RISK SET:
  # patients with SCr-defined AKI already established
  # by the landmark
  # ----------------------------------------------------------
  
  eligible_aki_i <- eligible_i[
    aki_by_landmark == TRUE
  ]
  
  setorder(eligible_aki_i, stay_id)
  
  if (anyDuplicated(eligible_aki_i$stay_id) > 0L) {
    stop(
      "Duplicate stay_id values remain in AKI-by-landmark risk set at ",
      map_i$manuscript_day
    )
  }
  
  landmark_datasets_aki_by_landmark[[map_i$manuscript_day]] <- eligible_aki_i
  
  # ----------------------------------------------------------
  # AUDIT
  # ----------------------------------------------------------
  
  audit_i <- data.table(
    manuscript_day = map_i$manuscript_day,
    landmark_h = map_i$landmark_h,
    
    # Original landmark population
    N_eligible = nrow(eligible_i),
    
    subsequent_ICU_deaths =
      sum(
        eligible_i$outcome_after_landmark == 1L,
        na.rm = TRUE
      ),
    
    subsequent_ICU_mortality =
      mean(
        eligible_i$outcome_after_landmark == 1L,
        na.rm = TRUE
      ),
    
    # NEW: SCr-AKI established by landmark
    N_AKI_by_landmark = nrow(eligible_aki_i),
    
    proportion_AKI_by_landmark =
      if (nrow(eligible_i) > 0L) {
        nrow(eligible_aki_i) / nrow(eligible_i)
      } else {
        NA_real_
      },
    
    N_not_yet_AKI_at_landmark =
      nrow(eligible_i) - nrow(eligible_aki_i),
    
    subsequent_ICU_deaths_AKI_by_landmark =
      sum(
        eligible_aki_i$outcome_after_landmark == 1L,
        na.rm = TRUE
      ),
    
    subsequent_ICU_mortality_AKI_by_landmark =
      if (nrow(eligible_aki_i) > 0L) {
        mean(
          eligible_aki_i$outcome_after_landmark == 1L,
          na.rm = TRUE
        )
      } else {
        NA_real_
      },
    
    # Existing renal comparators in the primary risk set
    old_KDIGO3_at_onset_N =
      sum(
        eligible_i$kdigo3_onset_by_landmark == TRUE,
        na.rm = TRUE
      ),
    
    RRT_by_landmark_N =
      sum(
        eligible_i$rrt_by_landmark == TRUE,
        na.rm = TRUE
      ),
    
    new_max_KDIGO3_by_landmark_N =
      sum(
        eligible_i$max_kdigo3_by_landmark == TRUE,
        na.rm = TRUE
      ),
    
    SCr_based_stage3_by_landmark_N =
      sum(
        eligible_i$scr_stage3_by_landmark == TRUE,
        na.rm = TRUE
      ),
    
    # Same renal comparators within the restricted
    # AKI-by-landmark population
    RRT_by_landmark_AKI_population_N =
      sum(
        eligible_aki_i$rrt_by_landmark == TRUE,
        na.rm = TRUE
      ),
    
    max_KDIGO3_by_landmark_AKI_population_N =
      sum(
        eligible_aki_i$max_kdigo3_by_landmark == TRUE,
        na.rm = TRUE
      )
  )
  
  audit_list[[map_i$manuscript_day]] <- audit_i
  
  # ----------------------------------------------------------
  # DISPLAY
  # ----------------------------------------------------------
  
  cat("\n============================================\n")
  cat(
    map_i$manuscript_day,
    " / H",
    map_i$landmark_h,
    "\n",
    sep = ""
  )
  cat("============================================\n")
  
  cat("\nPRIMARY LANDMARK RISK SET\n")
  cat(
    "N = ",
    nrow(eligible_i),
    "\n",
    sep = ""
  )
  cat(
    "Subsequent ICU deaths = ",
    sum(
      eligible_i$outcome_after_landmark == 1L,
      na.rm = TRUE
    ),
    "\n",
    sep = ""
  )
  cat(
    "Subsequent ICU mortality = ",
    sprintf(
      "%.1f%%",
      100 *
        mean(
          eligible_i$outcome_after_landmark == 1L,
          na.rm = TRUE
        )
    ),
    "\n",
    sep = ""
  )
  
  cat("\nSCr-AKI ESTABLISHED BY LANDMARK\n")
  cat(
    "N = ",
    nrow(eligible_aki_i),
    " / ",
    nrow(eligible_i),
    " (",
    sprintf(
      "%.1f%%",
      100 * nrow(eligible_aki_i) / nrow(eligible_i)
    ),
    ")\n",
    sep = ""
  )
  
  cat(
    "Not yet SCr-AKI at landmark = ",
    nrow(eligible_i) - nrow(eligible_aki_i),
    "\n",
    sep = ""
  )
  
  cat(
    "Subsequent ICU deaths = ",
    sum(
      eligible_aki_i$outcome_after_landmark == 1L,
      na.rm = TRUE
    ),
    "\n",
    sep = ""
  )
  
  cat(
    "Subsequent ICU mortality = ",
    sprintf(
      "%.1f%%",
      100 *
        mean(
          eligible_aki_i$outcome_after_landmark == 1L,
          na.rm = TRUE
        )
    ),
    "\n",
    sep = ""
  )
  
  cat("\nFULL AUDIT ROW\n")
  print(audit_i)
  
}

# ------------------------------------------------------------
# Combined audit table
# ------------------------------------------------------------

audit_table <- rbindlist(
  audit_list,
  use.names = TRUE,
  fill = TRUE
)

# ------------------------------------------------------------
# 6. POINT ESTIMATES
# ------------------------------------------------------------

# ============================================================
# 6A. Helper for AKI-by-landmark sensitivity analysis
# ============================================================

compute_aki_landmark_sensitivity_metrics <- function(dt) {
  
  dt <- as.data.table(copy(dt))
  
  landmark_name <- unique(dt$manuscript_day)
  landmark_h <- unique(dt$landmark_h)
  
  if (
    length(landmark_name) != 1L ||
    length(landmark_h) != 1L
  ) {
    stop(
      "Each input dataset must correspond to exactly one landmark."
    )
  }
  
  if (nrow(dt) == 0L) {
    stop(
      "Empty AKI-by-landmark risk set at ",
      landmark_name
    )
  }
  
  # ----------------------------------------------------------
  # Reference population:
  # SCr-defined AKI already established by landmark
  # ----------------------------------------------------------
  
  full_dt <- dt
  
  # ----------------------------------------------------------
  # Retrospective JLCM class benchmark
  # true_class is expected to be C1 / C2 / C3.
  # The additional "2" option makes this robust if labels
  # have been converted to numeric strings.
  # ----------------------------------------------------------
  
  progressive_dt <- full_dt[
    as.character(true_class) %in% c("C2", "2")
  ]
  
  # ----------------------------------------------------------
  # M2 model-based selections in the FULL restricted population
  # ----------------------------------------------------------
  
  m2_top10_dt <- select_top_fraction(
    full_dt,
    fraction = 0.10
  )
  
  m2_top20_dt <- select_top_fraction(
    full_dt,
    fraction = 0.20
  )
  
  # ----------------------------------------------------------
  # Conventional renal comparators
  # ----------------------------------------------------------
  
  rrt_dt <- full_dt[
    rrt_by_landmark == TRUE
  ]
  
  kdigo3_onset_dt <- full_dt[
    kdigo3_onset_by_landmark == TRUE
  ]
  
  kdigo3_max_dt <- full_dt[
    max_kdigo3_by_landmark == TRUE
  ]
  
  # ----------------------------------------------------------
  # M2 selection matched to exact RRT-selected N
  # ----------------------------------------------------------
  
  m2_matched_rrt_dt <- select_top_n(
    full_dt,
    n_select = nrow(rrt_dt)
  )
  
  # ----------------------------------------------------------
  # Non-RRT population
  # ----------------------------------------------------------
  
  non_rrt_dt <- full_dt[
    rrt_by_landmark == FALSE
  ]
  
  kdigo3_onset_nonrrt_dt <- non_rrt_dt[
    kdigo3_onset_by_landmark == TRUE
  ]
  
  kdigo3_max_nonrrt_dt <- non_rrt_dt[
    max_kdigo3_by_landmark == TRUE
  ]
  
  m2_nonrrt_top10_dt <- select_top_fraction(
    non_rrt_dt,
    fraction = 0.10
  )
  
  m2_nonrrt_top20_dt <- select_top_fraction(
    non_rrt_dt,
    fraction = 0.20
  )
  
  # ----------------------------------------------------------
  # Assemble results
  # ----------------------------------------------------------
  
  rbindlist(
    list(
      
      # ------------------------------------------------------
      # Retrospective benchmark
      # ------------------------------------------------------
      
      make_metric_row(
        landmark_name = landmark_name,
        landmark_h = landmark_h,
        strategy = "JLCM-assigned Progressive AKI archetype",
        reference_label =
          "Patients with SCr-defined AKI established by landmark",
        reference_dt = full_dt,
        selected_dt = progressive_dt
      ),
      
      # ------------------------------------------------------
      # M2 operational selections
      # ------------------------------------------------------
      
      make_metric_row(
        landmark_name = landmark_name,
        landmark_h = landmark_h,
        strategy = "M2 top 10%",
        reference_label =
          "Patients with SCr-defined AKI established by landmark",
        reference_dt = full_dt,
        selected_dt = m2_top10_dt
      ),
      
      make_metric_row(
        landmark_name = landmark_name,
        landmark_h = landmark_h,
        strategy = "M2 top 20%",
        reference_label =
          "Patients with SCr-defined AKI established by landmark",
        reference_dt = full_dt,
        selected_dt = m2_top20_dt
      ),
      
      # ------------------------------------------------------
      # Conventional renal criteria
      # ------------------------------------------------------
      
      make_metric_row(
        landmark_name = landmark_name,
        landmark_h = landmark_h,
        strategy = "KDIGO stage 3 at AKI onset",
        reference_label =
          "Patients with SCr-defined AKI established by landmark",
        reference_dt = full_dt,
        selected_dt = kdigo3_onset_dt
      ),
      
      make_metric_row(
        landmark_name = landmark_name,
        landmark_h = landmark_h,
        strategy =
          "Maximum SCr-based KDIGO stage 3 reached by landmark",
        reference_label =
          "Patients with SCr-defined AKI established by landmark",
        reference_dt = full_dt,
        selected_dt = kdigo3_max_dt
      ),
      
      make_metric_row(
        landmark_name = landmark_name,
        landmark_h = landmark_h,
        strategy = "RRT initiated by landmark",
        reference_label =
          "Patients with SCr-defined AKI established by landmark",
        reference_dt = full_dt,
        selected_dt = rrt_dt
      ),
      
      make_metric_row(
        landmark_name = landmark_name,
        landmark_h = landmark_h,
        strategy = "M2 matched to RRT-selected fraction",
        reference_label =
          "Patients with SCr-defined AKI established by landmark",
        reference_dt = full_dt,
        selected_dt = m2_matched_rrt_dt
      ),
      
      # ------------------------------------------------------
      # Non-RRT sensitivity analyses
      # ------------------------------------------------------
      
      make_metric_row(
        landmark_name = landmark_name,
        landmark_h = landmark_h,
        strategy =
          "KDIGO stage 3 at AKI onset among non-RRT patients",
        reference_label =
          "Patients with SCr-defined AKI established by landmark and not on RRT",
        reference_dt = non_rrt_dt,
        selected_dt = kdigo3_onset_nonrrt_dt
      ),
      
      make_metric_row(
        landmark_name = landmark_name,
        landmark_h = landmark_h,
        strategy =
          "Maximum SCr-based KDIGO stage 3 reached by landmark among non-RRT patients",
        reference_label =
          "Patients with SCr-defined AKI established by landmark and not on RRT",
        reference_dt = non_rrt_dt,
        selected_dt = kdigo3_max_nonrrt_dt
      ),
      
      make_metric_row(
        landmark_name = landmark_name,
        landmark_h = landmark_h,
        strategy = "M2 top 10% among non-RRT patients",
        reference_label =
          "Patients with SCr-defined AKI established by landmark and not on RRT",
        reference_dt = non_rrt_dt,
        selected_dt = m2_nonrrt_top10_dt
      ),
      
      make_metric_row(
        landmark_name = landmark_name,
        landmark_h = landmark_h,
        strategy = "M2 top 20% among non-RRT patients",
        reference_label =
          "Patients with SCr-defined AKI established by landmark and not on RRT",
        reference_dt = non_rrt_dt,
        selected_dt = m2_nonrrt_top20_dt
      )
    ),
    use.names = TRUE,
    fill = TRUE
  )
}


# ============================================================
# 6B. ORIGINAL analysis — unchanged
# ============================================================

point_results <- rbindlist(
  lapply(
    landmark_datasets,
    compute_refinement_metrics
  ),
  use.names = TRUE,
  fill = TRUE
)

cat("\n============================================\n")
cat("POINT ESTIMATES — ORIGINAL LANDMARK POPULATION\n")
cat("============================================\n")
print(point_results)


# ============================================================
# 6C. NEW sensitivity analysis:
# AKI already established by landmark
# ============================================================

point_results_aki_by_landmark <- rbindlist(
  lapply(
    landmark_datasets_aki_by_landmark,
    compute_aki_landmark_sensitivity_metrics
  ),
  use.names = TRUE,
  fill = TRUE
)

cat("\n============================================\n")
cat("POINT ESTIMATES — AKI ESTABLISHED BY LANDMARK\n")
cat("============================================\n")
print(point_results_aki_by_landmark)


# ------------------------------------------------------------
# 7. PATIENT-LEVEL BOOTSTRAP
# ------------------------------------------------------------

# ============================================================
# 7A. Generic bootstrap helper
# ============================================================

run_landmark_bootstrap <- function(
    dataset_list,
    metric_function,
    analysis_label,
    seed_offset = 0L
) {
  
  output_list <- list()
  
  for (landmark_name in names(dataset_list)) {
    
    dat_i <- dataset_list[[landmark_name]]
    
    landmark_h_i <- unique(dat_i$landmark_h)
    
    if (length(landmark_h_i) != 1L) {
      stop(
        "Invalid landmark_h in ",
        analysis_label,
        " / ",
        landmark_name
      )
    }
    
    set.seed(
      BOOT_SEED +
        landmark_h_i +
        seed_offset
    )
    
    cat(
      "\n============================================\n"
    )
    cat(
      "Bootstrap — ",
      analysis_label,
      " — ",
      landmark_name,
      " (",
      N_BOOT,
      " replicates)\n",
      sep = ""
    )
    cat(
      "============================================\n"
    )
    
    n_patients <- nrow(dat_i)
    
    if (n_patients == 0L) {
      stop(
        "Empty dataset supplied to bootstrap: ",
        analysis_label,
        " / ",
        landmark_name
      )
    }
    
    boot_i_list <- vector(
      "list",
      N_BOOT
    )
    
    for (b in seq_len(N_BOOT)) {
      
      idx <- sample.int(
        n = n_patients,
        size = n_patients,
        replace = TRUE
      )
      
      boot_sample <- dat_i[idx]
      
      res_b <- metric_function(
        boot_sample
      )
      
      res_b[, bootstrap_iteration := b]
      
      boot_i_list[[b]] <- res_b
      
      if (b %% 100L == 0L) {
        
        cat(
          analysis_label,
          " — ",
          landmark_name,
          ": completed ",
          b,
          " / ",
          N_BOOT,
          "\n",
          sep = ""
        )
      }
    }
    
    output_list[[landmark_name]] <- rbindlist(
      boot_i_list,
      use.names = TRUE,
      fill = TRUE
    )
  }
  
  rbindlist(
    output_list,
    use.names = TRUE,
    fill = TRUE
  )
}


# ============================================================
# 7B. ORIGINAL bootstrap
# ============================================================

bootstrap_results <- run_landmark_bootstrap(
  dataset_list = landmark_datasets,
  metric_function = compute_refinement_metrics,
  analysis_label = "Original landmark population",
  seed_offset = 0L
)


# ============================================================
# 7C. NEW AKI-by-landmark bootstrap
# ============================================================

bootstrap_results_aki_by_landmark <-
  run_landmark_bootstrap(
    dataset_list =
      landmark_datasets_aki_by_landmark,
    metric_function =
      compute_aki_landmark_sensitivity_metrics,
    analysis_label =
      "SCr-AKI established by landmark",
    seed_offset = 10000L
  )


# ------------------------------------------------------------
# 8. BOOTSTRAP 95% CIs
# ------------------------------------------------------------

# ============================================================
# 8A. Generic CI helper
# ============================================================

compute_bootstrap_ci <- function(bootstrap_dt) {
  
  bootstrap_dt[
    ,
    .(
      selected_fraction_low =
        safe_quantile(
          selected_fraction,
          0.025
        ),
      
      selected_fraction_high =
        safe_quantile(
          selected_fraction,
          0.975
        ),
      
      event_rate_reference_low =
        safe_quantile(
          event_rate_reference,
          0.025
        ),
      
      event_rate_reference_high =
        safe_quantile(
          event_rate_reference,
          0.975
        ),
      
      event_rate_selected_low =
        safe_quantile(
          event_rate_selected,
          0.025
        ),
      
      event_rate_selected_high =
        safe_quantile(
          event_rate_selected,
          0.975
        ),
      
      enrichment_factor_low =
        safe_quantile(
          enrichment_factor,
          0.025
        ),
      
      enrichment_factor_high =
        safe_quantile(
          enrichment_factor,
          0.975
        )
    ),
    by = .(
      manuscript_day,
      landmark_h,
      strategy,
      reference_population
    )
  ]
}


# ============================================================
# 8B. ORIGINAL bootstrap CIs
# ============================================================

bootstrap_ci <- compute_bootstrap_ci(
  bootstrap_results
)

final_results <- merge(
  point_results,
  bootstrap_ci,
  by = c(
    "manuscript_day",
    "landmark_h",
    "strategy",
    "reference_population"
  ),
  all.x = TRUE
)


# ============================================================
# 8C. NEW AKI-by-landmark bootstrap CIs
# ============================================================

bootstrap_ci_aki_by_landmark <-
  compute_bootstrap_ci(
    bootstrap_results_aki_by_landmark
  )

final_results_aki_by_landmark <- merge(
  point_results_aki_by_landmark,
  bootstrap_ci_aki_by_landmark,
  by = c(
    "manuscript_day",
    "landmark_h",
    "strategy",
    "reference_population"
  ),
  all.x = TRUE
)


# ------------------------------------------------------------
# 9. MANUSCRIPT-READY TABLES
# ------------------------------------------------------------

# ============================================================
# 9A. Generic table formatter
# ============================================================

build_manuscript_table <- function(
    results_dt,
    strategy_order
) {
  
  x <- copy(results_dt)
  
  x[
    ,
    landmark_order := match(
      manuscript_day,
      c("Day 1", "Day 2")
    )
  ]
  
  x[
    ,
    strategy_order_num := match(
      strategy,
      strategy_order
    )
  ]
  
  setorder(
    x,
    landmark_order,
    strategy_order_num
  )
  
  x[
    ,
    c(
      "landmark_order",
      "strategy_order_num"
    ) := NULL
  ]
  
  x[
    ,
    .(
      Landmark =
        paste0("H", landmark_h),
      
      Strategy =
        strategy,
      
      `Reference population` =
        reference_population,
      
      `N reference` =
        N_reference,
      
      `Subsequent ICU deaths, reference` =
        events_reference,
      
      `ICU mortality, reference (95% CI)` =
        fmt_pct_ci(
          event_rate_reference,
          event_rate_reference_low,
          event_rate_reference_high
        ),
      
      `N selected` =
        N_selected,
      
      `Selected fraction (95% CI)` =
        fmt_pct_ci(
          selected_fraction,
          selected_fraction_low,
          selected_fraction_high
        ),
      
      `Subsequent ICU deaths, selected` =
        events_selected,
      
      `ICU mortality, selected (95% CI)` =
        fmt_pct_ci(
          event_rate_selected,
          event_rate_selected_low,
          event_rate_selected_high
        ),
      
      `Enrichment factor (95% CI)` =
        fmt_num_ci(
          enrichment_factor,
          enrichment_factor_low,
          enrichment_factor_high,
          digits = 2
        )
    )
  ]
}


# ============================================================
# 9B. ORIGINAL table
# ============================================================

strategy_order_original <- c(
  "RRT initiated by landmark",
  "Maximum KDIGO stage 3 reached by landmark",
  "M2 matched to RRT-selected fraction",
  "KDIGO stage 3 at AKI onset among non-RRT patients",
  "Maximum SCr-based KDIGO stage 3 reached by landmark among non-RRT patients",
  "M2 top 10% among non-RRT patients",
  "M2 top 20% among non-RRT patients"
)

manuscript_table <- build_manuscript_table(
  final_results,
  strategy_order_original
)


# ============================================================
# 9C. NEW AKI-by-landmark table
# ============================================================

strategy_order_aki_by_landmark <- c(
  "JLCM-assigned Progressive AKI archetype",
  "M2 top 10%",
  "M2 top 20%",
  "KDIGO stage 3 at AKI onset",
  "Maximum SCr-based KDIGO stage 3 reached by landmark",
  "RRT initiated by landmark",
  "M2 matched to RRT-selected fraction",
  "KDIGO stage 3 at AKI onset among non-RRT patients",
  "Maximum SCr-based KDIGO stage 3 reached by landmark among non-RRT patients",
  "M2 top 10% among non-RRT patients",
  "M2 top 20% among non-RRT patients"
)

manuscript_table_aki_by_landmark <-
  build_manuscript_table(
    final_results_aki_by_landmark,
    strategy_order_aki_by_landmark
  )


# ------------------------------------------------------------
# 10. EXPORTS
# ------------------------------------------------------------

# ============================================================
# 10A. Audit
# ============================================================

fwrite(
  audit_table,
  file.path(
    OUT_DIR,
    "audit_landmarks_and_AKI_by_landmark.csv"
  )
)


# ============================================================
# 10B. ORIGINAL analysis
# ============================================================

fwrite(
  final_results,
  file.path(
    OUT_DIR,
    "reviewer2_enrichment_refinements_full.csv"
  )
)

fwrite(
  manuscript_table,
  file.path(
    OUT_DIR,
    "reviewer2_enrichment_refinements_manuscript_table.csv"
  )
)

saveRDS(
  final_results,
  file.path(
    OUT_DIR,
    "reviewer2_enrichment_refinements_full.rds"
  )
)


# ============================================================
# 10C. NEW AKI-by-landmark sensitivity analysis
# ============================================================

fwrite(
  final_results_aki_by_landmark,
  file.path(
    OUT_DIR,
    "reviewer2_AKI_by_landmark_sensitivity_full.csv"
  )
)

fwrite(
  manuscript_table_aki_by_landmark,
  file.path(
    OUT_DIR,
    "reviewer2_AKI_by_landmark_sensitivity_manuscript_table.csv"
  )
)

saveRDS(
  final_results_aki_by_landmark,
  file.path(
    OUT_DIR,
    "reviewer2_AKI_by_landmark_sensitivity_full.rds"
  )
)


# ============================================================
# 10D. Optional bootstrap replicates
# ============================================================

if (SAVE_BOOT_REPLICATES) {
  
  fwrite(
    bootstrap_results,
    file.path(
      OUT_DIR,
      "reviewer2_enrichment_refinements_bootstrap_replicates.csv"
    )
  )
  
  fwrite(
    bootstrap_results_aki_by_landmark,
    file.path(
      OUT_DIR,
      "reviewer2_AKI_by_landmark_sensitivity_bootstrap_replicates.csv"
    )
  )
}


# ------------------------------------------------------------
# 11. FINAL DISPLAY
# ------------------------------------------------------------

cat("\n============================================\n")
cat("AUDIT TABLE\n")
cat("============================================\n")
print(audit_table)


cat("\n============================================\n")
cat("ORIGINAL MANUSCRIPT-READY RESULTS\n")
cat("============================================\n")
print(manuscript_table)


cat("\n============================================\n")
cat("AKI-BY-LANDMARK SENSITIVITY RESULTS\n")
cat("============================================\n")
print(manuscript_table_aki_by_landmark)


cat("\nFiles saved in:\n")
cat(
  normalizePath(OUT_DIR),
  "\n"
)
