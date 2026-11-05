#!/usr/bin/env Rscript




suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# Enrichment analysis with patient-level bootstrap
# True clinical landmarks:
#   Day 1 = 24 h, using historical pred_dt day code 0
#   Day 2 = 48 h, using historical pred_dt day code 1
# Outcome: subsequent ICU death after each landmark
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

OUT_DIR <- file.path(BASE_DIR, "Enrichment_analysis", "Enrichment_core")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

MODEL_USE <- "M2"
TARGET_CLASS <- "C2"
TARGET_PROBABILITY <- "pC2"

# IMPORTANT: historical XGBoost day code is shifted by 24 h.
LANDMARK_MAP <- data.table(
  manuscript_day = c("Day 1", "Day 2"),
  landmark_h = c(24, 48),
  pred_day_code = c(0L, 1L),
  observation_window = c("0–24 h", "0–48 h")
)

TOP_PROPORTIONS <- c(0.10, 0.20)
RRR_GRID <- c(0.10, 0.15, 0.20, 0.25, 0.30)

ALPHA <- 0.05
POWER <- 0.80

N_BOOT <- 2000L
BOOT_SEED <- 20260728L
SAVE_BOOT_REPLICATES <- FALSE

# ------------------------------------------------------------
# 1. HELPERS
# ------------------------------------------------------------

to_posix <- function(x) {
  if (inherits(x, "POSIXt")) {
    return(as.POSIXct(x, tz = "UTC"))
  }
  if (is.numeric(x)) {
    return(as.POSIXct(x, origin = "1970-01-01", tz = "UTC"))
  }
  as.POSIXct(x, tz = "UTC")
}

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

sample_size_2prop_equal_alloc_rrr <- function(
    p_control,
    rrr,
    alpha = 0.05,
    power = 0.80
) {
  p1 <- as.numeric(p_control)
  p2 <- p1 * (1 - rrr)
  
  if (
    !is.finite(p1) || !is.finite(p2) || !is.finite(rrr) ||
    p1 <= 0 || p1 >= 1 || p2 <= 0 || p2 >= 1 || p2 >= p1
  ) {
    return(NA_real_)
  }
  
  z_alpha <- qnorm(1 - alpha / 2)
  z_beta <- qnorm(power)
  p_bar <- (p1 + p2) / 2
  
  n_per_arm <- (
    z_alpha * sqrt(2 * p_bar * (1 - p_bar)) +
      z_beta * sqrt(p1 * (1 - p1) + p2 * (1 - p2))
  )^2 / (p1 - p2)^2
  
  ceiling(n_per_arm)
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

fmt_integer_ci <- function(estimate, low, high) {
  if (any(is.na(c(estimate, low, high)))) return(NA_character_)
  sprintf(
    "%s [%s–%s]",
    format(round(estimate), big.mark = ",", scientific = FALSE, trim = TRUE),
    format(round(low), big.mark = ",", scientific = FALSE, trim = TRUE),
    format(round(high), big.mark = ",", scientific = FALSE, trim = TRUE)
  )
}

get_target_label <- function(target_id) {
  switch(
    target_id,
    C2 = "JLCM-assigned Progressive AKI archetype",
    Top10 = "Top 10% by predicted probability of Progressive AKI",
    Top20 = "Top 20% by predicted probability of Progressive AKI",
    target_id
  )
}

get_population_type <- function(target_id) {
  switch(
    target_id,
    C2 = "JLCM-assigned trajectory class",
    Top10 = "Operational model-based selection",
    Top20 = "Operational model-based selection",
    NA_character_
  )
}

# ------------------------------------------------------------
# 2. METRICS WITHIN ONE SAMPLE
# ------------------------------------------------------------

compute_enrichment_metrics <- function(dt) {
  dt <- as.data.table(copy(dt))
  
  if (nrow(dt) == 0L) stop("Empty dataset supplied.")
  if (any(!is.finite(dt[[TARGET_PROBABILITY]]))) {
    stop("Non-finite target probabilities found.")
  }
  if (any(is.na(dt$outcome_after_landmark))) {
    stop("Missing outcomes found.")
  }
  
  n_unselected <- nrow(dt)
  events_unselected <- sum(dt$outcome_after_landmark == 1L)
  rate_unselected <- mean(dt$outcome_after_landmark == 1L)
  
  target_list <- list()
  target_list[["C2"]] <- dt[true_class == TARGET_CLASS]
  
  ordered_dt <- copy(dt)
  setorderv(
    ordered_dt,
    cols = TARGET_PROBABILITY,
    order = -1L,
    na.last = TRUE
  )
  
  for (top_prop in TOP_PROPORTIONS) {
    target_id <- paste0("Top", round(100 * top_prop))
    n_selected <- max(1L, ceiling(n_unselected * top_prop))
    target_list[[target_id]] <- ordered_dt[seq_len(n_selected)]
  }
  
  group_results <- list()
  sample_size_results <- list()
  
  for (target_id in names(target_list)) {
    selected <- target_list[[target_id]]
    target_label <- get_target_label(target_id)
    population_type <- get_population_type(target_id)
    
    n_targeted <- nrow(selected)
    events_targeted <- sum(selected$outcome_after_landmark == 1L)
    rate_targeted <- mean(selected$outcome_after_landmark == 1L)
    
    enrichment_factor <- if (
      is.finite(rate_unselected) && rate_unselected > 0
    ) {
      rate_targeted / rate_unselected
    } else {
      NA_real_
    }
    
    probability_threshold <- if (target_id %in% c("Top10", "Top20")) {
      min(selected[[TARGET_PROBABILITY]], na.rm = TRUE)
    } else {
      NA_real_
    }
    
    group_results[[target_id]] <- data.table(
      target_id = target_id,
      target_population = target_label,
      population_type = population_type,
      N_unselected = n_unselected,
      events_unselected = events_unselected,
      event_rate_unselected = rate_unselected,
      N_targeted = n_targeted,
      events_targeted = events_targeted,
      event_rate_targeted = rate_targeted,
      enrichment_factor = enrichment_factor,
      probability_threshold = probability_threshold
    )
    
    sample_size_results[[target_id]] <- rbindlist(
      lapply(RRR_GRID, function(rrr) {
        n_per_arm_unselected <- sample_size_2prop_equal_alloc_rrr(
          p_control = rate_unselected,
          rrr = rrr,
          alpha = ALPHA,
          power = POWER
        )
        
        n_per_arm_targeted <- sample_size_2prop_equal_alloc_rrr(
          p_control = rate_targeted,
          rrr = rrr,
          alpha = ALPHA,
          power = POWER
        )
        
        total_n_unselected <- if (is.finite(n_per_arm_unselected)) {
          2 * n_per_arm_unselected
        } else {
          NA_real_
        }
        
        total_n_targeted <- if (is.finite(n_per_arm_targeted)) {
          2 * n_per_arm_targeted
        } else {
          NA_real_
        }
        
        sample_size_reduction <- if (
          is.finite(total_n_unselected) && total_n_unselected > 0 &&
          is.finite(total_n_targeted)
        ) {
          1 - total_n_targeted / total_n_unselected
        } else {
          NA_real_
        }
        
        data.table(
          target_id = target_id,
          target_population = target_label,
          population_type = population_type,
          assumed_RRR = rrr,
          event_rate_unselected = rate_unselected,
          event_rate_targeted = rate_targeted,
          enrichment_factor = enrichment_factor,
          total_n_unselected = total_n_unselected,
          total_n_targeted = total_n_targeted,
          sample_size_reduction = sample_size_reduction
        )
      }),
      use.names = TRUE,
      fill = TRUE
    )
  }
  
  list(
    group = rbindlist(group_results, use.names = TRUE, fill = TRUE),
    sample_size = rbindlist(sample_size_results, use.names = TRUE, fill = TRUE)
  )
}

# ------------------------------------------------------------
# 3. LOAD DATA
# ------------------------------------------------------------

stopifnot(file.exists(MIMIC_BF_RDS))
stopifnot(file.exists(PRED_RDS) || file.exists(PRED_CSV))

cat("=== Loading MIMIC-IV baseline data ===\n")
bf_M <- as.data.table(readRDS(MIMIC_BF_RDS))

required_baseline_columns <- c("stay_id", "intime", "outtime", "deathtime")
missing_baseline_columns <- setdiff(required_baseline_columns, names(bf_M))
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

# Vectorised ICU-death definition
bf_M[, death_icu := as.integer(
  !is.na(deathtime) &
    !is.na(intime) &
    !is.na(outtime) &
    deathtime >= intime &
    deathtime <= outtime
)]

cat("=== Loading saved XGBoost predictions ===\n")
if (file.exists(PRED_RDS)) {
  pred_dt <- as.data.table(readRDS(PRED_RDS))
} else {
  pred_dt <- fread(PRED_CSV)
}

required_prediction_columns <- c(
  "id", "model", "day", "true_class", "pC2"
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
# 4. LANDMARK DATASETS
# ------------------------------------------------------------

landmark_datasets <- list()
landmark_diagnostics <- list()

for (i in seq_len(nrow(LANDMARK_MAP))) {
  map_i <- LANDMARK_MAP[i]
  
  pred_slice <- pred_dt[
    model == MODEL_USE & day == map_i$pred_day_code,
    .(
      stay_id = as.integer(id[1]),
      true_class = as.character(true_class[1]),
      pC2 = mean(as.numeric(pC2), na.rm = TRUE)
    ),
    by = id
  ][, id := NULL]
  
  pred_slice <- unique(pred_slice, by = "stay_id")
  pred_slice[, true_class := factor(true_class, levels = c("C1", "C2", "C3"))]
  
  dat_i <- merge(
    pred_slice,
    bf_M[, .(stay_id, intime, outtime, deathtime, death_icu)],
    by = "stay_id",
    all.x = TRUE
  )
  
  dat_i[, landmark_time := intime + map_i$landmark_h * 3600]
  
  dat_i[, eligible_at_landmark :=
          !is.na(intime) &
          !is.na(outtime) &
          outtime > landmark_time &
          (is.na(deathtime) | deathtime > landmark_time)
  ]
  
  dat_i[, outcome_after_landmark := as.integer(
    death_icu == 1L &
      !is.na(deathtime) &
      deathtime > landmark_time
  )]
  
  eligible_i <- dat_i[
    eligible_at_landmark == TRUE &
      !is.na(true_class) &
      is.finite(pC2) &
      !is.na(outcome_after_landmark)
  ]
  
  setorder(eligible_i, stay_id)
  
  if (anyDuplicated(eligible_i$stay_id) > 0L) {
    stop("Duplicate stay_id values remain at ", map_i$manuscript_day)
  }
  
  n_i <- nrow(eligible_i)
  events_i <- sum(eligible_i$outcome_after_landmark == 1L)
  rate_i <- mean(eligible_i$outcome_after_landmark == 1L)
  
  cat("\n============================================\n")
  cat(map_i$manuscript_day, "\n", sep = "")
  cat("============================================\n")
  cat("Historical pred_dt day code: ", map_i$pred_day_code, "\n", sep = "")
  cat("Observation window: ", map_i$observation_window, "\n", sep = "")
  cat("N eligible: ", n_i, "\n", sep = "")
  cat("Subsequent ICU deaths: ", events_i, "\n", sep = "")
  cat("Subsequent ICU mortality: ", sprintf("%.2f%%", 100 * rate_i), "\n", sep = "")
  
  eligible_i[, `:=`(
    manuscript_day = map_i$manuscript_day,
    landmark_h = map_i$landmark_h,
    pred_day_code = map_i$pred_day_code,
    observation_window = map_i$observation_window
  )]
  
  landmark_datasets[[map_i$manuscript_day]] <- eligible_i
  
  landmark_diagnostics[[map_i$manuscript_day]] <- data.table(
    manuscript_day = map_i$manuscript_day,
    landmark_h = map_i$landmark_h,
    pred_day_code = map_i$pred_day_code,
    observation_window = map_i$observation_window,
    n_predictions = nrow(pred_slice),
    n_eligible = n_i,
    subsequent_icu_deaths = events_i,
    subsequent_icu_mortality = rate_i
  )
}

landmark_diagnostics_dt <- rbindlist(
  landmark_diagnostics,
  use.names = TRUE,
  fill = TRUE
)

# ------------------------------------------------------------
# 5. POINT ESTIMATES + BOOTSTRAP BY LANDMARK
# ------------------------------------------------------------

all_group_summary <- list()
all_sample_size_summary <- list()
all_manuscript_table <- list()
all_bootstrap_group <- list()
all_bootstrap_sample_size <- list()

for (landmark_name in names(landmark_datasets)) {
  dat_eligible <- landmark_datasets[[landmark_name]]
  
  landmark_h_i <- unique(dat_eligible$landmark_h)
  pred_day_code_i <- unique(dat_eligible$pred_day_code)
  observation_window_i <- unique(dat_eligible$observation_window)
  
  point_results <- compute_enrichment_metrics(dat_eligible)
  point_group <- point_results$group
  point_sample_size <- point_results$sample_size
  
  set.seed(BOOT_SEED + landmark_h_i)
  
  bootstrap_group_list <- vector("list", N_BOOT)
  bootstrap_sample_size_list <- vector("list", N_BOOT)
  n_patients <- nrow(dat_eligible)
  
  cat("\nBootstrap for ", landmark_name, ": ", N_BOOT, " iterations\n", sep = "")
  
  for (bootstrap_i in seq_len(N_BOOT)) {
    bootstrap_indices <- sample.int(
      n = n_patients,
      size = n_patients,
      replace = TRUE
    )
    
    bootstrap_sample <- dat_eligible[bootstrap_indices]
    bootstrap_results <- compute_enrichment_metrics(bootstrap_sample)
    
    group_i <- bootstrap_results$group
    group_i[, bootstrap_iteration := bootstrap_i]
    bootstrap_group_list[[bootstrap_i]] <- group_i
    
    sample_size_i <- bootstrap_results$sample_size
    sample_size_i[, bootstrap_iteration := bootstrap_i]
    bootstrap_sample_size_list[[bootstrap_i]] <- sample_size_i
    
    if (bootstrap_i %% 100L == 0L) {
      cat(
        landmark_name, ": completed ",
        bootstrap_i, " of ", N_BOOT, " iterations\n",
        sep = ""
      )
    }
  }
  
  bootstrap_group <- rbindlist(
    bootstrap_group_list,
    use.names = TRUE,
    fill = TRUE
  )
  
  bootstrap_sample_size <- rbindlist(
    bootstrap_sample_size_list,
    use.names = TRUE,
    fill = TRUE
  )
  
  group_ci <- bootstrap_group[
    ,
    .(
      event_rate_unselected_low = safe_quantile(event_rate_unselected, 0.025),
      event_rate_unselected_high = safe_quantile(event_rate_unselected, 0.975),
      event_rate_targeted_low = safe_quantile(event_rate_targeted, 0.025),
      event_rate_targeted_high = safe_quantile(event_rate_targeted, 0.975),
      enrichment_factor_low = safe_quantile(enrichment_factor, 0.025),
      enrichment_factor_high = safe_quantile(enrichment_factor, 0.975),
      probability_threshold_low = safe_quantile(probability_threshold, 0.025),
      probability_threshold_high = safe_quantile(probability_threshold, 0.975)
    ),
    by = target_id
  ]
  
  group_summary <- merge(
    point_group,
    group_ci,
    by = "target_id",
    all.x = TRUE
  )
  
  group_summary[, `:=`(
    manuscript_day = landmark_name,
    landmark_h = landmark_h_i,
    pred_day_code = pred_day_code_i,
    observation_window = observation_window_i
  )]
  
  sample_size_ci <- bootstrap_sample_size[
    ,
    .(
      total_n_unselected_low = safe_quantile(total_n_unselected, 0.025),
      total_n_unselected_high = safe_quantile(total_n_unselected, 0.975),
      total_n_targeted_low = safe_quantile(total_n_targeted, 0.025),
      total_n_targeted_high = safe_quantile(total_n_targeted, 0.975),
      sample_size_reduction_low = safe_quantile(sample_size_reduction, 0.025),
      sample_size_reduction_high = safe_quantile(sample_size_reduction, 0.975)
    ),
    by = .(target_id, assumed_RRR)
  ]
  
  sample_size_summary <- merge(
    point_sample_size,
    sample_size_ci,
    by = c("target_id", "assumed_RRR"),
    all.x = TRUE
  )
  
  sample_size_summary <- merge(
    sample_size_summary,
    group_summary[
      ,
      .(
        target_id,
        N_unselected,
        events_unselected,
        N_targeted,
        events_targeted,
        event_rate_unselected_low,
        event_rate_unselected_high,
        event_rate_targeted_low,
        event_rate_targeted_high,
        enrichment_factor_low,
        enrichment_factor_high,
        probability_threshold,
        probability_threshold_low,
        probability_threshold_high
      )
    ],
    by = "target_id",
    all.x = TRUE
  )
  
  sample_size_summary[, `:=`(
    manuscript_day = landmark_name,
    landmark_h = landmark_h_i,
    pred_day_code = pred_day_code_i,
    observation_window = observation_window_i
  )]
  
  manuscript_table <- sample_size_summary[
    ,
    .(
      Cohort = "MIMIC-IV external validation",
      Landmark = manuscript_day,
      `Cumulative observation window` = observation_window,
      `Historical prediction day code` = pred_day_code,
      Model = MODEL_USE,
      `Target population` = target_population,
      `Population type` = population_type,
      `N unselected` = N_unselected,
      `Subsequent ICU deaths, unselected` = events_unselected,
      `ICU mortality, unselected (95% CI)` = fmt_pct_ci(
        event_rate_unselected,
        event_rate_unselected_low,
        event_rate_unselected_high
      ),
      `N targeted` = N_targeted,
      `Subsequent ICU deaths, targeted` = events_targeted,
      `ICU mortality, targeted (95% CI)` = fmt_pct_ci(
        event_rate_targeted,
        event_rate_targeted_low,
        event_rate_targeted_high
      ),
      `Enrichment factor (95% CI)` = fmt_num_ci(
        enrichment_factor,
        enrichment_factor_low,
        enrichment_factor_high,
        digits = 2
      ),
      `Assumed RRR` = paste0(round(100 * assumed_RRR), "%"),
      `Total N, unselected (95% CI)` = fmt_integer_ci(
        total_n_unselected,
        total_n_unselected_low,
        total_n_unselected_high
      ),
      `Total N, targeted (95% CI)` = fmt_integer_ci(
        total_n_targeted,
        total_n_targeted_low,
        total_n_targeted_high
      ),
      `Sample-size reduction (95% CI)` = fmt_pct_ci(
        sample_size_reduction,
        sample_size_reduction_low,
        sample_size_reduction_high
      )
    )
  ]
  
  all_group_summary[[landmark_name]] <- group_summary
  all_sample_size_summary[[landmark_name]] <- sample_size_summary
  all_manuscript_table[[landmark_name]] <- manuscript_table
  
  if (SAVE_BOOT_REPLICATES) {
    bootstrap_group[, manuscript_day := landmark_name]
    bootstrap_sample_size[, manuscript_day := landmark_name]
    all_bootstrap_group[[landmark_name]] <- bootstrap_group
    all_bootstrap_sample_size[[landmark_name]] <- bootstrap_sample_size
  }
}

# ------------------------------------------------------------
# 6. COMBINE AND SORT
# ------------------------------------------------------------

group_summary_all <- rbindlist(
  all_group_summary,
  use.names = TRUE,
  fill = TRUE
)

sample_size_summary_all <- rbindlist(
  all_sample_size_summary,
  use.names = TRUE,
  fill = TRUE
)

manuscript_table_all <- rbindlist(
  all_manuscript_table,
  use.names = TRUE,
  fill = TRUE
)

group_summary_all[, landmark_order := match(manuscript_day, c("Day 1", "Day 2"))]
group_summary_all[, target_order := match(target_id, c("C2", "Top10", "Top20"))]
setorder(group_summary_all, landmark_order, target_order)
group_summary_all[, c("landmark_order", "target_order") := NULL]

sample_size_summary_all[, landmark_order := match(manuscript_day, c("Day 1", "Day 2"))]
sample_size_summary_all[, target_order := match(target_id, c("C2", "Top10", "Top20"))]
setorder(sample_size_summary_all, landmark_order, target_order, assumed_RRR)
sample_size_summary_all[, c("landmark_order", "target_order") := NULL]

manuscript_table_all[, landmark_order := match(Landmark, c("Day 1", "Day 2"))]
manuscript_table_all[, target_order := match(
  `Target population`,
  c(
    "JLCM-assigned Progressive AKI archetype",
    "Top 10% by predicted probability of Progressive AKI",
    "Top 20% by predicted probability of Progressive AKI"
  )
)]
manuscript_table_all[, rrr_order := as.numeric(sub("%", "", `Assumed RRR`, fixed = TRUE))]
setorder(manuscript_table_all, landmark_order, target_order, rrr_order)
manuscript_table_all[, c("landmark_order", "target_order", "rrr_order") := NULL]

# Compact observed-enrichment panel: one row per landmark and strategy
observed_enrichment_table <- group_summary_all[
  ,
  .(
    Cohort = "MIMIC-IV external validation",
    Landmark = manuscript_day,
    `Cumulative observation window` = observation_window,
    `Historical prediction day code` = pred_day_code,
    `Target population` = target_population,
    `N unselected` = N_unselected,
    `Subsequent ICU deaths, unselected` = events_unselected,
    `ICU mortality, unselected (95% CI)` = fmt_pct_ci(
      event_rate_unselected,
      event_rate_unselected_low,
      event_rate_unselected_high
    ),
    `N targeted` = N_targeted,
    `Subsequent ICU deaths, targeted` = events_targeted,
    `ICU mortality, targeted (95% CI)` = fmt_pct_ci(
      event_rate_targeted,
      event_rate_targeted_low,
      event_rate_targeted_high
    ),
    `Enrichment factor (95% CI)` = fmt_num_ci(
      enrichment_factor,
      enrichment_factor_low,
      enrichment_factor_high,
      digits = 2
    )
  )
]

# ------------------------------------------------------------
# 7. EXPORTS
# ------------------------------------------------------------

fwrite(
  landmark_diagnostics_dt,
  file.path(OUT_DIR, "A8_2_landmark_diagnostics_day1_day2.csv")
)

fwrite(
  observed_enrichment_table,
  file.path(OUT_DIR, "TableS13_panelA_observed_enrichment_day1_day2.csv")
)

fwrite(
  group_summary_all,
  file.path(OUT_DIR, "A8_2_enrichment_mortality_factor_day1_day2_with_CI.csv")
)

fwrite(
  sample_size_summary_all,
  file.path(OUT_DIR, "A8_2_sample_size_reduction_day1_day2_with_CI.csv")
)

fwrite(
  manuscript_table_all,
  file.path(OUT_DIR, "TableS13_enrichment_bootstrap_day1_day2.csv")
)

saveRDS(
  observed_enrichment_table,
  file.path(OUT_DIR, "TableS13_panelA_observed_enrichment_day1_day2.rds")
)

saveRDS(
  sample_size_summary_all,
  file.path(OUT_DIR, "A8_2_sample_size_reduction_day1_day2_with_CI.rds")
)

saveRDS(
  manuscript_table_all,
  file.path(OUT_DIR, "TableS13_enrichment_bootstrap_day1_day2.rds")
)

if (SAVE_BOOT_REPLICATES) {
  bootstrap_group_all <- rbindlist(
    all_bootstrap_group,
    use.names = TRUE,
    fill = TRUE
  )
  bootstrap_sample_size_all <- rbindlist(
    all_bootstrap_sample_size,
    use.names = TRUE,
    fill = TRUE
  )
  
  fwrite(
    bootstrap_group_all,
    file.path(OUT_DIR, "bootstrap_group_replicates_day1_day2.csv")
  )
  fwrite(
    bootstrap_sample_size_all,
    file.path(OUT_DIR, "bootstrap_sample_size_replicates_day1_day2.csv")
  )
}

# ------------------------------------------------------------
# 8. FINAL DISPLAY
# ------------------------------------------------------------

cat("\n============================================\n")
cat("LANDMARK DIAGNOSTICS\n")
cat("============================================\n")
print(landmark_diagnostics_dt)

cat("\n============================================\n")
cat("OBSERVED ENRICHMENT — DAY 1 AND DAY 2\n")
cat("============================================\n")
print(observed_enrichment_table)
