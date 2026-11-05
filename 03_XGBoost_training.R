# Historical source pipeline, submitted by investigators. NOT YET VALIDATED as the
# final landmark-aligned analysis; risk-set filtering must be checked.
#!/usr/bin/env Rscript
rm(list=ls()); gc(); graphics.off()

# À ajouter au début du script
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(xgboost)
  library(Matrix)
  library(pROC)
  library(ggplot2)
  library(parallel)      # NOUVEAU
  library(doParallel)    # NOUVEAU
  library(foreach)       # NOUVEAU
})

set.seed(20260130)

# ============================================================
# 0) PATHS + PARAMS
# ============================================================

BASE_DIR <- Sys.getenv("AKI_BASE_DIR")
if (!nzchar(BASE_DIR) || !dir.exists(BASE_DIR)) {
  stop("Set AKI_BASE_DIR to the local, access-controlled directory containing DF/ and Modele Joint/")
}

# Data
AUMC_BF_RDS <- file.path(BASE_DIR, "DF/AUMC/aumc_bf.rds")
AUMC_TS_RDS <- file.path(BASE_DIR, "DF/AUMC/aumc_ts.rds")

MIMIC_BF_RDS <- file.path(BASE_DIR, "DF/MIMIC/mimic31_aki_admission_base_v3.rds")
MIMIC_TS_RDS <- file.path(BASE_DIR, "DF/MIMIC/timeseries_0_14d_core_long.rds")

# Labels
AUMC_MODEL_RDS <- file.path(BASE_DIR, "Modele Joint/AUMC/AKISCr_J14/AUMC_joint_G3_AUMC_SCr_H14_quad_w99_cr_H14d_time_s.rds")
MIMIC_PP_CSV   <- file.path(BASE_DIR, "Modele Joint/AUMC/AKISCr_J14/external_validation_MIMIC/tables/ppMIMIC_AUMC_joint_G3_AUMC_SCr_H14_quad_w99_cr_H14d_time_s.csv")

# Output
OUT_DIR <- Sys.getenv("AKI_XGB_OUT_DIR", unset = file.path(BASE_DIR, "XGB_external_AUMC_train__MIMIC_G3"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Landmark days
DAYS_TARGET <- 0:7
MAX_DAY     <- 14L

# XGB params (Dominic-like)
XGB_MC <- list(
  booster          = "gbtree",
  objective        = "multi:softprob",
  eval_metric      = "mlogloss",
  eta              = 0.05,
  max_depth        = 4,
  subsample        = 0.7,
  colsample_bytree = 1.0
)

XGB_BIN <- list(
  booster          = "gbtree",
  objective        = "binary:logistic",
  eval_metric      = "auc",
  eta              = 0.05,
  max_depth        = 4,
  subsample        = 0.7,
  colsample_bytree = 1.0
)

NROUND_MAX  <- 300
NROUND_ESTOP <- 20
NFOLD <- 5

# Optional: recalibration MIMIC 80/20
DO_RECAL_MIMIC <- TRUE
RECAL_PROP <- 0.80

msg <- function(...) cat(sprintf(...), "\n")

# ============================================================
# 1) SMALL HELPERS
# ============================================================

detect_value_col <- function(dt) {
  if ("value" %in% names(dt)) return("value")
  if ("med"   %in% names(dt)) return("med")
  stop("No value column found (expected 'value' or 'med').")
}

normalize_gender <- function(x) {
  x <- as.character(x)
  x <- ifelse(x %in% c("M","Male","MALE"), "Male",
              ifelse(x %in% c("F","Female","FEMALE"), "Female",
                     ifelse(is.na(x) | x %in% c("U","Unknown","UNK",""), "U", x)))
  x
}

# Hard leakage guard
assert_no_leakage <- function(nms) {
  bad_patterns <- c("los", "expire", "death", "discharge", "event_", "duration",
                    "delay", "onset", "aki_stage", "tte_", "outcome")
  bad <- nms[grepl(paste(bad_patterns, collapse="|"), nms, ignore.case=TRUE)]
  if (length(bad) > 0) stop("Leakage-like cols detected: ", paste(bad, collapse=", "))
  invisible(TRUE)
}

# Softmax safety (xgb already outputs probs for multi:softprob, but keep safe)
safe_prob <- function(p, eps = 1e-15) pmin(pmax(p, eps), 1-eps)

# Multiclass balanced accuracy = mean recall per class
balanced_acc <- function(y_true_int, y_pred_int, K) {
  rec <- sapply(0:(K-1), function(k) {
    den <- sum(y_true_int == k)
    if (den == 0) return(NA_real_)
    sum(y_true_int == k & y_pred_int == k) / den
  })
  mean(rec, na.rm = TRUE)
}

# OvR AUC macro/weighted from multiclass probability matrix
one_vs_rest_auc <- function(y_true_int, p_mat, k) {
  y_bin <- as.integer(y_true_int == k)
  if (length(unique(y_bin)) < 2) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(response = y_bin, predictor = p_mat[, k+1], quiet = TRUE)))
}

multiclass_auc_ovr <- function(y_true_int, p_mat, K) {
  auc_k <- sapply(0:(K-1), function(k) one_vs_rest_auc(y_true_int, p_mat, k))
  prev  <- sapply(0:(K-1), function(k) mean(y_true_int == k))
  list(
    auc_by_class = auc_k,
    auc_macro    = mean(auc_k, na.rm = TRUE),
    auc_weighted = sum(auc_k * prev, na.rm = TRUE)
  )
}

# Multiclass Brier (one-hot)
brier_multiclass <- function(y_true_int, p_mat, K) {
  n <- length(y_true_int)
  Y <- matrix(0, nrow = n, ncol = K)
  Y[cbind(1:n, y_true_int + 1L)] <- 1
  mean(rowSums((p_mat - Y)^2))
}

# ECE for binary probs
# ECE pour binaire - CORRIGÉ
ece_binary <- function(prob, truth01, k = 10) {
  prob <- safe_prob(prob)
  cuts <- seq(0, 1, length.out = k + 1)
  bin  <- cut(prob, cuts, include.lowest = TRUE, labels = FALSE)
  
  # Calculer les statistiques par bin
  exp  <- tapply(prob, bin, mean, na.rm = TRUE)
  obs  <- tapply(truth01, bin, mean, na.rm = TRUE)
  nbin <- tapply(truth01, bin, length)
  
  # Remplacer les NA par 0
  exp[is.na(exp)] <- 0
  obs[is.na(obs)] <- 0
  nbin[is.na(nbin)] <- 0
  
  # Calculer ECE
  ece <- if (sum(nbin) > 0) sum(abs(obs-exp) * nbin) / sum(nbin) else NA_real_
  
  # Créer le data frame SEULEMENT pour les bins non vides
  bins_present <- as.integer(names(exp))
  df <- data.frame(
    bin = bins_present,
    exp = as.numeric(exp),
    obs = as.numeric(obs),
    n = as.numeric(nbin)
  )
  
  list(ece = ece, df = df)
}

# ECE multiclasse poolé - CORRIGÉ
ece_multiclass_pooled <- function(p_mat, y_true, k = 10) {
  # Pooler toutes les prédictions
  probs_pooled <- as.vector(p_mat)
  
  # Créer les vraies valeurs binaires correspondantes
  n_samples <- length(y_true)
  n_classes <- ncol(p_mat)
  y_binary <- rep(0, n_samples * n_classes)
  
  for (i in 1:n_samples) {
    for (j in 1:n_classes) {
      idx <- (i - 1) * n_classes + j
      y_binary[idx] <- as.integer(y_true[i] == (j - 1))
    }
  }
  
  # Créer les bins
  breaks <- seq(0, 1, length.out = k + 1)
  bins <- cut(probs_pooled, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  
  # Calculer les statistiques par bin en utilisant tapply
  exp  <- tapply(probs_pooled, bins, mean, na.rm = TRUE)
  obs  <- tapply(y_binary, bins, mean, na.rm = TRUE)
  nbin <- tapply(y_binary, bins, length)
  
  # Remplacer les NA par 0
  exp[is.na(exp)] <- 0
  obs[is.na(obs)] <- 0
  nbin[is.na(nbin)] <- 0
  
  # Calculer l'ECE
  ece <- if (sum(nbin) > 0) sum(abs(exp - obs) * nbin) / sum(nbin) else NA_real_
  
  # Créer le data frame SEULEMENT pour les bins présents
  bins_present <- as.integer(names(exp))
  df <- data.frame(
    bin = bins_present,
    exp = as.numeric(exp),
    obs = as.numeric(obs),
    n = as.numeric(nbin)
  )
  
  list(ece = ece, df = df)
}

save_reliability <- function(cal_df, title, path_png) {
  p <- ggplot(cal_df %>% dplyr::filter(n > 0),
              aes(x = exp, y = obs, size = n)) +
    geom_point(alpha = 0.7) +
    geom_abline(linetype = "dashed") +
    coord_equal(xlim = c(0,1), ylim = c(0,1)) +
    theme_bw(base_size = 12) +
    labs(title = title, x = "Predicted probability", y = "Observed frequency", size = "n")
  ggsave(path_png, p, width = 5, height = 5, dpi = 200)
}

# ============================================================
# BOOTSTRAP CONFIDENCE INTERVALS - VERSION PARALLÈLE
# ============================================================

# Configuration parallèle globale (à appeler une seule fois)
setup_parallel <- function(n_cores = NULL) {
  if (is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)  # Garder 1 core libre
  }
  msg("Setting up parallel processing with %d cores", n_cores)
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  return(cl)
}

stop_parallel <- function(cl) {
  if (!is.null(cl)) {
    stopCluster(cl)
    msg("Parallel processing stopped")
  }
}

# Bootstrap pour balanced accuracy (parallélisé)
bootstrap_balanced_acc <- function(y_true_int, y_pred_int, K, n_boot = 1000, seed = 123) {
  set.seed(seed)
  n <- length(y_true_int)
  
  # Générer tous les seeds à l'avance pour reproductibilité
  boot_seeds <- sample.int(.Machine$integer.max, n_boot)
  
  boot_stats <- foreach(i = 1:n_boot, .combine = c, 
                        .packages = NULL,
                        .export = c("balanced_acc")) %dopar% {
                          set.seed(boot_seeds[i])
                          idx <- sample(1:n, replace = TRUE)
                          balanced_acc(y_true_int[idx], y_pred_int[idx], K)
                        }
  
  list(
    mean = mean(boot_stats, na.rm = TRUE),
    ci_lower = quantile(boot_stats, 0.025, na.rm = TRUE),
    ci_upper = quantile(boot_stats, 0.975, na.rm = TRUE),
    se = sd(boot_stats, na.rm = TRUE)
  )
}

# Bootstrap pour AUC binaire (parallélisé)
bootstrap_auc_binary <- function(y_true_bin, prob, n_boot = 1000, seed = 123) {
  if (length(unique(y_true_bin)) < 2) {
    return(list(mean = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, se = NA_real_))
  }
  
  set.seed(seed)
  n <- length(y_true_bin)
  boot_seeds <- sample.int(.Machine$integer.max, n_boot)
  
  boot_stats <- foreach(i = 1:n_boot, .combine = c,
                        .packages = c("pROC"),
                        .export = NULL) %dopar% {
                          set.seed(boot_seeds[i])
                          idx <- sample(1:n, replace = TRUE)
                          y_boot <- y_true_bin[idx]
                          p_boot <- prob[idx]
                          
                          if (length(unique(y_boot)) < 2) return(NA_real_)
                          
                          roc_obj <- tryCatch(
                            pROC::roc(y_boot, p_boot, quiet = TRUE),
                            error = function(e) NULL
                          )
                          if (is.null(roc_obj)) return(NA_real_)
                          as.numeric(pROC::auc(roc_obj))
                        }
  
  list(
    mean = mean(boot_stats, na.rm = TRUE),
    ci_lower = quantile(boot_stats, 0.025, na.rm = TRUE),
    ci_upper = quantile(boot_stats, 0.975, na.rm = TRUE),
    se = sd(boot_stats, na.rm = TRUE)
  )
}

# Bootstrap pour AUC multiclasse (parallélisé)
bootstrap_auc_macro <- function(y_true_int, p_mat, K, n_boot = 1000, seed = 123) {
  set.seed(seed)
  n <- length(y_true_int)
  boot_seeds <- sample.int(.Machine$integer.max, n_boot)
  
  boot_stats <- foreach(i = 1:n_boot, .combine = c,
                        .packages = c("pROC"),
                        .export = c("multiclass_auc_ovr", "one_vs_rest_auc")) %dopar% {
                          set.seed(boot_seeds[i])
                          idx <- sample(1:n, replace = TRUE)
                          y_boot <- y_true_int[idx]
                          p_boot <- p_mat[idx, , drop = FALSE]
                          
                          auc_res <- multiclass_auc_ovr(y_boot, p_boot, K)
                          auc_res$auc_macro
                        }
  
  list(
    mean = mean(boot_stats, na.rm = TRUE),
    ci_lower = quantile(boot_stats, 0.025, na.rm = TRUE),
    ci_upper = quantile(boot_stats, 0.975, na.rm = TRUE),
    se = sd(boot_stats, na.rm = TRUE)
  )
}

# Bootstrap pour Brier score (parallélisé)
bootstrap_brier_mc <- function(y_true_int, p_mat, K, n_boot = 1000, seed = 123) {
  set.seed(seed)
  n <- length(y_true_int)
  boot_seeds <- sample.int(.Machine$integer.max, n_boot)
  
  boot_stats <- foreach(i = 1:n_boot, .combine = c,
                        .packages = NULL,
                        .export = c("brier_multiclass")) %dopar% {
                          set.seed(boot_seeds[i])
                          idx <- sample(1:n, replace = TRUE)
                          y_boot <- y_true_int[idx]
                          p_boot <- p_mat[idx, , drop = FALSE]
                          brier_multiclass(y_boot, p_boot, K)
                        }
  
  list(
    mean = mean(boot_stats, na.rm = TRUE),
    ci_lower = quantile(boot_stats, 0.025, na.rm = TRUE),
    ci_upper = quantile(boot_stats, 0.975, na.rm = TRUE),
    se = sd(boot_stats, na.rm = TRUE)
  )
}

# Bootstrap pour ECE (parallélisé)
bootstrap_ece <- function(y_true_int, p_mat, K, n_boot = 1000, seed = 123, k_bins = 10) {
  set.seed(seed)
  n <- length(y_true_int)
  boot_seeds <- sample.int(.Machine$integer.max, n_boot)
  
  boot_stats <- foreach(i = 1:n_boot, .combine = c,
                        .packages = NULL,
                        .export = c("ece_multiclass_pooled", "safe_prob")) %dopar% {
                          set.seed(boot_seeds[i])
                          idx <- sample(1:n, replace = TRUE)
                          y_boot <- y_true_int[idx]
                          p_boot <- p_mat[idx, , drop = FALSE]
                          
                          cal <- ece_multiclass_pooled(p_boot, y_boot, k = k_bins)
                          cal$ece
                        }
  
  list(
    mean = mean(boot_stats, na.rm = TRUE),
    ci_lower = quantile(boot_stats, 0.025, na.rm = TRUE),
    ci_upper = quantile(boot_stats, 0.975, na.rm = TRUE),
    se = sd(boot_stats, na.rm = TRUE)
  )
}
# ============================================================
# DECISION CURVE ANALYSIS
# ============================================================

# DCA pour classification binaire
dca_binary <- function(y_true_bin, prob, 
                       thresholds = seq(0.01, 0.99, by = 0.01),
                       label = "Model") {
  
  n <- length(y_true_bin)
  prevalence <- mean(y_true_bin)
  
  # Calculer le net benefit pour chaque seuil
  nb_model <- sapply(thresholds, function(pt) {
    # Décision basée sur le modèle
    pred <- as.integer(prob >= pt)
    
    # True positives et False positives
    tp <- sum(pred == 1 & y_true_bin == 1)
    fp <- sum(pred == 1 & y_true_bin == 0)
    
    # Net benefit
    (tp / n) - (fp / n) * (pt / (1 - pt))
  })
  
  # Strategy "Treat all"
  nb_all <- sapply(thresholds, function(pt) {
    prevalence - (1 - prevalence) * (pt / (1 - pt))
  })
  
  # Strategy "Treat none"
  nb_none <- rep(0, length(thresholds))
  
  # Data frame pour ggplot
  df <- data.frame(
    threshold = rep(thresholds, 3),
    net_benefit = c(nb_model, nb_all, nb_none),
    strategy = rep(c(label, "Treat All", "Treat None"), each = length(thresholds))
  )
  
  list(df = df, prevalence = prevalence)
}

# DCA pour classification multiclasse (One-vs-Rest)
# Pour chaque classe, on calcule le net benefit de prédire cette classe
dca_multiclass_ovr <- function(y_true_int, p_mat, K, 
                               thresholds = seq(0.01, 0.99, by = 0.01),
                               class_labels = NULL) {
  
  if (is.null(class_labels)) class_labels <- paste0("C", 1:K)
  
  dca_results <- lapply(0:(K-1), function(k) {
    y_bin <- as.integer(y_true_int == k)
    prob_k <- p_mat[, k + 1]
    
    dca_res <- dca_binary(y_bin, prob_k, thresholds, label = class_labels[k + 1])
    dca_res$df$class <- class_labels[k + 1]
    dca_res$df
  })
  
  rbindlist(dca_results)
}

# Plot DCA (binaire)
plot_dca_binary <- function(dca_df, title = "Decision Curve Analysis") {
  ggplot(dca_df, aes(x = threshold, y = net_benefit, color = strategy, linetype = strategy)) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_color_manual(values = c("Model" = "#E69F00", "Treat All" = "#56B4E9", "Treat None" = "#999999")) +
    scale_linetype_manual(values = c("Model" = "solid", "Treat All" = "dashed", "Treat None" = "dotted")) +
    theme_bw(base_size = 12) +
    labs(
      title = title,
      x = "Threshold Probability",
      y = "Net Benefit",
      color = "Strategy",
      linetype = "Strategy"
    ) +
    theme(legend.position = "bottom") +
    xlim(0, 1)
}

# Plot DCA multiclasse (facettes par classe)
plot_dca_multiclass <- function(dca_mc_df, title = "Decision Curve Analysis (Multiclass)") {
  # Filtrer seulement les modèles (pas "Treat All" / "Treat None" pour lisibilité)
  model_only <- dca_mc_df %>% filter(strategy != "Treat All" & strategy != "Treat None")
  
  ggplot(dca_mc_df, aes(x = threshold, y = net_benefit, color = strategy, linetype = strategy)) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    facet_wrap(~ class, ncol = 2, scales = "free_y") +
    theme_bw(base_size = 10) +
    labs(
      title = title,
      x = "Threshold Probability",
      y = "Net Benefit",
      color = "Strategy",
      linetype = "Strategy"
    ) +
    theme(legend.position = "bottom") +
    xlim(0, 1)
}

# ============================================================
# 2) LOAD LABELS
# ============================================================

msg("== Load AUMC model labels ==")
mod_joint_A <- readRDS(AUMC_MODEL_RDS)
pp_A <- as.data.frame(mod_joint_A$pprob) %>%
  mutate(stay_id = id,
         class = factor(class, levels = 1:3, labels = paste0("C", 1:3))) %>%
  select(stay_id, class)

msg("== Load MIMIC labels (CSV) ==")
pp_raw <- suppressMessages(readr::read_csv(MIMIC_PP_CSV, show_col_types = FALSE)) %>% as.data.frame()
if ("stay_id" %in% names(pp_raw)) {
  pp_M <- pp_raw %>% mutate(stay_id = as.integer(stay_id))
} else if ("id" %in% names(pp_raw)) {
  pp_M <- pp_raw %>% mutate(stay_id = as.integer(id))
} else {
  stop("MIMIC pp CSV must contain 'stay_id' or 'id'.")
}

pp_M <- pp_M %>%
  mutate(class = factor(class, levels = 1:3, labels = paste0("C", 1:3))) %>%
  select(stay_id, class)

class_levels <- levels(pp_A$class)
K <- length(class_levels)

# ============================================================
# 3) LOAD DATA (BF + TS)
# ============================================================

msg("== Load BF/TS ==")
bf_A <- readRDS(AUMC_BF_RDS) %>% as.data.frame()
ts_A <- readRDS(AUMC_TS_RDS) %>% as.data.frame()

bf_M <- readRDS(MIMIC_BF_RDS) %>% as.data.frame()
ts_M <- readRDS(MIMIC_TS_RDS) %>% as.data.frame()

setDT(ts_A); setDT(ts_M)

bf_cls_A <- bf_A %>% inner_join(pp_A, by = "stay_id") %>% mutate(cohort = "AUMC")
bf_cls_M <- bf_M %>% inner_join(pp_M, by = "stay_id") %>% mutate(cohort = "MIMIC")

msg("AUMC labelled stays: %d", nrow(bf_cls_A))
msg("MIMIC labelled stays: %d", nrow(bf_cls_M))

# ============================================================
# 4) HARMONISE BF (minimal)
# ============================================================

bf_cls_A <- bf_cls_A %>%
  mutate(
    sepsis_adm = as.integer((sepsis_icu %in% 1) | (group_sepsis_adm %in% 1) | (sepsis_at_admission %in% 1)),
    surgery_adm = as.integer(grepl("SURG", admission_type, ignore.case = TRUE)),
    exp_acei_arb_72h = as.integer((exp_acei_72h %in% 1) | (exp_arb_72h %in% 1)),
    exp_nephrotox_72h = as.integer(exp_atb_nephrotox_72h %in% 1),
    exp_abx_72h = as.integer((exp_vanco %in% 1) | (exp_genta %in% 1) | (exp_tobra %in% 1) |
                               (exp_amika %in% 1) | (exp_bactrim %in% 1) | (exp_atb_nephrotox_72h %in% 1))
  )

bf_cls_M <- bf_cls_M %>%
  mutate(
    sepsis_adm = as.integer(group_sepsis_adm %in% 1),
    surgery_adm = as.integer(grepl("SURG", admission_type, ignore.case = TRUE)),
    exp_acei_arb_72h = as.integer((exp_acei_72h %in% 1) | (exp_arb_72h %in% 1)),
    exp_nephrotox_72h = as.integer(exp_atb_nephrotox_72h %in% 1),
    exp_abx_72h = as.integer(exp_atb_72h %in% 1)
  )

bf_cls_A <- bf_cls_A %>% mutate(ckd = as.integer(ckd))
bf_cls_M <- bf_cls_M %>% mutate(ckd = as.integer(ckd))

bf_cls_A$gender_std <- normalize_gender(bf_cls_A$gender_std)
bf_cls_M$gender_std <- normalize_gender(bf_cls_M$gender_std)

bf_vars_keep <- c(
  "agegroup","gender_std","bmi_admit",
  "egfr_baseline_mlmin","ckd",
  "admission_type", "sepsis_adm", "surgery_adm",
  "sofa_total_adm","sofa_resp_adm","sofa_cardio_adm","sofa_renal_adm",
  "sofa_cns_adm","sofa_coag_adm","sofa_liver_adm",
  "exp_abx_72h","exp_nsaid_72h","exp_nephrotox_72h","exp_acei_arb_72h"
)

missA <- setdiff(bf_vars_keep, names(bf_cls_A))
missM <- setdiff(bf_vars_keep, names(bf_cls_M))
if (length(missA) > 0) stop("Missing BF vars in AUMC: ", paste(missA, collapse=", "))
if (length(missM) > 0) stop("Missing BF vars in MIMIC: ", paste(missM, collapse=", "))

# ============================================================
# 5) TS MAP + RRT derivation (same as your pipeline)
# ============================================================

ts_map <- data.table(
  canon = c(
    "scr","bun","uo_mlkghr_6hr",
    "hr","sbp","dbp","map","temp",
    "fio2","pao2","pf_ratio","peep","vt_observed","vt_set",
    "flag_any_vaso","rate_norepinephrine","rate_epinephrine","rate_phenylephrine",
    "na","k","lactate","plt","wbc",
    "flag_vm_invasive"
  ),
  aumc = c(
    "scr","bun","uo_mlkghr_6hr",
    "hr","sbp","dbp","map","temp",
    "fio2","pao2","pao2_fio2","peep","vt_observed","vt_set",
    "flag_any_vaso","rate_norepinephrine","rate_epinephrine","rate_phenylephrine",
    "na","k","lactate_bg","plt","wbc",
    "flag_vm_invasive"
  ),
  mimic = c(
    "scr","bun","uo_mlkghr_6hr",
    "hr","sbp","dbp","map","temp",
    "fio2","pao2","pao2_fio2","peep","vt_observed","vt_set",
    "flag_any_vaso","rate_norepinephrine","rate_epinephrine","rate_phenylephrine",
    "na","k","lactate_bg","plt","wbc",
    "flag_vm_invasive"
  )
)

canonize_ts <- function(ts_dt, map_dt, side = c("aumc","mimic")) {
  side <- match.arg(side)
  m <- copy(map_dt)[, .(canon, src = get(side))]
  setkey(m, src)
  ts_dt <- ts_dt[var %in% m$src]
  setkey(ts_dt, var)
  ts_dt <- m[ts_dt, on = .(src = var), nomatch = 0L]
  ts_dt[, var := canon][, c("canon","src") := NULL]
  ts_dt
}

aggregate_daily <- function(ts_dt, val_col, max_day = 14L) {
  ts_dt[
    ,
    .(
      mean = mean(get(val_col), na.rm = TRUE),
      last = {
        v <- get(val_col); i <- which(!is.na(v))
        if (length(i) == 0) NA_real_ else v[max(i)]
      }
    ),
    by = .(stay_id, var, t_day = pmin(max_day, floor(pmax(0, t_hours)/24L)))
  ]
}

pivot_daily_wide <- function(daily_dt) {
  m1 <- dcast(daily_dt, stay_id + t_day ~ paste0(var, "_mean"), value.var = "mean")
  m2 <- dcast(daily_dt, stay_id + t_day ~ paste0(var, "_last"), value.var = "last")
  setkey(m1, stay_id, t_day); setkey(m2, stay_id, t_day)
  m1[m2]
}

collapse_window <- function(daily_wide, k) {
  dt <- daily_wide[t_day <= k]
  feat_cols <- setdiff(names(dt), c("stay_id","t_day"))
  win_mean <- dt[, lapply(.SD, function(x) mean(x, na.rm = TRUE)),
                 by = stay_id, .SDcols = feat_cols]
  setnames(win_mean, feat_cols, paste0(feat_cols, "_wmean_d0_", k))
  dt_k <- daily_wide[t_day == k]
  win_last <- if (nrow(dt_k) == 0) data.table(stay_id = unique(daily_wide$stay_id))
  else dt_k[, c("stay_id", feat_cols), with=FALSE]
  setnames(win_last, feat_cols, paste0(feat_cols, "_d", k))
  setkey(win_mean, stay_id); setkey(win_last, stay_id)
  win_mean[win_last]
}

prepare_dataset <- function(bf_cls, ts_dt, ts_map, side,
                            bf_vars_keep, days_target,
                            rrt_day_full, max_day = 14L) {
  bf_dt <- as.data.table(bf_cls)[, c("stay_id","class","cohort", bf_vars_keep), with=FALSE]
  setkey(bf_dt, stay_id)
  
  ts_dt <- copy(ts_dt)
  val_col <- detect_value_col(ts_dt)
  
  ts_dt <- ts_dt[var %in% ts_map[[side]]]
  ts_dt <- canonize_ts(ts_dt, as.data.table(ts_map), side)
  
  daily_long <- aggregate_daily(ts_dt, val_col, max_day)
  daily_wide <- pivot_daily_wide(daily_long)
  
  setkey(daily_wide, stay_id, t_day)
  setkey(rrt_day_full, stay_id, t_day)
  daily_wide <- daily_wide[rrt_day_full]
  daily_wide[is.na(flag_rrt), flag_rrt := 0L]
  
  out <- rbindlist(lapply(days_target, function(k) {
    Xk <- collapse_window(daily_wide, k)
    setkey(Xk, stay_id)
    ds <- bf_dt[Xk]
    ds[, day_window := as.integer(k)]
    ds
  }), use.names=TRUE, fill=TRUE)
  
  out
}

# ---- Build RRT grids
ids_A <- unique(bf_cls_A$stay_id)
ids_M <- unique(bf_cls_M$stay_id)

vars_A_need <- unique(c(ts_map$aumc, "flag_rrt"))
vars_M_need <- unique(c(ts_map$mimic, "crrt_uf_ml_hour"))

ts_cols_keep_A <- intersect(c("stay_id","t_hours","var","value","med"), names(ts_A))
ts_cols_keep_M <- intersect(c("stay_id","t_hours","var","value","med"), names(ts_M))

ts_A_s <- ts_A[stay_id %in% ids_A & var %in% vars_A_need, ..ts_cols_keep_A]
ts_M_s <- ts_M[stay_id %in% ids_M & var %in% vars_M_need, ..ts_cols_keep_M]

# MIMIC RRT flag
rrt_day_M <- ts_M_s[
  var == "crrt_uf_ml_hour",
  .(flag_rrt = 1L),
  by = .(stay_id, t_day = pmin(MAX_DAY, floor(pmax(0, t_hours)/24L)))
]
rrt_day_M_full <- CJ(stay_id = ids_M, t_day = 0:MAX_DAY, unique = TRUE)
setkey(rrt_day_M_full, stay_id, t_day)
setkey(rrt_day_M, stay_id, t_day)
rrt_day_M_full[rrt_day_M, flag_rrt := i.flag_rrt]
rrt_day_M_full[is.na(flag_rrt), flag_rrt := 0L]

# AUMC RRT flag
val_col_A <- detect_value_col(ts_A_s)
rrt_day_A <- ts_A_s[
  var == "flag_rrt",
  .(flag_rrt = as.integer(any(get(val_col_A) == 1, na.rm = TRUE))),
  by = .(stay_id, t_day = pmin(MAX_DAY, floor(pmax(0, t_hours)/24L)))
]
rrt_day_A_full <- CJ(stay_id = ids_A, t_day = 0:MAX_DAY, unique = TRUE)
setkey(rrt_day_A_full, stay_id, t_day)
setkey(rrt_day_A, stay_id, t_day)
rrt_day_A_full[rrt_day_A, flag_rrt := i.flag_rrt]
rrt_day_A_full[is.na(flag_rrt), flag_rrt := 0L]

# ============================================================
# 6) BUILD DATASETS (AUMC train | MIMIC test) per day_window
# ============================================================

msg("== Build AUMC dataset ==")
ds_A <- prepare_dataset(
  bf_cls = bf_cls_A,
  ts_dt  = ts_A_s[var != "flag_rrt"],
  ts_map = ts_map,
  side   = "aumc",
  bf_vars_keep = bf_vars_keep,
  days_target = DAYS_TARGET,
  rrt_day_full = rrt_day_A_full,
  max_day = MAX_DAY
)
msg("AUMC DS: %d rows | %d cols", nrow(ds_A), ncol(ds_A))

msg("== Build MIMIC dataset ==")
ds_M <- prepare_dataset(
  bf_cls = bf_cls_M,
  ts_dt  = ts_M_s[var != "crrt_uf_ml_hour"],
  ts_map = ts_map,
  side   = "mimic",
  bf_vars_keep = bf_vars_keep,
  days_target = DAYS_TARGET,
  rrt_day_full = rrt_day_M_full,
  max_day = MAX_DAY
)
msg("MIMIC DS: %d rows | %d cols", nrow(ds_M), ncol(ds_M))

assert_no_leakage(names(ds_A))
assert_no_leakage(names(ds_M))

# ============================================================
# 7) DEFINE FEATURE SETS (M0..M4 + S0_SCR)
# ============================================================

setDT(ds_A); setDT(ds_M)

# Factor handling: learn levels from AUMC ONLY, apply to MIMIC
for (cc in c("agegroup","gender_std","admission_type")) {
  ds_A[is.na(get(cc)), (cc) := "UNK"]
  ds_M[is.na(get(cc)), (cc) := "UNK"]
  ds_A[, (cc) := factor(as.character(get(cc)))]
  lev <- levels(ds_A[[cc]])
  ds_M[, (cc) := factor(as.character(get(cc)), levels = lev)]
  ds_M[is.na(get(cc)), (cc) := "UNK"]  # fallback (if "UNK" not in lev, you'll see NAs; add UNK in AUMC if needed)
}

levels_y <- class_levels

fixed_vars <- c(
  "agegroup","gender_std","admission_type",
  "bmi_admit","egfr_baseline_mlmin","ckd",
  "sepsis_adm","surgery_adm","sofa_total_adm"
)

renal_ts <- grep("^scr_|^uo_|^flag_rrt", names(ds_A), value = TRUE)
hemo_ts  <- grep("^sbp_|^dbp_|^map_|^temp_|rate_norepinephrine", names(ds_A), value = TRUE)
vent_ts  <- grep("^pao2_fio2|^peep_|^flag_vm", names(ds_A), value = TRUE)
bio_ts   <- grep("^na_|^k_|^bun_|^lactate_|^plt_|^wbc_", names(ds_A), value = TRUE)
struct_vars <- c("day_window")

model_raw_cols <- list(
  M0 = c(fixed_vars, struct_vars),
  M1 = c(fixed_vars, renal_ts, struct_vars),
  M2 = c(fixed_vars, renal_ts, hemo_ts, struct_vars),
  M3 = c(fixed_vars, renal_ts, hemo_ts, vent_ts, struct_vars),
  M4 = c(fixed_vars, renal_ts, hemo_ts, vent_ts, bio_ts, struct_vars)
)

scr_cols <- grep("^scr_", names(ds_A), value = TRUE)
model_raw_cols$S0_SCR <- c(scr_cols, "day_window")

# ============================================================
# 8) SPARSE MATRIX BUILDER (train terms fixed on AUMC)
# ============================================================

build_sparse_AUMC_train__MIMIC_test <- function(train_dt, test_dt, cols_keep) {
  cols_keep <- intersect(cols_keep, names(train_dt))
  cols_keep <- intersect(cols_keep, names(test_dt))
  stopifnot(length(cols_keep) > 0)
  
  tr <- train_dt[, ..cols_keep]; te <- test_dt[, ..cols_keep]
  setDF(tr); setDF(te)
  
  # Build design matrix with AUMC terms, then apply to MIMIC
  fml <- as.formula(paste0("~ 0 + ", paste(cols_keep, collapse=" + ")))
  
  mf_tr <- model.frame(fml, data = tr, na.action = na.pass)
  trm   <- attr(mf_tr, "terms")
  ctr   <- attr(mf_tr, "contrasts")
  
  mm_tr <- model.matrix(trm, data = mf_tr, contrasts.arg = ctr)
  mf_te <- model.frame(trm, data = te, na.action = na.pass)
  mm_te <- model.matrix(trm, data = mf_te, contrasts.arg = ctr)
  
  Xtr <- Matrix(mm_tr, sparse = TRUE)
  Xte <- Matrix(mm_te, sparse = TRUE)
  
  # Align columns (should already match, but be paranoid)
  all_cn <- union(colnames(Xtr), colnames(Xte))
  add_missing <- function(X, all_cn) {
    miss <- setdiff(all_cn, colnames(X))
    if (length(miss) == 0) return(X[, all_cn, drop=FALSE])
    Z <- Matrix(0, nrow=nrow(X), ncol=length(miss), sparse=TRUE)
    colnames(Z) <- miss
    X2 <- cbind(X, Z)
    X2[, all_cn, drop=FALSE]
  }
  Xtr <- add_missing(Xtr, all_cn)
  Xte <- add_missing(Xte, all_cn)
  stopifnot(identical(colnames(Xtr), colnames(Xte)))
  
  list(X_train = Xtr, X_test = Xte, feature_names = colnames(Xtr))
}

# ============================================================
# 9) TRAIN/EVAL - VERSION PARALLÈLE (bootstraps seulement)
# ============================================================

# Setup parallélisation
cl <- setup_parallel()  # Détecte automatiquement le nombre de cores

N_BOOT <- 500

results_mc  <- list()
results_ovr <- list()
dca_results <- list()
pred_results <- list()      # <--- AJOUTE
pred_ovr_results <- list()  # <--- (optionnel, si tu sauves OvR aussi)

total_iterations <- length(names(model_raw_cols)) * length(DAYS_TARGET)
current_iteration <- 0

for (mn in names(model_raw_cols)) {
  msg("\n==============================")
  msg("MODEL: %s", mn)
  msg("==============================")
  
  for (d in DAYS_TARGET) {
    current_iteration <- current_iteration + 1
    msg("-> Day window: %d [%d/%d]", d, current_iteration, total_iterations)
    
    tr_d <- ds_A[day_window == d]
    te_d <- ds_M[day_window == d]
    
    tr_d <- tr_d[!is.na(class)]
    te_d <- te_d[!is.na(class)]
    
    if (nrow(tr_d) < 50 || nrow(te_d) < 50) {
      msg("   Skip (too few rows): train=%d test=%d", nrow(tr_d), nrow(te_d))
      next
    }
    
    y_tr <- as.integer(factor(tr_d$class, levels = levels_y)) - 1L
    y_te <- as.integer(factor(te_d$class, levels = levels_y)) - 1L
    
    design <- build_sparse_AUMC_train__MIMIC_test(tr_d, te_d, model_raw_cols[[mn]])
    dtr <- xgb.DMatrix(design$X_train, label = y_tr, missing = NA)
    dte <- xgb.DMatrix(design$X_test,  label = y_te, missing = NA)
    
    # ---- Multiclass training (séquentiel, rapide)
    param_mc <- c(XGB_MC, list(num_class = K))
    cv_mc <- xgb.cv(
      params = param_mc,
      data = dtr,
      nrounds = NROUND_MAX,
      nfold = NFOLD,
      early_stopping_rounds = NROUND_ESTOP,
      verbose = 0
    )
    bst_mc <- xgb.train(
      params = param_mc,
      data = dtr,
      nrounds = cv_mc$best_iteration,
      verbose = 0
    )
    
    p <- predict(bst_mc, dte)
    pmat <- matrix(p, ncol = K, byrow = TRUE)
    pmat <- apply(pmat, 2, safe_prob)
    pmat <- pmat / rowSums(pmat)
    pred_int <- max.col(pmat) - 1L
    
    # ---- SAVE patient-level multiclass probabilities (for enrichment + outcome DCA)
    # Needs an ID. Try to detect it.
    id_col <- intersect(c("stay_id","icustay_id","hadm_id","subject_id"), names(te_d))
    if (length(id_col) == 0) stop("No patient ID column found in te_d (stay_id/icustay_id/hadm_id/subject_id).")
    id_col <- id_col[1]
    
    pred_pat <- data.table(
      id  = te_d[[id_col]],
      model = mn,
      day = d,
      true_class = te_d$class,
      pC1 = pmat[,1],
      pC2 = pmat[,2],
      pC3 = pmat[,3]
    )
    pred_pat[, pHigh12 := pC1 + pC2]
    
    # IMPORTANT: if duplicates (same id appears multiple times in te_d), aggregate
    # (should not happen, but safer)
    pred_pat <- pred_pat[, .(
      model = model[1],
      day = day[1],
      true_class = true_class[1],
      pC1 = mean(pC1),
      pC2 = mean(pC2),
      pC3 = mean(pC3),
      pHigh12 = mean(pHigh12)
    ), by = id]
    
    pred_results[[paste(mn,d,sep="__")]] <- pred_pat
    
    # Métriques de base
    auc_res <- multiclass_auc_ovr(y_te, pmat, K)
    brier_mc <- brier_multiclass(y_te, pmat, K)
    balacc <- balanced_acc(y_te, pred_int, K)
    cal_mc <- ece_multiclass_pooled(pmat, y_te, k = 10)
    
    save_reliability(
      cal_mc$df,
      sprintf("Frozen multiclass %s | day %d (AUMC->MIMIC)", mn, d),
      file.path(OUT_DIR, sprintf("cal_mc_frozen_%s_day%02d.png", mn, d))
    )
    
    # ---- Bootstrap CI (PARALLÈLE) ⚡
    msg("   Computing bootstrap CI (n=%d, parallel)...", N_BOOT)
    seed_base <- 20260206 + as.integer(factor(mn, levels = names(model_raw_cols))) * 1000 + d
    
    # Les fonctions bootstrap utilisent maintenant %dopar% en interne
    boot_balacc <- bootstrap_balanced_acc(y_te, pred_int, K, n_boot = N_BOOT, seed = seed_base + 1)
    boot_auc <- bootstrap_auc_macro(y_te, pmat, K, n_boot = N_BOOT, seed = seed_base + 2)
    boot_brier <- bootstrap_brier_mc(y_te, pmat, K, n_boot = N_BOOT, seed = seed_base + 3)
    boot_ece <- bootstrap_ece(y_te, pmat, K, n_boot = N_BOOT, seed = seed_base + 4)
    
    # ---- DCA
    msg("   Computing DCA...")
    dca_mc <- tryCatch({
      dca_res <- dca_multiclass_ovr(y_te, pmat, K, class_labels = levels_y)
      dca_res$model <- mn
      dca_res$day <- d
      dca_res
    }, error = function(e) {
      msg("   Warning: DCA failed - %s", e$message)
      NULL
    })
    
    if (!is.null(dca_mc)) {
      dca_results[[paste(mn, d, sep = "__")]] <- dca_mc
      p_dca <- plot_dca_multiclass(dca_mc, sprintf("DCA: %s | Day %d", mn, d))
      ggsave(file.path(OUT_DIR, sprintf("dca_mc_frozen_%s_day%02d.png", mn, d)),
             p_dca, width = 8, height = 6, dpi = 200)
    }
    
    results_mc[[paste(mn,d,sep="__")]] <- data.table(
      model = mn, day = d,
      n_train = nrow(tr_d), n_test = nrow(te_d),
      balanced_acc = balacc,
      balanced_acc_ci_lower = boot_balacc$ci_lower,
      balanced_acc_ci_upper = boot_balacc$ci_upper,
      auc_macro = auc_res$auc_macro,
      auc_macro_ci_lower = boot_auc$ci_lower,
      auc_macro_ci_upper = boot_auc$ci_upper,
      auc_weighted = auc_res$auc_weighted,
      brier_multiclass = brier_mc,
      brier_mc_ci_lower = boot_brier$ci_lower,
      brier_mc_ci_upper = boot_brier$ci_upper,
      ece_pooled = cal_mc$ece,
      ece_ci_lower = boot_ece$ci_lower,
      ece_ci_upper = boot_ece$ci_upper,
      best_nround = cv_mc$best_iteration
    )
    
    # ---- OvR binaries
    prev <- table(factor(y_tr, levels=0:(K-1)))
    w_pos <- sapply(0:(K-1), function(k) {
      npos <- as.numeric(prev[k+1])
      nneg <- length(y_tr) - npos
      if (npos == 0) return(1)
      nneg / npos
    })
    
    for (k in 0:(K-1)) {
      cls <- levels_y[k+1]
      
      ytr_bin <- as.integer(y_tr == k)
      yte_bin <- as.integer(y_te == k)
      
      dtrb <- xgb.DMatrix(design$X_train, label = ytr_bin, missing = NA)
      w <- rep(1, length(ytr_bin))
      w[ytr_bin == 1] <- w_pos[k+1]
      setinfo(dtrb, "weight", w)
      
      cv_b <- xgb.cv(
        params = XGB_BIN,
        data = dtrb,
        nrounds = NROUND_MAX,
        nfold = NFOLD,
        early_stopping_rounds = NROUND_ESTOP,
        verbose = 0
      )
      bst_b <- xgb.train(
        params = XGB_BIN,
        data = dtrb,
        nrounds = cv_b$best_iteration,
        verbose = 0
      )
      
      prob <- safe_prob(predict(bst_b, dte))
      pred <- as.integer(prob > 0.5)
      
      roc_o <- tryCatch(pROC::roc(yte_bin, prob, quiet = TRUE), error = function(e) NULL)
      auc <- if (is.null(roc_o)) NA_real_ else as.numeric(pROC::auc(roc_o))
      acc <- mean(pred == yte_bin)
      sens <- if (sum(yte_bin==1) > 0) sum(pred==1 & yte_bin==1) / sum(yte_bin==1) else NA_real_
      spec <- if (sum(yte_bin==0) > 0) sum(pred==0 & yte_bin==0) / sum(yte_bin==0) else NA_real_
      brier <- mean((prob - yte_bin)^2)
      
      # Bootstrap OvR (PARALLÈLE) ⚡
      boot_auc_ovr <- bootstrap_auc_binary(yte_bin, prob, n_boot = N_BOOT, seed = seed_base + 10 + k)
      
      cal_b <- ece_binary(prob, yte_bin, k = 10)
      save_reliability(
        cal_b$df,
        sprintf("Frozen OvR %s %s | day %d", mn, cls, d),
        file.path(OUT_DIR, sprintf("cal_ovr_frozen_%s_%s_day%02d.png", mn, cls, d))
      )
      
      dca_bin <- tryCatch({
        dca_binary(yte_bin, prob, label = sprintf("%s_%s", mn, cls))
      }, error = function(e) NULL)
      
      if (!is.null(dca_bin)) {
        p_dca_bin <- plot_dca_binary(dca_bin$df, sprintf("DCA: %s %s | Day %d", mn, cls, d))
        ggsave(file.path(OUT_DIR, sprintf("dca_ovr_frozen_%s_%s_day%02d.png", mn, cls, d)),
               p_dca_bin, width = 6, height = 4, dpi = 200)
      }
      
      results_ovr[[paste(mn,d,cls,sep="__")]] <- data.table(
        model = mn, day = d, class = cls,
        n_train = nrow(tr_d), n_test = nrow(te_d),
        auc = auc,
        auc_ci_lower = boot_auc_ovr$ci_lower,
        auc_ci_upper = boot_auc_ovr$ci_upper,
        accuracy = acc, 
        sensitivity = sens, 
        specificity = spec,
        brier = brier, 
        ece = cal_b$ece,
        best_nround = cv_b$best_iteration
      )
    }
    
    rm(design, dtr, dte, bst_mc, cv_mc); gc()
  }
}

# Arrêter la parallélisation
stop_parallel(cl)

# Reste du code identique...
mc_dt  <- rbindlist(results_mc,  fill = TRUE)
ovr_dt <- rbindlist(results_ovr, fill = TRUE)
dca_dt <- rbindlist(dca_results, fill = TRUE)
pred_dt <- rbindlist(pred_results, fill = TRUE)


fwrite(mc_dt,  file.path(OUT_DIR, "FROZEN_multiclass_AUMC_train__MIMIC_test_byDay_withCI_bis.csv"))
fwrite(ovr_dt, file.path(OUT_DIR, "FROZEN_OvR_AUMC_train__MIMIC_test_byDay_withCI_bis.csv"))
fwrite(dca_dt, file.path(OUT_DIR, "DCA_results_frozen_bis.csv"))
saveRDS(pred_dt, file.path(OUT_DIR, "pred_dt.rds"))
fwrite(pred_dt, file.path(OUT_DIR, "pred_dt.csv"))

# Simple curves
p1 <- ggplot(mc_dt, aes(x = day, y = balanced_acc, color = model, group = model)) +
  geom_line() + geom_point() + theme_bw() +
  scale_x_continuous(breaks = DAYS_TARGET) +
  ylim(0,1) +
  labs(title = "Frozen external validation (AUMC->MIMIC): balanced accuracy", x = "Day", y = "Balanced accuracy")
ggsave(file.path(OUT_DIR, "curve_frozen_balanced_acc_bis.png"), p1, width = 8, height = 4, dpi = 200)

p2 <- ggplot(mc_dt, aes(x = day, y = brier_multiclass, color = model, group = model)) +
  geom_line() + geom_point() + theme_bw() +
  scale_x_continuous(breaks = DAYS_TARGET) +
  labs(title = "Frozen external validation (AUMC->MIMIC): multiclass Brier", x = "Day", y = "Brier")
ggsave(file.path(OUT_DIR, "curve_frozen_brier_bis.png"), p2, width = 8, height = 4, dpi = 200)


# ============================================================
# 10) OPTIONAL: RECALIBRATION / SITE-ADAPTATION ON MIMIC (80/20)
# ============================================================

if (DO_RECAL_MIMIC) {
  msg("\n== RECALIBRATION ON MIMIC (site adaptation), 80/20 by stay_id, per day ==")
  
  # patient-level split on MIMIC stays (same split used for all days to be clean)
  mim_ids <- unique(ds_M$stay_id)
  mim_lbl <- unique(ds_M[, .(stay_id, class)])
  setkey(mim_lbl, stay_id)
  
  # stratified split by class
  split_ids_strat <- function(stay_id, y, train_prop = 0.8, seed = 1) {
    set.seed(seed)
    dt <- data.table(stay_id = stay_id, y = y)
    dt <- unique(dt, by = "stay_id")
    tr_ids <- dt[, {
      n <- .N
      n_tr <- max(1L, floor(n * train_prop))
      .(stay_id = sample(stay_id, n_tr))
    }, by = y]$stay_id
    list(train = tr_ids, test = setdiff(dt$stay_id, tr_ids))
  }
  
  spM <- split_ids_strat(mim_lbl$stay_id, mim_lbl$class, train_prop = RECAL_PROP, seed = 20260206)
  mim_train_ids <- spM$train
  mim_test_ids  <- spM$test
  
  recal_mc <- list()
  recal_ovr <- list()
  
  for (mn in names(model_raw_cols)) {
    for (d in DAYS_TARGET) {
      tr_d <- ds_M[stay_id %in% mim_train_ids & day_window == d]
      te_d <- ds_M[stay_id %in% mim_test_ids  & day_window == d]
      if (nrow(tr_d) < 50 || nrow(te_d) < 50) next
      
      y_tr <- as.integer(factor(tr_d$class, levels = levels_y)) - 1L
      y_te <- as.integer(factor(te_d$class, levels = levels_y)) - 1L
      
      design <- build_sparse_AUMC_train__MIMIC_test(tr_d, te_d, model_raw_cols[[mn]])
      dtr <- xgb.DMatrix(design$X_train, label = y_tr, missing = NA)
      dte <- xgb.DMatrix(design$X_test,  label = y_te, missing = NA)
      
      param_mc <- c(XGB_MC, list(num_class = K))
      cv_mc <- xgb.cv(param_mc, dtr, NROUND_MAX, NFOLD, early_stopping_rounds = NROUND_ESTOP, verbose = 0)
      bst_mc <- xgb.train(param_mc, dtr, nrounds = cv_mc$best_iteration, verbose = 0)
      
      p <- predict(bst_mc, dte)
      pmat <- matrix(p, ncol = K, byrow = TRUE)
      pmat <- apply(pmat, 2, safe_prob); pmat <- pmat / rowSums(pmat)
      pred_int <- max.col(pmat) - 1L
      
      auc_res <- multiclass_auc_ovr(y_te, pmat, K)
      brier_mc <- brier_multiclass(y_te, pmat, K)
      balacc <- balanced_acc(y_te, pred_int, K)
      cal_mc <- ece_multiclass_pooled(pmat, y_te, k = 10)
      
      recal_mc[[paste(mn,d,sep="__")]] <- data.table(
        model = mn, day = d,
        n_train = nrow(tr_d), n_test = nrow(te_d),
        balanced_acc = balacc,
        auc_macro = auc_res$auc_macro,
        auc_weighted = auc_res$auc_weighted,
        brier_multiclass = brier_mc,
        ece_pooled = cal_mc$ece,
        best_nround = cv_mc$best_iteration
      )
      
      # OvR
      prev <- table(factor(y_tr, levels=0:(K-1)))
      w_pos <- sapply(0:(K-1), function(k) {
        npos <- as.numeric(prev[k+1]); nneg <- length(y_tr) - npos
        if (npos == 0) return(1); nneg / npos
      })
      
      for (k in 0:(K-1)) {
        cls <- levels_y[k+1]
        ytr_bin <- as.integer(y_tr == k)
        yte_bin <- as.integer(y_te == k)
        
        dtrb <- xgb.DMatrix(design$X_train, label = ytr_bin, missing = NA)
        w <- rep(1, length(ytr_bin)); w[ytr_bin==1] <- w_pos[k+1]
        setinfo(dtrb, "weight", w)
        
        cv_b <- xgb.cv(XGB_BIN, dtrb, NROUND_MAX, NFOLD, early_stopping_rounds = NROUND_ESTOP, verbose = 0)
        bst_b <- xgb.train(XGB_BIN, dtrb, nrounds = cv_b$best_iteration, verbose = 0)
        
        prob <- safe_prob(predict(bst_b, dte))
        pred <- as.integer(prob > 0.5)
        
        roc_o <- tryCatch(pROC::roc(yte_bin, prob, quiet = TRUE), error = function(e) NULL)
        auc <- if (is.null(roc_o)) NA_real_ else as.numeric(pROC::auc(roc_o))
        acc <- mean(pred == yte_bin)
        sens <- if (sum(yte_bin==1) > 0) sum(pred==1 & yte_bin==1) / sum(yte_bin==1) else NA_real_
        spec <- if (sum(yte_bin==0) > 0) sum(pred==0 & yte_bin==0) / sum(yte_bin==0) else NA_real_
        brier <- mean((prob - yte_bin)^2)
        cal_b <- ece_binary(prob, yte_bin, k = 10)
        
        recal_ovr[[paste(mn,d,cls,sep="__")]] <- data.table(
          model = mn, day = d, class = cls,
          n_train = nrow(tr_d), n_test = nrow(te_d),
          auc = auc, accuracy = acc, sensitivity = sens, specificity = spec,
          brier = brier, ece = cal_b$ece,
          best_nround = cv_b$best_iteration
        )
      }
      rm(design, dtr, dte, bst_mc, cv_mc); gc()
    }
  }
  
  recal_mc_dt <- rbindlist(recal_mc, fill=TRUE)
  recal_ovr_dt <- rbindlist(recal_ovr, fill=TRUE)
  fwrite(recal_mc_dt,  file.path(OUT_DIR, "RECAL_MIMIC_multiclass_80_20_byDay.csv"))
  fwrite(recal_ovr_dt, file.path(OUT_DIR, "RECAL_MIMIC_OvR_80_20_byDay.csv"))
  
  # Compare frozen vs recal (multiclass)
  mc_dt$setting <- "Frozen AUMC->MIMIC"
  recal_mc_dt$setting <- "Recal MIMIC 80/20"
  cmp <- rbind(mc_dt, recal_mc_dt, fill=TRUE)
  
  p_cmp <- ggplot(cmp, aes(day, balanced_acc, color = setting, group = setting)) +
    geom_line() + geom_point() + theme_bw() +
    scale_x_continuous(breaks = DAYS_TARGET) + ylim(0,1) +
    labs(title="Balanced accuracy over time", x="Day", y="Balanced accuracy")
  ggsave(file.path(OUT_DIR, "compare_balanced_acc_frozen_vs_recal.png"), p_cmp, width=7, height=4, dpi=200)
}

msg("\n✅ DONE. Outputs in: %s", OUT_DIR)

#### Ajout
# ============================================================
# GRAPHIQUES AVEC INTERVALLES DE CONFIANCE
# ============================================================

# Balanced accuracy avec IC
p_balacc_ci <- ggplot(mc_dt, aes(x = day, y = balanced_acc, color = model, group = model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = balanced_acc_ci_lower, ymax = balanced_acc_ci_upper, fill = model),
              alpha = 0.2, color = NA) +
  theme_bw(base_size = 12) +
  scale_x_continuous(breaks = DAYS_TARGET) +
  ylim(0, 1) +
  labs(
    title = "Frozen external validation (AUMC→MIMIC): Balanced Accuracy with 95% CI",
    x = "Day",
    y = "Balanced Accuracy",
    caption = sprintf("Bootstrap CI (n=%d)", N_BOOT)
  ) +
  theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "curve_frozen_balanced_acc_withCI.png"), 
       p_balacc_ci, width = 10, height = 6, dpi = 300)

# AUC macro avec IC
p_auc_ci <- ggplot(mc_dt, aes(x = day, y = auc_macro, color = model, group = model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = auc_macro_ci_lower, ymax = auc_macro_ci_upper, fill = model),
              alpha = 0.2, color = NA) +
  theme_bw(base_size = 12) +
  scale_x_continuous(breaks = DAYS_TARGET) +
  ylim(0, 1) +
  labs(
    title = "Frozen external validation: AUC Macro with 95% CI",
    x = "Day",
    y = "AUC (Macro)",
    caption = sprintf("Bootstrap CI (n=%d)", N_BOOT)
  ) +
  theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "curve_frozen_auc_macro_withCI.png"), 
       p_auc_ci, width = 10, height = 6, dpi = 300)

# Brier score avec IC
p_brier_ci <- ggplot(mc_dt, aes(x = day, y = brier_multiclass, color = model, group = model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = brier_mc_ci_lower, ymax = brier_mc_ci_upper, fill = model),
              alpha = 0.2, color = NA) +
  theme_bw(base_size = 12) +
  scale_x_continuous(breaks = DAYS_TARGET) +
  labs(
    title = "Frozen external validation: Brier Score with 95% CI",
    x = "Day",
    y = "Brier Score (lower is better)",
    caption = sprintf("Bootstrap CI (n=%d)", N_BOOT)
  ) +
  theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "curve_frozen_brier_withCI.png"), 
       p_brier_ci, width = 10, height = 6, dpi = 300)

# DCA comparaison entre modèles (pour un jour donné, ex: day 3)
if (nrow(dca_dt) > 0) {
  dca_day3 <- dca_dt %>% filter(day == 3, strategy != "Treat All", strategy != "Treat None")
  
  if (nrow(dca_day3) > 0) {
    p_dca_compare <- ggplot(dca_day3, aes(x = threshold, y = net_benefit, color = model)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      facet_wrap(~ class, ncol = 2) +
      theme_bw(base_size = 11) +
      labs(
        title = "Decision Curve Analysis - Day 3 (Model Comparison)",
        x = "Threshold Probability",
        y = "Net Benefit",
        color = "Model"
      ) +
      theme(legend.position = "bottom") +
      xlim(0, 1)
    
    ggsave(file.path(OUT_DIR, "dca_comparison_day3.png"), 
           p_dca_compare, width = 10, height = 8, dpi = 300)
  } else {
    msg("Warning: No DCA data for day 3")
  }
} else {
  msg("Warning: No DCA results to plot")
}

ggsave(file.path(OUT_DIR, "dca_comparison_day3.png"), 
       p_dca_compare, width = 10, height = 8, dpi = 300)

# ============================================================
# TABLEAU FORMATÉ POUR PUBLICATION
# ============================================================

format_ci <- function(mean, lower, upper, digits = 3) {
  sprintf("%.*f (%.*f-%.*f)", digits, mean, digits, lower, digits, upper)
}

# Table pour Day 3 (exemple)
mc_table <- mc_dt %>%
  filter(day == 2) %>%
  mutate(
    `Balanced Accuracy` = format_ci(balanced_acc, balanced_acc_ci_lower, balanced_acc_ci_upper),
    `AUC Macro` = format_ci(auc_macro, auc_macro_ci_lower, auc_macro_ci_upper),
    `Brier Score` = format_ci(brier_multiclass, brier_mc_ci_lower, brier_mc_ci_upper),
    `ECE` = format_ci(ece_pooled, ece_ci_lower, ece_ci_upper)
  ) %>%
  select(Model = model, `Balanced Accuracy`, `AUC Macro`, `Brier Score`, `ECE`, 
         `N Train` = n_train, `N Test` = n_test)

fwrite(mc_table, file.path(OUT_DIR, "table1_performance_day2_withCI.csv"))

# Version LaTeX
library(xtable)
xt <- xtable(mc_table, 
             caption = "External validation performance at Day 3 (AUMC→MIMIC) with 95\\% bootstrap confidence intervals",
             label = "tab:performance")
print(xt, file = file.path(OUT_DIR, "table1_latex.tex"), include.rownames = FALSE)




library(data.table)

# results_mc: list of data.table (1 row each)
tab_mc <- rbindlist(results_mc, use.names = TRUE, fill = TRUE)
setorder(tab_mc, model, day)

# (optionnel) ordre des modèles
tab_mc[, model := factor(model, levels = c("M0","S0_SCR","M1","M2","M3","M4"))]
setorder(tab_mc, model, day)

# tables “prêtes à plotter”
tab_auc_macro    <- tab_mc[, .(model, day, auc_macro)]
tab_auc_weighted <- tab_mc[, .(model, day, auc_weighted)]
tab_balacc       <- tab_mc[, .(model, day, balanced_acc)]
tab_brier        <- tab_mc[, .(model, day, brier_multiclass)]
tab_ece          <- tab_mc[, .(model, day, ece_pooled)]

# si tu veux exporter
fwrite(tab_mc, file.path(OUT_DIR, "tab_mc_AUMC_to_MIMIC.xsls"))
fwrite(tab_auc_macro, file.path(OUT_DIR, "tab_auc_macro.xsls"))
fwrite(tab_auc_weighted, file.path(OUT_DIR, "tab_auc_weighted.xsls"))

library(data.table)

# results_ovr: list of data.table (1 row each)
tab_ovr <- rbindlist(results_ovr, use.names = TRUE, fill = TRUE)
setorder(tab_ovr, model, class, day)

# (optionnel) ordre modèles + classes
tab_ovr[, model := factor(model, levels = c("M0","S0_SCR","M1","M2","M3","M4"))]
tab_ovr[, class := factor(class, levels = c("C1","C2","C3","C4"))]
setorder(tab_ovr, model, class, day)

# tables “prêtes à plotter”
tab_auc_by_class <- tab_ovr[, .(model, day, class, auc)]
tab_brier_by_class <- tab_ovr[, .(model, day, class, brier)]
tab_ece_by_class <- tab_ovr[, .(model, day, class, ece)]
tab_sensspec_by_class <- tab_ovr[, .(model, day, class, sensitivity, specificity)]

# export si besoin
fwrite(tab_ovr, file.path(OUT_DIR, "tab_ovr_AUMC_to_MIMIC.csv"))
fwrite(tab_auc_by_class, file.path(OUT_DIR, "tab_auc_by_class.csv"))

tab_auc_long <- melt(
  tab_mc,
  id.vars = c("model","day"),
  measure.vars = c("auc_macro","auc_weighted"),
  variable.name = "metric",
  value.name = "auc"
)
# metric = "auc_macro" ou "auc_weighted"

# ============================================================
# CALIBRATION PLOTS
# ============================================================
library(magick)

OUT_DIR <- file.path(Sys.getenv("AKI_BASE_DIR"), "XGB_external_AUMC_train__MIMIC_G3")

DAYS   <- 0:7
MODELS <- c("M0","M1","M2","M3","M4","S0_SCR")
CLASS  <- "C3"   # changer ici : "C1", "C2", "C3", "C4"

img_path <- function(model, class, day) {
  file.path(OUT_DIR, sprintf("cal_ovr_frozen_%s_%s_day%02d.png", model, class, day))
}

read_or_placeholder <- function(path, label) {
  if (file.exists(path)) {
    img <- image_read(path)
  } else {
    img <- image_blank(width = 700, height = 700, color = "white")
    img <- image_annotate(img, paste0("MISSING\n", label),
                          size = 40, color = "red", gravity = "center")
  }
  image_resize(img, "700x700!")
}

# Construire la matrice d’images
img_mat <- lapply(DAYS, function(d) {
  lapply(MODELS, function(m) {
    p <- img_path(m, CLASS, d)
    read_or_placeholder(p, paste(m, CLASS, d))
  })
})

# Ajouter titres colonnes (modèles)
annotate_top <- function(img, text) {
  header <- image_blank(width = image_info(img)$width,
                        height = 80, color = "white") |>
    image_annotate(text, size = 40, color = "black", gravity = "center")
  image_append(c(header, img), stack = TRUE)
}

img_mat[[1]] <- Map(function(img, m) annotate_top(img, m),
                    img_mat[[1]], MODELS)

# Ajouter titres lignes (jours)
annotate_left <- function(row_imgs, text) {
  strip <- image_blank(width = 120,
                       height = image_info(row_imgs[[1]])$height,
                       color = "white") |>
    image_annotate(text, size = 45, color = "black",
                   gravity = "center", degrees = -90)
  c(list(strip), row_imgs)
}

img_mat <- Map(function(row_imgs, d)
  annotate_left(row_imgs, paste0("Day ", d)),
  img_mat, DAYS)

# Assembler
row_imgs <- lapply(img_mat, function(row_list)
  image_append(do.call(c, row_list), stack = FALSE))

grid_img <- image_append(do.call(c, row_imgs), stack = TRUE)

# Sauvegarde
out_png <- file.path(OUT_DIR, sprintf("CAL_OVR_%s_grid.png", CLASS))
out_pdf <- file.path(OUT_DIR, sprintf("CAL_OVR_%s_grid.pdf", CLASS))

image_write(grid_img, path = out_png, format = "png")
image_write(grid_img, path = out_pdf, format = "pdf")

cat("Saved:\n", out_png, "\n", out_pdf, "\n")




##### Sauvegarde finale
saveRDS(mc_dt,  file.path(OUT_DIR, "mc_dt.rds"))
saveRDS(ovr_dt, file.path(OUT_DIR, "ovr_dt.rds"))
saveRDS(dca_dt, file.path(OUT_DIR, "dca_dt.rds"))

writeLines(capture.output(sessionInfo()),
           con = file.path(OUT_DIR, "sessionInfo.txt"))

saveRDS(list(
  seed = 20260130,
  DAYS_TARGET = DAYS_TARGET,
  model_raw_cols = model_raw_cols,
  XGB_MC = XGB_MC,
  XGB_BIN = XGB_BIN,
  NROUND_MAX = NROUND_MAX,
  NROUND_ESTOP = NROUND_ESTOP,
  NFOLD = NFOLD,
  N_BOOT = N_BOOT,
  DO_RECAL_MIMIC = DO_RECAL_MIMIC,
  RECAL_PROP = RECAL_PROP
), file.path(OUT_DIR, "run_config.rds"))
