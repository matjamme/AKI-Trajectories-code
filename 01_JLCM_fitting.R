#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  # Uses the active R library paths; no user-specific HPC library path.
  library(lcmm)
  library(survival)
  library(data.table)
})

cat("==== Jointlcmm (quadratic) + winsor + competing risks + horizon variable ====\n")

# ----------------------------
# Args
# ----------------------------
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  res <- list(
    db_path    = NULL,
    dataset    = "DATASET",   # label used in filenames/logs: AUMC / MIMIC / etc
    out_prefix = "CR2_quad",
    out_dir    = NULL,
    
    n_cores    = 16L,
    reps       = 30L,
    maxiter    = 200L,
    g_min      = 2L,
    g_max      = 5L,
    
    horizon_d  = 14L,         # 7 or 14
    time_var   = "time_s",    # "time_s" (0-1) or "time_day" (days)
    with_rrt   = 0L,          # 0 = Surv~1, 1 = Surv~rrt_td
    
    seed       = 2025L,
    winsor_p   = 0.99,
    maxiter_g1 = 500L
  )
  
  if (length(args) > 0L) {
    for (a in args) {
      if (!grepl("^--", a)) next
      split <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
      key <- split[1]
      val <- if (length(split) > 1) split[2] else NA_character_
      
      if (key == "db_path")    res$db_path <- val
      if (key == "dataset")    res$dataset <- val
      if (key == "out_prefix") res$out_prefix <- val
      if (key == "out_dir")    res$out_dir <- val
      
      if (key == "n_cores")    res$n_cores <- as.integer(val)
      if (key == "reps")       res$reps <- as.integer(val)
      if (key == "maxiter")    res$maxiter <- as.integer(val)
      if (key == "g_min")      res$g_min <- as.integer(val)
      if (key == "g_max")      res$g_max <- as.integer(val)
      
      if (key == "horizon_d")  res$horizon_d <- as.integer(val)
      if (key == "time_var")   res$time_var <- as.character(val)
      if (key == "with_rrt")   res$with_rrt <- as.integer(val)
      
      if (key == "seed")       res$seed <- as.integer(val)
      if (key == "winsor_p")   res$winsor_p <- as.numeric(val)
      if (key == "maxiter_g1") res$maxiter_g1 <- as.integer(val)
    }
  }
  
  res
}

cfg <- parse_args()
if (is.null(cfg$db_path) || cfg$db_path == "") stop("Argument --db_path manquant")
if (is.null(cfg$out_dir) || cfg$out_dir == "") cfg$out_dir <- getwd()
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

cat("[PKG] lcmm version:", as.character(utils::packageVersion("lcmm")), "\n")
cat("[CFG] dataset     :", cfg$dataset, "\n")
cat("[CFG] db_path     :", cfg$db_path, "\n")
cat("[CFG] out_dir     :", cfg$out_dir, "\n")
cat("[CFG] out_prefix  :", cfg$out_prefix, "\n")
cat("[CFG] n_cores     :", cfg$n_cores, "\n")
cat("[CFG] reps        :", cfg$reps, "\n")
cat("[CFG] maxiter     :", cfg$maxiter, "\n")
cat("[CFG] g_min..g_max:", cfg$g_min, "->", cfg$g_max, "\n")
cat("[CFG] horizon_d   :", cfg$horizon_d, "days\n")
cat("[CFG] time_var    :", cfg$time_var, "\n")
cat("[CFG] with_rrt    :", cfg$with_rrt, "(0=Surv~1, 1=Surv~rrt_td)\n")
cat("[CFG] winsor_p    :", cfg$winsor_p, "\n")
cat("[CFG] seed        :", cfg$seed, "\n")

# ----------------------------
# Threads + seed
# ----------------------------
Sys.setenv(
  OMP_NUM_THREADS      = cfg$n_cores,
  MKL_NUM_THREADS      = cfg$n_cores,
  OPENBLAS_NUM_THREADS = cfg$n_cores
)
set.seed(cfg$seed)

# ----------------------------
# Load
# ----------------------------
cat("[LOAD] Reading:", cfg$db_path, "\n")
if (!file.exists(cfg$db_path)) stop("Fichier RDS introuvable: ", cfg$db_path)

dt <- as.data.table(readRDS(cfg$db_path))

# ----------------------------
# Validate columns
# ----------------------------
req <- c("id", "scr_mgdl", "tte_pos", "event2", "time_day", "time_s")
if (cfg$with_rrt == 1L) req <- c(req, "rrt_td")

miss <- setdiff(req, names(dt))
if (length(miss) > 0L) stop("Variables manquantes : ", paste(miss, collapse = ", "))

# enforce types
dt[, id := as.integer(id)]
dt[, tte_pos := as.numeric(tte_pos)]
dt[, event2 := as.integer(event2)]
dt[, time_day := as.numeric(time_day)]
dt[, time_s := as.numeric(time_s)]
dt[, scr_mgdl := as.numeric(scr_mgdl)]
if (cfg$with_rrt == 1L) dt[, rrt_td := as.numeric(rrt_td)]

cat("[DATA] rows     :", nrow(dt), "\n")
cat("[DATA] patients :", uniqueN(dt$id), "\n")

pt0 <- unique(dt[, .(id, event2)])
cat("[DATA] event2 distribution (patient-level raw):\n")
print(prop.table(table(pt0$event2, useNA="ifany")))

# ----------------------------
# Winsorize SCr
# ----------------------------
q <- quantile(dt$scr_mgdl, probs = cfg$winsor_p, na.rm = TRUE)
cat("[WINSOR] p=", cfg$winsor_p, " quantile(sScr)=", as.numeric(q), "\n", sep="")
dt[, scr_w99 := pmin(scr_mgdl, as.numeric(q))]

# ----------------------------
# Apply horizon consistently
# - Longitudinal: keep time_day <= horizon_d
# - Endpoints: censor at horizon_d
# ----------------------------
dt <- dt[is.finite(scr_w99) & is.finite(time_day) & is.finite(time_s) & is.finite(tte_pos) & !is.na(event2)]
dt <- dt[time_day >= 0 & time_day <= cfg$horizon_d]

# IMPORTANT: preserve raw tte for event recoding
dt[, tte_raw := tte_pos]
dt[, tte_pos := pmin(tte_raw, cfg$horizon_d)]
dt[, event2  := ifelse(tte_raw <= cfg$horizon_d, event2, 0L)]
dt[, tte_raw := NULL]

# choose time covariate
if (!cfg$time_var %in% c("time_s", "time_day")) stop("time_var must be time_s or time_day")
tvar <- cfg$time_var

# minimal modeling columns
keep_cols <- c("id", tvar, "scr_w99", "tte_pos", "event2")
if (cfg$with_rrt == 1L) keep_cols <- c(keep_cols, "rrt_td")
dt <- dt[, ..keep_cols]

setorder(dt, id, get(tvar))

cat("[DATA] After horizon filter: rows=", nrow(dt), " patients=", uniqueN(dt$id), "\n", sep="")
pt <- unique(dt[, .(id, event2)])
cat("[DATA] event2 distribution (patient-level after censor):\n")
print(prop.table(table(pt$event2, useNA="ifany")))

dat <- as.data.frame(dt)

# ----------------------------
# Survival formula
# ----------------------------
surv_form <- if (cfg$with_rrt == 1L) {
  Surv(tte_pos, event2) ~ rrt_td
} else {
  Surv(tte_pos, event2) ~ 1
}

# ----------------------------
# Fit G=1 seed
# ----------------------------
cat("\n[SEED] Fitting G=1 quadratic | time=", tvar,
    " | H=", cfg$horizon_d, " | with_rrt=", cfg$with_rrt, "...\n", sep="")

fixed_form  <- as.formula(sprintf("scr_w99 ~ %s + I(%s^2)", tvar, tvar))
random_form <- as.formula(sprintf("~ %s", tvar))

G1 <- Jointlcmm(
  fixed      = fixed_form,
  random     = random_form,
  subject    = "id",
  idiag      = FALSE,
  survival   = surv_form,
  hazard     = c("Weibull","Weibull"),
  hazardtype = "Specific",
  link       = "linear",
  logscale   = TRUE,
  ng         = 1,
  data       = dat,
  na.action  = 1,
  maxiter    = cfg$maxiter_g1,
  verbose    = TRUE
)

cat("[G1] conv=", G1$conv, " logLik=", G1$loglik, " BIC=", G1$BIC, " AIC=", G1$AIC, "\n")

# Only use a converged G1 model to initialise the grid search.
if (G1$conv != 1L) {
  stop("G1 did not converge (conv=", G1$conv,
       "). Check the model fit before running gridsearch.")
}

g1_file <- file.path(cfg$out_dir,
                     sprintf("%s_joint_G1_%s_quad_w99_cr_H%dd_%s.rds",
                             cfg$dataset, cfg$out_prefix, cfg$horizon_d, tvar))
saveRDS(G1, g1_file)
cat("[SAVE] ", g1_file, "\n", sep="")

# ----------------------------
# Gridsearch G = g_min..g_max
# ----------------------------
res_list <- list()
res_list[["G1"]] <- data.frame(G=1L, logLik=G1$loglik, BIC=G1$BIC, AIC=G1$AIC, conv=G1$conv)

for (G in cfg$g_min:cfg$g_max) {
  cat("\n[GRID] G=", G, " | rep=", cfg$reps, " | maxiter=", cfg$maxiter, "\n", sep="")
  if (G == 1L) next
  
  jm_G <- tryCatch(
    {
      gridsearch(
        rep     = cfg$reps,
        maxiter = cfg$maxiter,
        minit   = G1,
        m = Jointlcmm(
          fixed      = fixed_form,
          random     = random_form,
          mixture    = as.formula(sprintf("~ %s + I(%s^2)", tvar, tvar)),
          subject    = "id",
          idiag      = FALSE,
          survival   = surv_form,
          hazard     = c("Weibull","Weibull"),
          hazardtype = "Specific",
          link       = "linear",
          logscale   = TRUE,
          ng         = G,
          data       = dat,
          na.action  = 1
        )
      )
    },
    error = function(e) {
      cat("[GRID] ERROR G=", G, " : ", e$message, "\n", sep="")
      NULL
    }
  )
  if (is.null(jm_G)) next
  
  cat("[G", G, "] conv=", jm_G$conv, " logLik=", jm_G$loglik, " BIC=", jm_G$BIC, " AIC=", jm_G$AIC, "\n", sep="")
  
  out_rds <- file.path(cfg$out_dir,
                       sprintf("%s_joint_G%d_%s_quad_w99_cr_H%dd_%s.rds",
                               cfg$dataset, G, cfg$out_prefix, cfg$horizon_d, tvar))
  saveRDS(jm_G, out_rds)
  cat("[SAVE] ", out_rds, "\n", sep="")
  
  res_list[[paste0("G", G)]] <- data.frame(G=G, logLik=jm_G$loglik, BIC=jm_G$BIC, AIC=jm_G$AIC, conv=jm_G$conv)
}

sel_tab <- do.call(rbind, res_list)
sel_tab <- sel_tab[order(sel_tab$G), ]

cat("\n[MODEL SELECTION]\n")
print(sel_tab)

out_csv <- file.path(cfg$out_dir,
                     sprintf("%s_joint_selection_%s_quad_w99_cr_H%dd_%s.csv",
                             cfg$dataset, cfg$out_prefix, cfg$horizon_d, tvar))
fwrite(sel_tab, out_csv)
cat("[SAVE] ", out_csv, "\n", sep="")

cat("\n==== END gridsearch ====\n")
