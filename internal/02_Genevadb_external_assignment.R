## ============================================================
## EXTERNAL VALIDATION — AUMC trained Jointlcmm -> GENEVE data
## Robust:
##  - scr_w99 harmonised using q99(AUMC)
##  - robust TIMEVAR detection (predictY) on GENEVE
##  - ensure TIMEVAR column exists in BOTH AUMC & GENEVE
## Outputs (per model):
##  - ppGENEVE + class sizes
##  - heatmap mean posterior probs
##  - traj GENEVE observed (median) vs predicted (model)
##  - CIF GENEVE (empirical) by predicted class
##  - CIF compare AUMC (solid) vs GENEVE (dashed), by class & event
## ============================================================
rm(list = ls())
graphics.off()

suppressPackageStartupMessages({
  library(lcmm)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(cmprsk)
  library(tibble)
})

## -------------------------
## CONFIG
## -------------------------
DAY_MAX <- 14
Y_LIM   <- c(0, 5)

## AUMC data (reference dataset for q99 + CIF compare)
DB_AUMC_RDS <- Sys.getenv("AKI_AUMC_RDS")

## GENEVE data (external dataset)
DB_GENEVE_RDS <- Sys.getenv("AKI_EXTERNAL_RDS")

## Directory containing TRAINED models (AUMC fits)
MODELS_DIR <- Sys.getenv("AKI_MODELS_DIR")

## Optional: load only models matching pattern
MODEL_PATTERN <- if (nzchar(Sys.getenv("AKI_MODEL_PATTERN"))) Sys.getenv("AKI_MODEL_PATTERN") else NULL

## Output
OUT_DIR <- Sys.getenv("AKI_OUT_DIR")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "tables"),  showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "figures/traj"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "figures/cif"),  showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "figures/cif_compare"),  showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "figures/heatmap"), showWarnings = FALSE, recursive = TRUE)

## -------------------------
## 0) Load datasets
## -------------------------
dtA <- as.data.table(readRDS(DB_AUMC_RDS))
dtG <- as.data.table(readRDS(DB_GENEVE_RDS))

# garder une copie de l'id original
dtA[, id_chr := as.character(id)]
dtG[, id_chr := as.character(id)]

# construire un id numérique pour lcmm/predictClass
all_ids <- unique(c(dtA$id_chr, dtG$id_chr))
id_map <- data.table(
  id_chr = all_ids,
  id_num = seq_along(all_ids)
)

setkey(id_map, id_chr)
setkey(dtA, id_chr)
setkey(dtG, id_chr)

dtA <- id_map[dtA]
dtG <- id_map[dtG]

# lcmm veut un subject numérique
dtA[, id := as.integer(id_num)]
dtG[, id := as.integer(id_num)]

stopifnot("scr_mgdl" %in% names(dtA), "scr_mgdl" %in% names(dtG))
stopifnot(all(c("id","tte_pos","event2") %in% names(dtA)))
stopifnot(all(c("id","tte_pos","event2") %in% names(dtG)))

## -------------------------
## 1) Harmonise scr_w99 using q99(AUMC)
## -------------------------
q99_A <- quantile(dtA$scr_mgdl, probs = 0.99, na.rm = TRUE)
cat("[SCr] q99(AUMC scr_mgdl) =", as.numeric(q99_A), "mg/dL\n")

dtA[, scr_w99 := pmin(scr_mgdl, as.numeric(q99_A))]
dtG[, scr_w99 := pmin(scr_mgdl, as.numeric(q99_A))]

dtG[, event2 := fifelse(tte_pos > 14, 0L, event2)]
dtG[, tte_pos := pmin(tte_pos, 14)]

dtG[, .N, by = event2][order(event2)]
summary(dtG$tte_pos)

dtG[, .N, by = event2][order(event2)]

dtG[, .(
  min_tte = min(tte_pos, na.rm = TRUE),
  med_tte = median(tte_pos, na.rm = TRUE),
  max_tte = max(tte_pos, na.rm = TRUE)
), by = event2][order(event2)]


## -------------------------
## 2) Harmonise / ensure time variables exist
## -------------------------
normalize_time_s_to_01 <- function(dt, day_max = 14, label = "dt") {
  if (!("time_s" %in% names(dt))) return(dt)
  x <- dt$time_s
  x <- x[is.finite(x)]
  if (!length(x)) return(dt)
  
  if (max(x, na.rm = TRUE) > 1.5) {
    den <- day_max * 24 * 3600
    cat("[TIME]", label, "time_s not in [0,1] -> assuming seconds; dividing by", den, "\n")
    dt[, time_s := as.numeric(time_s) / den]
  } else {
    cat("[TIME]", label, "time_s already ~[0,1]\n")
    dt[, time_s := as.numeric(time_s)]
  }
  
  dt[, time_s := pmin(pmax(time_s, 0), 1)]
  dt
}

dtA <- normalize_time_s_to_01(dtA, day_max = DAY_MAX, label = "AUMC")
dtG <- normalize_time_s_to_01(dtG, day_max = DAY_MAX, label = "GENEVE")

if (!("time_day" %in% names(dtA)) && ("time_s" %in% names(dtA))) dtA[, time_day := time_s * DAY_MAX]
if (!("time_day" %in% names(dtG)) && ("time_s" %in% names(dtG))) dtG[, time_day := time_s * DAY_MAX]

if (!("time_use" %in% names(dtA))) {
  if ("time_s" %in% names(dtA)) dtA[, time_use := time_s] else
    if ("time_day" %in% names(dtA)) dtA[, time_use := time_day] else
      stop("AUMC: cannot create time_use (no time_s/time_day).")
}

if (!("time_use" %in% names(dtG))) {
  if ("time_s" %in% names(dtG)) dtG[, time_use := time_s] else
    if ("time_day" %in% names(dtG)) dtG[, time_use := time_day] else
      stop("GENEVE: cannot create time_use (no time_s/time_day).")
}

## -------------------------
## 3) Load AUMC-trained models
## -------------------------
rds_all <- list.files(MODELS_DIR, pattern = "\\.rds$", full.names = TRUE)
if (!is.null(MODEL_PATTERN)) rds_all <- rds_all[grepl(MODEL_PATTERN, basename(rds_all))]

fits <- list()
for (f in rds_all) {
  obj <- try(readRDS(f), silent = TRUE)
  if (inherits(obj, "try-error")) next
  if (!inherits(obj, "Jointlcmm")) next
  nm <- tools::file_path_sans_ext(basename(f))
  fits[[nm]] <- obj
}
stopifnot(length(fits) >= 1)

fits <- fits[order(vapply(fits, function(x) x$ng, numeric(1)))]
cat("[MODELS] loaded:", paste(sprintf("%s(G=%d)", names(fits), vapply(fits, `[[`, numeric(1), "ng")), collapse = ", "), "\n")

## -------------------------
## 4) Detect TIMEVAR on GENEVE via predictY()
## -------------------------
detect_timevar_by_predictY <- function(fit, dt,
                                       candidates = c("time_use","time_s","time_day","t_hours","time")) {
  candidates <- intersect(candidates, names(dt))
  if (!length(candidates)) stop("No time candidates in dt. Columns: ", paste(names(dt), collapse = ", "))
  
  for (tv in candidates) {
    x <- dt[[tv]]
    x <- x[is.finite(x)]
    if (length(x) < 10) next
    tgrid <- sort(unique(quantile(x, probs = c(0.05,0.25,0.5,0.75,0.95), na.rm = TRUE)))
    newd <- data.frame(tmp = tgrid); names(newd) <- tv
    ok <- try(suppressWarnings(predictY(fit, newdata = newd, var.time = tv, draws = FALSE)), silent = TRUE)
    if (!inherits(ok, "try-error")) return(tv)
  }
  
  stop("TIMEVAR detection failed by predictY() on candidates: ", paste(candidates, collapse = ", "))
}

ensure_timevar <- function(dt, timevar) {
  if (timevar %in% names(dt)) return(dt)
  
  if (timevar == "time_use") {
    if ("time_s" %in% names(dt)) { dt[, time_use := time_s]; return(dt) }
    if ("time_day" %in% names(dt)) { dt[, time_use := time_day]; return(dt) }
  }
  
  stop("Cannot build missing timevar=", timevar, " from columns: ", paste(names(dt), collapse = ", "))
}

TIMEVAR <- detect_timevar_by_predictY(fits[[1]], dtG)
cat("[TIMEVAR] detected:", TIMEVAR, "\n")

dtG <- ensure_timevar(dtG, TIMEVAR)
dtA <- ensure_timevar(dtA, TIMEVAR)
stopifnot(TIMEVAR %in% names(dtA), TIMEVAR %in% names(dtG))

## -------------------------
## 5) Build standard long tables for plots
## -------------------------
datG_long <- dtG[
  is.finite(scr_w99) & is.finite(get(TIMEVAR)) & is.finite(tte_pos) & !is.na(event2),
  .(id = as.integer(id),
    time = as.numeric(get(TIMEVAR)),
    scr_w99 = as.numeric(scr_w99),
    tte_pos = as.numeric(tte_pos),
    event2  = as.integer(event2))
] %>% as.data.frame()

datA_long <- dtA[
  is.finite(scr_w99) & is.finite(get(TIMEVAR)) & is.finite(tte_pos) & !is.na(event2),
  .(id = as.integer(id),
    time = as.numeric(get(TIMEVAR)),
    scr_w99 = as.numeric(scr_w99),
    tte_pos = as.numeric(tte_pos),
    event2  = as.integer(event2))
] %>% as.data.frame()

t_min <- min(datG_long$time, na.rm = TRUE)
t_max <- max(datG_long$time, na.rm = TRUE)
time_grid <- seq(t_min, t_max, length.out = 60)

to_day <- function(time_val, ref_time_vec) {
  if (TIMEVAR == "time_day") return(time_val)
  
  if (TIMEVAR %in% c("time_s","time_use")) {
    return(time_val * DAY_MAX)
  }
  
  rng <- range(ref_time_vec, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[2] <= rng[1]) return(time_val)
  (time_val - rng[1]) / (rng[2] - rng[1]) * DAY_MAX
}

## -------------------------
## 6) Core helpers
## -------------------------
predict_pp_external <- function(fit, dt, timevar) {
  needed <- c("id", timevar, "scr_w99", "tte_pos", "event2")
  miss <- setdiff(needed, names(dt))
  if (length(miss)) stop("Missing columns for predictClass: ", paste(miss, collapse = ", "))
  
  newd <- as.data.frame(dt[, ..needed])
  pp <- suppressWarnings(predictClass(fit, newdata = newd, subject = "id")) %>% as.data.frame()
  pp$id <- as.integer(pp$id)
  pp$class <- as.integer(pp$class)
  names(pp) <- sub("^probYT", "prob", names(pp))
  pp
}

heatmap_mat <- function(pp, G) {
  pcols <- grep("^prob\\d+$", names(pp), value = TRUE)
  stopifnot(length(pcols) >= G)
  pcols <- pcols[1:G]
  
  M <- do.call(rbind, lapply(1:G, function(k) {
    sub <- pp[pp$class == k, pcols, drop = FALSE]
    if (!nrow(sub)) return(rep(NA_real_, G))
    colMeans(sub, na.rm = TRUE)
  }))
  colnames(M) <- paste0("target_", 1:G)
  
  data.frame(assigned = 1:G, M, check.names = FALSE) %>%
    pivot_longer(-assigned, names_to = "target", values_to = "pmean") %>%
    mutate(target = as.integer(gsub("target_", "", target)))
}

plot_heat <- function(df, G, title) {
  ggplot(df, aes(x = target, y = assigned, fill = pmean)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(is.na(pmean), "", sprintf("%.2f", pmean))), size = 3.2) +
    scale_x_continuous(breaks = 1:G) +
    scale_y_continuous(breaks = 1:G) +
    coord_equal() +
    labs(title = title, x = "Target class j", y = "Assigned class k", fill = "P̄") +
    theme_minimal(base_size = 12)
}

predicted_traj_df <- function(fit, time_grid, ref_time_vec) {
  G <- fit$ng
  newd <- data.frame(tmp = time_grid); names(newd) <- TIMEVAR
  
  pr  <- predictY(fit, newdata = newd, var.time = TIMEVAR, draws = FALSE)
  mat <- as.data.frame(pr$pred)
  cn  <- colnames(mat)
  
  keep <- grep("^Ypred_class[0-9]+$", cn, value = TRUE)
  if (length(keep) < G) keep <- grep("^Ypred_[0-9]+$", cn, value = TRUE)
  if (length(keep) < G) {
    num_cols <- vapply(mat, is.numeric, TRUE)
    keep <- cn[num_cols][seq_len(G)]
  }
  keep <- keep[seq_len(G)]
  
  df <- data.frame(time_val = time_grid, mat[, keep, drop = FALSE], check.names = FALSE) %>%
    pivot_longer(cols = all_of(keep), names_to = "ClassRaw", values_to = "scr_pred")
  
  df$Class <- factor(match(df$ClassRaw, keep), levels = 1:G, labels = paste0("Class ", 1:G))
  df$day   <- to_day(df$time_val, ref_time_vec)
  df[, c("day","Class","scr_pred")]
}

plot_traj_obs_vs_pred_external <- function(fit, dat_long, pp, title) {
  G <- fit$ng
  pp_h <- pp %>% transmute(id = as.integer(id), class = as.integer(class))
  
  long_df <- dat_long %>%
    inner_join(pp_h, by = "id") %>%
    mutate(day = floor(to_day(time, dat_long$time))) %>%
    filter(day >= 0, day <= DAY_MAX, is.finite(scr_w99), !is.na(class))
  
  obs <- long_df %>%
    group_by(class, day) %>%
    summarise(med = median(scr_w99), .groups = "drop") %>%
    mutate(Class = factor(paste0("Class ", class), levels = paste0("Class ", 1:G)))
  
  pred <- predicted_traj_df(fit, time_grid, dat_long$time)
  
  ggplot() +
    geom_line(data = pred, aes(x = day, y = scr_pred, colour = Class), linewidth = 1.1, linetype = 2) +
    geom_line(data = obs,  aes(x = day, y = med, colour = Class), linewidth = 1.1) +
    facet_wrap(~ Class, ncol = 2) +
    scale_x_continuous("Days", breaks = 0:DAY_MAX, limits = c(0, DAY_MAX)) +
    ylab("SCr (mg/dL, winsorized using q99 AUMC)") +
    coord_cartesian(ylim = Y_LIM) +
    labs(title = title, colour = "Class") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
}

to_day_dt <- function(time_val, dt_ref, timevar, day_max = 14) {
  if (timevar == "time_day") return(time_val)
  if (timevar %in% c("time_s","time_use")) return(time_val * day_max)
  
  rng <- range(dt_ref[[timevar]], na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[2] <= rng[1]) return(time_val)
  (time_val - rng[1]) / (rng[2] - rng[1]) * day_max
}

traj_compare_plot <- function(dtA, dtG, ppA, ppG, timevar, day_max = 14,
                              y_lim = c(0, 5),
                              title = "SCr(w99) trajectories — AUMC (solid) vs GENEVE (dashed)") {
  
  ppA_h <- ppA %>% transmute(id = as.integer(id), class = as.integer(class))
  ppG_h <- ppG %>% transmute(id = as.integer(id), class = as.integer(class))
  
  dA <- as.data.frame(dtA[, .(id = as.integer(id),
                              time = as.numeric(get(timevar)),
                              scr_w99 = as.numeric(scr_w99))]) %>%
    inner_join(ppA_h, by = "id") %>%
    mutate(dataset = "AUMC",
           day = floor(to_day_dt(time, dtA, timevar, day_max)))
  
  dG <- as.data.frame(dtG[, .(id = as.integer(id),
                              time = as.numeric(get(timevar)),
                              scr_w99 = as.numeric(scr_w99))]) %>%
    inner_join(ppG_h, by = "id") %>%
    mutate(dataset = "GENEVE",
           day = floor(to_day_dt(time, dtG, timevar, day_max)))
  
  df <- bind_rows(dA, dG) %>%
    filter(is.finite(scr_w99), is.finite(day), day >= 0, day <= day_max, !is.na(class))
  
  G <- max(df$class, na.rm = TRUE)
  
  stats <- df %>%
    group_by(dataset, class, day) %>%
    summarise(
      med = median(scr_w99, na.rm = TRUE),
      q1  = quantile(scr_w99, 0.25, na.rm = TRUE, type = 2),
      q3  = quantile(scr_w99, 0.75, na.rm = TRUE, type = 2),
      .groups = "drop"
    ) %>%
    tidyr::complete(dataset, class, day = 0:day_max) %>%
    mutate(Class = factor(paste0("Class ", class), levels = paste0("Class ", 1:G)))
  
  stats_A  <- stats %>% filter(dataset == "AUMC", is.finite(med))
  stats_G  <- stats %>% filter(dataset == "GENEVE", is.finite(med))
  ribbon_G <- stats %>% filter(dataset == "GENEVE", is.finite(q1), is.finite(q3))
  
  ggplot() +
    geom_ribbon(
      data = ribbon_G,
      aes(x = day, ymin = q1, ymax = q3, fill = Class, group = Class),
      alpha = 0.15, colour = NA, show.legend = FALSE
    ) +
    geom_line(
      data = stats_A,
      aes(x = day, y = med, colour = Class, group = Class),
      linewidth = 1.1, linetype = "solid"
    ) +
    geom_line(
      data = stats_G,
      aes(x = day, y = med, colour = Class, group = Class),
      linewidth = 1.1, linetype = "dashed"
    ) +
    facet_wrap(~ Class, ncol = 2) +
    scale_x_continuous("Days", breaks = 0:day_max, limits = c(0, day_max)) +
    ylab("SCr (mg/dL, winsorized using q99 AUMC)") +
    coord_cartesian(ylim = y_lim) +
    labs(title = title, colour = "Class") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}

cif_onecause <- function(sub, cause) {
  ci <- cmprsk::cuminc(ftime = sub$tte_pos, fstatus = sub$event2, cencode = 0)
  nms <- names(ci)
  
  key <- NULL
  if (as.character(cause) %in% nms) key <- as.character(cause)
  if (is.null(key)) {
    cand <- grep(paste0("^", cause, "(\\s|$)"), nms, value = TRUE)
    if (length(cand)) key <- cand[1]
  }
  if (is.null(key)) return(NULL)
  
  obj <- ci[[key]]
  if (is.null(obj$time) || is.null(obj$est)) return(NULL)
  tibble::tibble(time = obj$time, cif = obj$est)
}

cif_by_class_both <- function(dt, pp_hard, G) {
  df <- as.data.frame(unique(dt[, .(id, tte_pos, event2)])) %>%
    inner_join(pp_hard, by = "id") %>%
    filter(is.finite(tte_pos), is.finite(event2), !is.na(class_hard))
  
  out <- lapply(1:G, function(k) {
    sub <- df[df$class_hard == k, c("tte_pos","event2")]
    d1 <- cif_onecause(sub, 1)
    d2 <- cif_onecause(sub, 2)
    
    if (!is.null(d1)) d1 <- mutate(d1, Class = paste0("Class ", k), Event = "ICU death")
    if (!is.null(d2)) d2 <- mutate(d2, Class = paste0("Class ", k), Event = "ICU discharge alive")
    
    bind_rows(d1, d2)
  })
  
  bind_rows(out) %>%
    mutate(
      Class = factor(Class, levels = paste0("Class ", 1:G)),
      Event = factor(Event, levels = c("ICU death", "ICU discharge alive"))
    )
}

plot_cif_geneve <- function(dtG, ppG, G, title) {
  pp_h <- ppG %>% transmute(id = as.integer(id), class_hard = as.integer(class))
  df <- cif_by_class_both(dtG, pp_h, G) %>% filter(is.finite(time), is.finite(cif), time >= 0, time <= DAY_MAX)
  
  ggplot(df, aes(time, cif, colour = Class)) +
    geom_line(linewidth = 1.1) +
    facet_wrap(~ Event, nrow = 1) +
    scale_x_continuous("Days", limits = c(0, DAY_MAX)) +
    scale_y_continuous("Cumulative incidence", limits = c(0, 1)) +
    labs(title = title, colour = "Class") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}

plot_cif_compare <- function(dtA, dtG, ppA, ppG, G, admin = 14,
                             title = "Empirical CIF — AUMC (solid) vs GENEVE (dashed)") {
  
  ppA_h <- ppA %>% transmute(id = as.integer(id), class_hard = factor(as.integer(class), levels = 1:G)) %>% filter(!is.na(class_hard))
  ppG_h <- ppG %>% transmute(id = as.integer(id), class_hard = factor(as.integer(class), levels = 1:G)) %>% filter(!is.na(class_hard))
  
  cifA <- cif_by_class_both(dtA, ppA_h, G) %>% mutate(dataset = "AUMC")
  cifG <- cif_by_class_both(dtG, ppG_h, G) %>% mutate(dataset = "GENEVE")
  
  df <- bind_rows(cifA, cifG) %>%
    filter(is.finite(time), is.finite(cif), time >= 0, time <= admin)
  
  ggplot(df, aes(x = time, y = cif, colour = Class, linetype = dataset)) +
    geom_line(linewidth = 1.1) +
    facet_grid(Event ~ Class) +
    scale_linetype_manual(values = c(AUMC = "solid", GENEVE = "dashed")) +
    scale_x_continuous("Days", limits = c(0, admin)) +
    scale_y_continuous("Cumulative incidence", limits = c(0, 1)) +
    labs(title = title, colour = "Class", linetype = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}

post_summary <- function(pp, G) {
  pcols <- paste0("prob", 1:G)
  stopifnot(all(pcols %in% names(pp)))
  
  P <- as.matrix(pp[, pcols, drop = FALSE])
  pmax <- apply(P, 1, max, na.rm = TRUE)
  
  if (G == 1) {
    ent_mean <- NA_real_
  } else {
    eps <- 1e-15
    P2 <- pmin(pmax(P, eps), 1 - eps)
    ent <- 1 - (-rowSums(P2 * log(P2)) / log(G))
    ent_mean <- mean(ent, na.rm = TRUE)
  }
  
  data.frame(
    n = nrow(pp),
    Mean_posterior_probability = mean(pmax, na.rm = TRUE),
    P_ge_0.70 = mean(pmax >= 0.70, na.rm = TRUE) * 100,
    P_ge_0.90 = mean(pmax >= 0.90, na.rm = TRUE) * 100,
    Mean_entropy = ent_mean
  )
}

## -------------------------
## 7) RUN over all models
## -------------------------
summary_tab <- list()

for (nm in names(fits)) {
  fit <- fits[[nm]]
  G   <- fit$ng
  cat("\n=== ", nm, " | G=", G, " ===\n", sep = "")
  
  ## Predict on GENEVE
  ppG <- predict_pp_external(fit, dtG, TIMEVAR)
  write_csv(ppG, file.path(OUT_DIR, "tables", paste0("ppGENEVE_", nm, ".csv")))
  
  ## Predict on AUMC
  ppA <- predict_pp_external(fit, dtA, TIMEVAR)
  write_csv(ppA, file.path(OUT_DIR, "tables", paste0("ppAUMC_", nm, ".csv")))
  
  ## Posterior uncertainty
  sG <- post_summary(ppG, G)
  sA <- post_summary(ppA, G)
  
  summ_unc <- dplyr::bind_rows(
    dplyr::mutate(as.data.frame(sG), dataset = "GENEVE"),
    dplyr::mutate(as.data.frame(sA), dataset = "AUMC")
  ) %>%
    dplyr::mutate(
      model = nm,
      G = G,
      Mean_posterior_probability = round(Mean_posterior_probability, 3),
      P_ge_0.70 = round(P_ge_0.70, 1),
      P_ge_0.90 = round(P_ge_0.90, 1),
      Mean_entropy = round(Mean_entropy, 3)
    )
  
  readr::write_csv(summ_unc, file.path(OUT_DIR, "tables", paste0("posterior_summary_", nm, ".csv")))
  
  ## Class sizes (GENEVE)
  tab <- as.data.frame(table(ppG$class))
  names(tab) <- c("class", "n")
  tab$class <- as.integer(as.character(tab$class))
  tab$pct <- round(100 * tab$n / sum(tab$n), 2)
  write_csv(tab, file.path(OUT_DIR, "tables", paste0("class_sizes_GENEVE_", nm, ".csv")))
  print(tab)
  
  ## Heatmap (GENEVE)
  hm <- heatmap_mat(ppG, G)
  ggsave(file.path(OUT_DIR, "figures/heatmap", paste0("heatmap_GENEVE_", nm, ".png")),
         plot_heat(hm, G, paste0("GENEVE under AUMC model — mean posterior probs (", nm, ")")),
         width = 7, height = 6, dpi = 200)
  
  ## Heatmap (AUMC)
  hmA <- heatmap_mat(ppA, G)
  ggsave(file.path(OUT_DIR, "figures/heatmap", paste0("heatmap_AUMC_", nm, ".png")),
         plot_heat(hmA, G, paste0("AUMC under same model — mean posterior probs (", nm, ")")),
         width = 7, height = 6, dpi = 200)
  
  ## Traj obs vs pred (GENEVE)
  p_tr <- plot_traj_obs_vs_pred_external(
    fit, datG_long, ppG,
    title = paste0("GENEVE — observed (median) vs predicted SCr — ", nm)
  )
  ggsave(file.path(OUT_DIR, "figures/traj", paste0("traj_obs_pred_GENEVE_", nm, ".png")),
         p_tr, width = 10, height = 7, dpi = 200)
  
  ## Traj compare AUMC vs GENEVE
  p_tr_cmp <- traj_compare_plot(
    dtA = dtA, dtG = dtG,
    ppA = ppA, ppG = ppG,
    timevar = TIMEVAR,
    day_max = DAY_MAX,
    y_lim = Y_LIM,
    title = paste0("SCr(w99) — AUMC (solid) vs GENEVE (dashed) — ", nm)
  )
  ggsave(file.path(OUT_DIR, "figures/traj", paste0("traj_AUMC_vs_GENEVE_", nm, ".png")),
         p_tr_cmp, width = 10, height = 7, dpi = 200)
  
  ## CIF empirical on GENEVE
  p_cif <- plot_cif_geneve(dtG, ppG, G, paste0("GENEVE — empirical CIF by predicted class — ", nm))
  ggsave(file.path(OUT_DIR, "figures/cif", paste0("cif_GENEVE_", nm, ".png")),
         p_cif, width = 10, height = 6, dpi = 200)
  
  ## CIF compare AUMC vs GENEVE
  p_cif_cmp <- plot_cif_compare(
    dtA = dtA, dtG = dtG, ppA = ppA, ppG = ppG, G = G, admin = DAY_MAX,
    title = paste0("Empirical CIF — AUMC (solid) vs GENEVE (dashed) — ", nm)
  )
  ggsave(file.path(OUT_DIR, "figures/cif_compare", paste0("cif_AUMC_vs_GENEVE_", nm, ".png")),
         p_cif_cmp, width = 14, height = 7.5, dpi = 200)
  
  summary_tab[[nm]] <- data.frame(
    model = nm, G = G,
    n_patients = length(unique(ppG$id)),
    min_class_n = min(tab$n),
    min_class_pct = min(tab$pct),
    stringsAsFactors = FALSE
  )
}

summ <- bind_rows(summary_tab)
write_csv(summ, file.path(OUT_DIR, "tables", "summary_external_GENEVE.csv"))

unc_files <- list.files(file.path(OUT_DIR, "tables"), pattern = "^posterior_summary_.*\\.csv$", full.names = TRUE)
unc_all <- dplyr::bind_rows(lapply(unc_files, readr::read_csv, show_col_types = FALSE))
readr::write_csv(unc_all, file.path(OUT_DIR, "tables", "posterior_summary_ALLmodels.csv"))

cat("\n[SAVED] ", OUT_DIR, "\n", sep = "")
cat("✓ External validation done. TIMEVAR=", TIMEVAR, " | DAY_MAX=", DAY_MAX, "\n", sep = "")

## ============================================================
## EXPORT MAPPING old_id / new_id / class for G3 model on GENEVE
## ============================================================

fitG3_name <- "AUMC_joint_G3_AUMC_SCr_H14_quad_w99_cr_H14d_time_s"

## 1) lire les posterior probs déjà exportées
ppG3_file <- file.path(OUT_DIR, "tables", paste0("ppGENEVE_", fitG3_name, ".csv"))
ppG3 <- readr::read_csv(ppG3_file, show_col_types = FALSE)

## 2) table de correspondance old_id / new_id depuis dtG
mapG <- unique(as.data.frame(dtG[, .(new_id = id, old_id = id_chr)]))
mapG$new_id <- as.integer(mapG$new_id)

stopifnot(!anyDuplicated(mapG$new_id))
stopifnot(!anyDuplicated(mapG$old_id))

## 3) harmoniser le nom de l'id dans ppG3
stopifnot("id" %in% names(ppG3))
ppG3$new_id <- as.integer(ppG3$id)

## 4) joindre les classes aux anciens IDs
ppG3_map <- ppG3 %>%
  left_join(mapG, by = "new_id") %>%
  relocate(old_id, new_id, class)

## 5) version patient-level unique
ppG3_map_patient <- ppG3_map %>%
  distinct(old_id, new_id, .keep_all = TRUE)

map_class_G3 <- ppG3_map_patient %>%
  dplyr::select(old_id, new_id, class)

## 6) exports
readr::write_csv(
  ppG3_map,
  file.path(OUT_DIR, "tables", paste0("ppGENEVE_G3_with_old_new_id_", fitG3_name, ".csv"))
)

readr::write_csv(
  map_class_G3,
  file.path(OUT_DIR, "tables", paste0("mapping_oldid_newid_class_G3_", fitG3_name, ".csv"))
)

## 7) checks
cat("N rows ppG3_map: ", nrow(ppG3_map), "\n", sep = "")
cat("N unique old_id: ", dplyr::n_distinct(ppG3_map$old_id), "\n", sep = "")
cat("N unique new_id: ", dplyr::n_distinct(ppG3_map$new_id), "\n", sep = "")
cat("Missing old_id: ", sum(is.na(ppG3_map$old_id)), "\n", sep = "")

print(table(map_class_G3$class, useNA = "ifany"))
print(head(map_class_G3))
