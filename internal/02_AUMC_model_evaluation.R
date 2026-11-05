## ============================================================
## EVALUATION — Jointlcmm quadratic + CR (AUMC/MIMIC/...)
## Robuste: auto-détection time variable + auto-load models
## Outputs: metrics + elbow + traj + obs vs pred + CIF
## ============================================================

rm(list = ls())
graphics.off()

suppressPackageStartupMessages({
  library(lcmm)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(readr)
  library(stringr)
})

## -------------------------
## CONFIG (à adapter seulement ici)
## -------------------------
DB_RDS     <- Sys.getenv("AKI_AUMC_RDS")
MODELS_DIR <- Sys.getenv("AKI_MODELS_DIR")

## Horizon clinique (pour les plots en jours). Si tu ne sais pas: laisse NA et on infère depuis time_s si possible.
DAY_MAX <- 14

## Pattern optionnel pour ne charger qu’un sous-ensemble de modèles
## ex: "VENT_noCKD_H14_quad" ou "MIMIC" etc. Laisse NULL pour tout.
MODEL_PATTERN <- if (nzchar(Sys.getenv("AKI_MODEL_PATTERN"))) Sys.getenv("AKI_MODEL_PATTERN") else NULL

## -------------------------
## 0. Dossiers & données
## -------------------------
setwd(MODELS_DIR)

dir.create("figures",          showWarnings = FALSE)
dir.create("figures/traj",     showWarnings = FALSE)
dir.create("figures/cif",      showWarnings = FALSE)
dir.create("figures/elbow",    showWarnings = FALSE)
dir.create("tables",           showWarnings = FALSE)

dt0 <- as.data.table(readRDS(DB_RDS))

## Winsorisation SCr à 99% (si scr_w99 absent ou si tu veux recalculer proprement)
if (!("scr_w99" %in% names(dt0))) {
  stopifnot("scr_mgdl" %in% names(dt0))
  q99 <- quantile(dt0$scr_mgdl, probs = 0.99, na.rm = TRUE)
  dt0[, scr_w99 := pmin(scr_mgdl, as.numeric(q99))]
} else {
  ## option: refaire quand même pour être sûr que c'est cohérent
  q99 <- quantile(dt0$scr_mgdl, probs = 0.99, na.rm = TRUE)
  dt0[, scr_w99 := pmin(scr_mgdl, as.numeric(q99))]
}

stopifnot(all(c("id","tte_pos","event2","scr_w99") %in% names(dt0)))

## -------------------------
## 1. Chargement auto des modèles
## -------------------------
rds_all <- list.files(MODELS_DIR, pattern = "\\.rds$", full.names = TRUE)

if (!is.null(MODEL_PATTERN)) {
  rds_all <- rds_all[grepl(MODEL_PATTERN, basename(rds_all))]
}

fits <- list()
for (f in rds_all) {
  obj <- try(readRDS(f), silent = TRUE)
  if (inherits(obj, "try-error")) next
  if (!inherits(obj, "Jointlcmm")) next
  nm <- tools::file_path_sans_ext(basename(f))
  fits[[nm]] <- obj
}

stopifnot(length(fits) >= 1)

if (!("time_use" %in% names(dt0))) {
  if ("time_s" %in% names(dt0)) {
    dt0[, time_use := time_s]
  } else if ("time_day" %in% names(dt0)) {
    dt0[, time_use := time_day]
  }
}


## Tri par ng
fits <- fits[order(vapply(fits, function(x) x$ng, numeric(1)))]
message("Modèles chargés : ", paste(sprintf("%s(G=%d)", names(fits), vapply(fits, `[[`, numeric(1), "ng")), collapse = ", "))


## -------------------------
## 2. Détection ROBUSTE de la variable temps du modèle
##   -> on teste predictY() sur plusieurs candidats
## -------------------------
detect_timevar_by_predictY <- function(fit, dt0,
                                       candidates = c("time_use","time_s","time_day","t_hours","time")) {
  candidates <- intersect(candidates, names(dt0))
  if (!length(candidates)) {
    stop("Aucun candidat temps trouvé dans dt0. Colonnes dt0: ", paste(names(dt0), collapse = ", "))
  }
  
  ## mini-grille de temps (dans l’échelle de chaque candidat)
  for (tv in candidates) {
    x <- dt0[[tv]]
    x <- x[is.finite(x)]
    if (length(x) < 10) next
    
    tgrid <- sort(unique(quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)))
    newd  <- data.frame(tmp = tgrid)
    names(newd) <- tv
    
    ok <- try({
      suppressWarnings(
        lcmm::predictY(fit, newdata = newd, var.time = tv, draws = FALSE)
      )
    }, silent = TRUE)
    
    if (!inherits(ok, "try-error")) return(tv)
  }
  
  ## Si rien ne marche, on affiche le message d'erreur le plus informatif
  last_err <- NULL
  for (tv in candidates) {
    x <- dt0[[tv]]
    x <- x[is.finite(x)]
    if (length(x) < 10) next
    tgrid <- sort(unique(quantile(x, probs = c(0.1, 0.5, 0.9), na.rm = TRUE)))
    newd  <- data.frame(tmp = tgrid); names(newd) <- tv
    last_err <- tryCatch({
      lcmm::predictY(fit, newdata = newd, var.time = tv, draws = FALSE)
      NULL
    }, error = function(e) e$message)
    if (is.null(last_err)) return(tv)
  }
  
  stop("Impossible de détecter la variable temps via predictY(). Dernière erreur: ", last_err)
}

TIMEVAR <- detect_timevar_by_predictY(fits[[1]], dt0)
message("[TIMEVAR] détecté par predictY(): ", TIMEVAR)
stopifnot(TIMEVAR %in% names(dt0))


## Dataset pour éval / plots
dat <- dt0[is.finite(scr_w99) & is.finite(get(TIMEVAR)) & is.finite(tte_pos) & !is.na(event2),
           .(id, time = get(TIMEVAR), scr_w99, tte_pos, event2)]
dat <- as.data.frame(dat)

## -------------------------
## 3. Conversion en jours pour l'axe X
## -------------------------
## On privilégie time_day si présent (c'est le "vrai" jour ICU).
has_time_day <- "time_day" %in% names(dt0)
has_time_s   <- "time_s" %in% names(dt0)

infer_daymax <- function(dt0, default = DAY_MAX) {
  if (is.finite(default)) return(default)
  if (has_time_s) {
    mx <- max(dt0$time_s, na.rm = TRUE)
    ## si c'est un time_s normalisé, max ~ 1
    if (mx <= 1.05) return(14)  # fallback raisonnable
  }
  14
}

DAY_MAX_USE <- infer_daymax(dt0, DAY_MAX)

to_day <- function(time_val, dt0) {
  ## cas 1: on a time_day et TIMEVAR est time_day ou time_s/time_use
  if (has_time_day && TIMEVAR == "time_day") return(time_val)
  if (has_time_day && TIMEVAR %in% c("time_s","time_use")) {
    ## si time_s est dispo et ressemble à time_day/DAY_MAX
    if (has_time_s) {
      ## approx: days = time_s * DAY_MAX
      return(time_val * DAY_MAX_USE)
    }
  }
  ## fallback linéaire sur [min,max]
  tmin <- min(dt0[[TIMEVAR]], na.rm = TRUE)
  tmax <- max(dt0[[TIMEVAR]], na.rm = TRUE)
  if (!is.finite(tmin) || !is.finite(tmax) || tmax <= tmin) return(time_val)
  (time_val - tmin) / (tmax - tmin) * DAY_MAX_USE
}

## time grid dans l'unité du modèle
t_min <- min(dat$time, na.rm = TRUE)
t_max <- max(dat$time, na.rm = TRUE)
time_grid <- seq(t_min, t_max, length.out = 60)

## -------------------------
## 4. Fonctions utilitaires (pprob, metrics, plots)
## -------------------------
pprob_std <- function(fit) {
  G <- fit$ng
  pp <- try({
    if (!is.null(fit$pprob)) fit$pprob else {
      tmp <- lcmm::postprob(fit)
      if (is.list(tmp) && !is.null(tmp$pprob)) tmp$pprob else as.data.frame(tmp)
    }
  }, silent = TRUE)
  
  if (inherits(pp, "try-error") || !is.data.frame(pp)) stop("Impossible d’extraire pprob pour G=", G)
  
  if (!("id" %in% names(pp))) {
    if ("subject" %in% names(pp)) names(pp)[names(pp) == "subject"] <- "id"
    if (!("id" %in% names(pp))) stop("Pas de colonne id/subject dans pprob (G=", G, ")")
  }
  
  cnc <- gsub("\\s+","", tolower(names(pp)))
  idx <- grep("^prob\\d+$", cnc)
  
  if (length(idx) < G) {
    cand <- which(vapply(pp, is.numeric, TRUE))
    cand <- cand[vapply(pp[cand], function(x){
      r <- range(x, na.rm = TRUE)
      is.finite(r[1]) && is.finite(r[2]) && r[1] >= -1e-9 && r[2] <= 1 + 1e-9
    }, TRUE)]
    stopifnot(length(cand) >= G)
    idx <- cand[seq_len(G)]
  } else {
    ord <- order(as.integer(gsub("^prob(\\d+)$","\\1", cnc[idx])))
    idx <- idx[ord][seq_len(G)]
  }
  
  names(pp)[idx] <- paste0("prob", 1:G)
  
  if (!("class" %in% names(pp))) {
    P <- as.matrix(pp[, paste0("prob", 1:G), drop = FALSE])
    pp$class <- max.col(P, ties.method = "first")
  }
  
  pp[, c("id", paste0("prob", 1:G), "class")]
}

class_counts <- function(fit) {
  pp <- pprob_std(fit)
  G  <- fit$ng
  hard <- factor(pp$class, levels = 1:G)
  tab  <- as.data.frame(table(hard))
  n    <- sum(tab$Freq)
  tab %>% transmute(G = G, Class = as.integer(hard), n = as.integer(Freq), pct = round(100 * Freq / n, 2))
}

cols_for_G <- function(G, lvl = NULL) {
  cols <- hcl.colors(G, "Dark3")
  if (!is.null(lvl)) setNames(cols, lvl) else cols
}

metrics_from_fit <- function(fit) {
  G  <- fit$ng
  pp <- pprob_std(fit)
  n  <- dplyr::n_distinct(pp$id)
  P  <- as.matrix(pp[, paste0("prob", 1:G), drop = FALSE])
  eps <- 1e-12
  
  entropy_total <- -sum(P * log(pmax(P, eps)))
  entropy_norm  <- entropy_total / (n * log(G))
  separation    <- 1 - entropy_norm
  
  hard <- pp$class
  tab  <- table(base::factor(hard, levels = 1:G))
  min_n   <- as.integer(min(tab))
  min_pct <- 100 * min_n / n
  
  max_post <- apply(P, 1, max, na.rm = TRUE)
  
  BIC_val <- suppressWarnings(if (!is.null(fit$BIC)) as.numeric(fit$BIC) else NA_real_)
  ICL_val <- suppressWarnings(if (!is.null(fit$ICL)) as.numeric(fit$ICL) else NA_real_)
  if (is.na(ICL_val) && is.finite(BIC_val)) {
    ICL_val <- BIC_val + 2 * entropy_total
  }
  
  logLik_val <- suppressWarnings(if (!is.null(fit$loglik)) as.numeric(fit$loglik) else NA_real_)
  
  data.frame(
    G = G,
    n = n,
    logLik = logLik_val,
    BIC = BIC_val,
    ICL = ICL_val,
    Entropy_total = entropy_total,
    Entropy_norm  = entropy_norm,
    Separation    = separation,
    MinClassN     = min_n,
    MinClassPct   = round(min_pct, 2),
    MeanMaxPost   = round(mean(max_post, na.rm = TRUE), 4),
    MedianMaxPost = round(median(max_post, na.rm = TRUE), 4),
    stringsAsFactors = FALSE
  )
}

class_metrics_from_fit <- function(fit) {
  G  <- fit$ng
  pp <- pprob_std(fit)
  n  <- dplyr::n_distinct(pp$id)
  
  out <- lapply(1:G, function(k){
    sub <- pp %>% filter(class == k)
    nk  <- nrow(sub)
    pct <- 100 * nk / n
    pcols <- paste0("prob", 1:G)
    diag_mean <- if (nk > 0) mean(sub[[paste0("prob", k)]], na.rm = TRUE) else NA_real_
    diag_med  <- if (nk > 0) median(sub[[paste0("prob", k)]], na.rm = TRUE) else NA_real_
    max_post  <- if (nk > 0) apply(as.matrix(sub[, pcols, drop = FALSE]), 1, max, na.rm = TRUE) else numeric()
    data.frame(
      G = G, Class = k,
      n_class = nk, Pct_class = round(pct, 2),
      DiagonalProb_mean = round(diag_mean, 4),
      DiagonalProb_median = round(diag_med, 4),
      MaxPost_mean = if (length(max_post)) round(mean(max_post), 4) else NA_real_,
      MaxPost_median = if (length(max_post)) round(median(max_post), 4) else NA_real_
    )
  })
  do.call(rbind, out)
}

plot_elbow <- function(df, chosen_G = NA) {
  df <- df %>% arrange(G)
  
  gline <- function(p, xlab, ylab) {
    p + geom_line(linewidth = 0.8) +
      geom_point(size = 2) +
      xlab(xlab) + ylab(ylab) +
      theme_minimal(base_size = 12) +
      theme(panel.grid.minor = element_blank())
  }
  
  pA <- gline(ggplot(df, aes(G, BIC)), "Number of classes (G)", "BIC") + ggtitle("A  Elbow — BIC")
  pB <- gline(ggplot(df, aes(G, ICL)), "Number of classes (G)", "ICL") + ggtitle("B  Elbow — ICL")
  pC <- gline(ggplot(df, aes(G, Separation)), "Number of classes (G)", "Separation") + ggtitle("C  Separation")
  pD <- gline(ggplot(df, aes(G, MinClassPct)), "Number of classes (G)", "Min class size (%)") +
    geom_hline(yintercept = 5, linetype = "dotted") + ggtitle("D  Min class size")
  
  if (is.finite(chosen_G)) {
    add_v <- function(p) p + geom_vline(xintercept = chosen_G, linetype = "dashed")
    pA <- add_v(pA); pB <- add_v(pB); pC <- add_v(pC); pD <- add_v(pD)
  }
  
  cowplot::plot_grid(pA, pB, pC, pD, ncol = 2)
}

predicted_traj_df <- function(fit, time_grid) {
  G <- fit$ng
  newd <- data.frame(time = time_grid)
  names(newd) <- TIMEVAR
  
  pr  <- lcmm::predictY(fit, newdata = newd, var.time = TIMEVAR, draws = FALSE)
  mat <- as.data.frame(pr$pred)
  cn  <- colnames(mat)
  
  keep <- grep("^Ypred_class[0-9]+$", cn, value = TRUE)
  if (length(keep) < G) keep <- grep("^Ypred_[0-9]+$", cn, value = TRUE)
  if (length(keep) < G) {
    keep <- grep("pred", cn, value = TRUE)
    keep <- keep[!grepl("var|se|sd|low|up|inf|sup|ci", keep, ignore.case = TRUE)]
  }
  if (length(keep) < G) {
    num_cols <- vapply(mat, is.numeric, TRUE)
    stopifnot(sum(num_cols) >= G)
    keep <- cn[num_cols][seq_len(G)]
  }
  keep <- keep[seq_len(G)]
  
  df <- data.frame(time_val = time_grid, mat[, keep, drop = FALSE], check.names = FALSE) %>%
    tidyr::pivot_longer(cols = tidyselect::all_of(keep),
                        names_to = "ClassRaw", values_to = "scr_pred")
  
  df$Class <- factor(match(df$ClassRaw, keep),
                     levels = seq_len(G),
                     labels = paste0("Class ", seq_len(G)))
  df$time_day <- to_day(df$time_val, dt0)
  
  df[, c("time_day", "Class", "scr_pred")]
}

plot_traj_pred <- function(fit, time_grid, title_prefix = "MODEL") {
  df   <- predicted_traj_df(fit, time_grid)
  G    <- fit$ng
  cols <- cols_for_G(G, levels(df$Class))
  
  ggplot(df, aes(time_day, scr_pred, colour = Class)) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = cols, name = NULL) +
    labs(
      title = sprintf("%s — Predicted SCr trajectories (G=%d, quadratic)", title_prefix, G),
      x = sprintf("Days since reference (0–%d)", DAY_MAX_USE),
      y = "Predicted SCr (mg/dL, winsorized 99%)"
    ) +
    coord_cartesian(ylim = c(0, 5)) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "top",
          panel.grid.minor = element_blank())
}

plot_traj_obs_vs_pred <- function(fit, dat_long, day_break = 1, title_prefix = "MODEL") {
  G  <- fit$ng
  pp <- pprob_std(fit)
  hard <- pp %>% select(id, class)
  
  long_df <- dat_long %>%
    inner_join(hard, by = "id") %>%
    filter(is.finite(scr_w99), is.finite(time)) %>%
    mutate(time_day = to_day(time, dt0),
           time_bin = floor(time_day / day_break) * day_break)
  
  obs_summary <- long_df %>%
    group_by(class, time_bin) %>%
    summarise(med_scr = median(scr_w99), .groups = "drop") %>%
    mutate(Class = factor(paste0("Class ", class),
                          levels = paste0("Class ", 1:G)))
  
  df_pred <- predicted_traj_df(fit, time_grid)
  cols <- cols_for_G(G, levels(obs_summary$Class))
  
  ggplot() +
    geom_line(data = df_pred,
              aes(x = time_day, y = scr_pred, colour = Class),
              linewidth = 1.1, linetype = 2) +
    geom_line(data = obs_summary,
              aes(x = time_bin, y = med_scr, colour = Class),
              linewidth = 1.1) +
    facet_wrap(~ Class, ncol = ceiling(G / 2)) +
    scale_color_manual(values = cols, guide = "none") +
    scale_x_continuous("Days since reference",
                       breaks = seq(0, DAY_MAX_USE, by = 2),
                       limits = c(0, DAY_MAX_USE)) +
    ylab("SCr (mg/dL, winsorized 99%)") +
    labs(title = sprintf("%s — Observed (median) vs predicted SCr (G=%d)", title_prefix, G)) +
    coord_cartesian(ylim = c(0, 5)) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank())
}

cuminc_safe <- function(mdl, admin = 14, by = 1,
                        alt_bys = c(1, 0.5, 0.25, 0.1),
                        min_pts = 6, max_consec_fail = 10) {
  for (b in alt_bys) {
    tgrid <- seq(1e-6, admin, by = b)
    res <- try(lcmm::cuminc(mdl, time = tgrid), silent = TRUE)
    if (!inherits(res, "try-error")) return(list(ci = res, used_admin = max(tgrid), used_by = b, mode = "bulk"))
  }
  for (b in alt_bys) {
    tgrid <- seq(1e-6, admin, by = b)
    acc <- list(); used_t <- c(); ok_any <- FALSE; consec_fail <- 0
    for (tt in tgrid) {
      rr <- try(lcmm::cuminc(mdl, time = tt), silent = TRUE)
      if (inherits(rr, "try-error")) {
        consec_fail <- consec_fail + 1
        if (consec_fail >= max_consec_fail && length(used_t) >= min_pts) break
        next
      } else {
        ok_any <- TRUE; consec_fail <- 0
        acc <- c(acc, list(rr)); used_t <- c(used_t, tt)
      }
    }
    if (ok_any && length(used_t) >= min_pts) {
      n_events <- length(acc[[1]])
      merged <- vector("list", n_events)
      for (e in seq_len(n_events)) {
        merged[[e]] <- do.call(rbind, lapply(acc, function(x) as.data.frame(x[[e]])))
      }
      return(list(ci = merged, used_admin = max(used_t), used_by = b, mode = "piecemeal"))
    }
  }
  stop("cuminc() échoue malgré les tentatives.")
}

plot_cif <- function(mdl, admin = 14, title_prefix = "MODEL") {
  G <- mdl$ng
  out <- cuminc_safe(mdl, admin = admin, by = 1)
  ci <- out$ci
  
  dfs <- lapply(seq_along(ci), function(e) {
    df <- as.data.frame(ci[[e]])
    df$Event <- factor(df$event, levels = c(1, 2),
                       labels = c("ICU death", "ICU discharge alive"))
    df
  })
  
  df <- bind_rows(dfs) %>%
    pivot_longer(starts_with("class"), names_to = "Class", values_to = "CIF") %>%
    mutate(Class = factor(gsub("^class","Class ", Class),
                          levels = paste0("Class ", 1:G)))
  
  cols <- cols_for_G(G, levels(df$Class))
  
  ggplot(df, aes(time, CIF, colour = Class)) +
    geom_line(linewidth = 1.1) +
    facet_wrap(~ Event, nrow = 1) +
    scale_color_manual(values = cols, name = "Class") +
    scale_x_continuous("time (same unit as tte_pos)",
                       breaks = pretty(c(0, admin)),
                       limits = c(0, admin)) +
    scale_y_continuous("Cumulative incidence", limits = c(0, 1)) +
    ggtitle(sprintf("%s — CIF by class (G=%d)", title_prefix, G)) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank())
}

## -------------------------
## 5. RUN
## -------------------------
met_global  <- do.call(rbind, lapply(fits, metrics_from_fit)) %>% arrange(G)
met_byclass <- do.call(rbind, lapply(fits, class_metrics_from_fit)) %>% arrange(G, Class)

write_csv(met_global,  sprintf("tables/metrics_global_quadratic_%s.csv", basename(MODELS_DIR)))
write_csv(met_byclass, sprintf("tables/metrics_byclass_quadratic_%s.csv", basename(MODELS_DIR)))

print(met_global)
print(head(met_byclass, 10))

best_G_BIC <- met_global$G[which.min(met_global$BIC)]
p_elbow <- plot_elbow(met_global, chosen_G = best_G_BIC)
ggsave("figures/elbow/elbow_quadratic.png", p_elbow, width = 10, height = 7, dpi = 300)

## dat_long pour obs vs pred (il faut la colonne "time" standardisée)
dat_long <- dt0[is.finite(scr_w99) & is.finite(get(TIMEVAR)) & is.finite(tte_pos) & !is.na(event2),
                .(id, time = get(TIMEVAR), scr_w99, tte_pos, event2)] %>%
  as.data.frame()

for (nm in names(fits)) {
  fit <- fits[[nm]]
  G   <- fit$ng
  
  tab_classes <- class_counts(fit)
  write_csv(tab_classes, sprintf("tables/class_sizes_%s.csv", nm))
  
  p_pred <- plot_traj_pred(fit, time_grid, title_prefix = nm)
  ggsave(sprintf("figures/traj/SCr_pred_%s.png", nm), p_pred, width = 8, height = 5, dpi = 300)
  
  p_obs_pred <- plot_traj_obs_vs_pred(fit, dat_long, day_break = 1, title_prefix = nm)
  ggsave(sprintf("figures/traj/SCr_obs_vs_pred_%s.png", nm), p_obs_pred, width = 10, height = 6, dpi = 300)
  
  p_cif <- plot_cif(fit, admin = DAY_MAX_USE, title_prefix = nm)
  ggsave(sprintf("figures/cif/CIF_%s.png", nm), p_cif, width = 10, height = 5, dpi = 300)
}

message("✓ Evaluation terminée (robuste TIMEVAR=", TIMEVAR, ", DAY_MAX=", DAY_MAX_USE, ").")







# ============================================================
# Number/proportion of patients contributing observed data
# across follow-up, by trajectory class
# ============================================================

observability_by_class <- function(fit, dt0, dat_long, days = 0:14) {
  
  pp <- pprob_std(fit) %>%
    dplyr::select(id, class)
  
  # One row per patient for event/censoring time
  patient_level <- dt0 %>%
    dplyr::select(id, tte_pos, event2) %>%
    dplyr::distinct(id, .keep_all = TRUE) %>%
    dplyr::inner_join(pp, by = "id")
  
  # Total class sizes
  class_n <- patient_level %>%
    dplyr::count(class, name = "class_total")
  
  # ----------------------------------------------------------
  # 1. Patients still under ICU observation at each day
  # ----------------------------------------------------------
  at_risk <- tidyr::crossing(
    class = sort(unique(patient_level$class)),
    day = days
  ) %>%
    dplyr::left_join(class_n, by = "class") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      n_under_observation =
        sum(patient_level$class == class &
              patient_level$tte_pos >= day,
            na.rm = TRUE),
      
      pct_under_observation =
        100 * n_under_observation / class_total,
      
      pct_event_truncated =
        100 - pct_under_observation
    ) %>%
    dplyr::ungroup()
  
  # ----------------------------------------------------------
  # 2. Patients with >=1 actual SCr measurement during that day
  # ----------------------------------------------------------
  scr_contributors <- dat_long %>%
    dplyr::inner_join(pp, by = "id") %>%
    dplyr::mutate(
      time_day = to_day(time, dt0),
      day = floor(time_day)
    ) %>%
    dplyr::filter(day %in% days) %>%
    dplyr::distinct(id, class, day) %>%
    dplyr::count(class, day, name = "n_with_scr") %>%
    dplyr::right_join(
      tidyr::crossing(
        class = sort(unique(pp$class)),
        day = days
      ),
      by = c("class", "day")
    ) %>%
    dplyr::mutate(n_with_scr = tidyr::replace_na(n_with_scr, 0L)) %>%
    dplyr::left_join(class_n, by = "class") %>%
    dplyr::mutate(
      pct_with_scr = 100 * n_with_scr / class_total
    )
  
  out <- at_risk %>%
    dplyr::left_join(
      scr_contributors %>%
        dplyr::select(class, day, n_with_scr, pct_with_scr),
      by = c("class", "day")
    )
  
  return(out)
}

fit3 <- fits[[which(vapply(fits, function(x) x$ng, numeric(1)) == 3)[1]]]

obs_table <- observability_by_class(
  fit = fit3,
  dt0 = dt0,
  dat_long = dat_long,
  days = 0:14
)

print(obs_table)

write_csv(
  obs_table,
  "tables/trajectory_observability_by_class.csv"
)
