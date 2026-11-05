#!/usr/bin/env Rscript
# Landmark enrichment of subsequent ICU mortality in MIMIC-IV.
# M2 top-decile/top-quintile selection, trial sample-size projections,
# and comparisons with KDIGO stage and renal replacement therapy.

argv <- commandArgs(trailingOnly = TRUE)
opt <- list()
for (a in argv) {
  if (!startsWith(a, "--") || !grepl("=", a, fixed = TRUE)) {
    stop("Use --key=value: ", a)
  }
  kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
  opt[[kv[[1]]]] <- paste(kv[-1], collapse = "=")
}

if (!is.null(opt$base_dir)) Sys.setenv(AKI_BASE_DIR = opt$base_dir)
base_dir <- Sys.getenv("AKI_BASE_DIR")
if (!nzchar(base_dir)) stop("Set AKI_BASE_DIR or provide --base_dir=/path/to/project")
if (!dir.exists(base_dir)) stop("Data directory does not exist: ", base_dir)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run with Rscript 06_Enrichment_analysis.R")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))

scripts <- c(
  file.path(script_dir, "internal", "06_enrichment_core.R"),
  file.path(script_dir, "internal", "06_enrichment_comparators.R")
)
for (script in scripts) {
  if (!file.exists(script)) stop("Missing analysis module: ", script)
  cat("\nRunning enrichment module: ", basename(script), "\n", sep = "")
  source(script, local = new.env(parent = globalenv()), echo = FALSE)
}
