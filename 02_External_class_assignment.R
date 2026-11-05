#!/usr/bin/env Rscript
# Dispatcher for the original AUMC evaluation and two external-assignment
# scripts kept under internal/. No patient-level data are included.
`%||%` <- function(x, y) if (is.null(x)) y else x
args <- commandArgs(trailingOnly = TRUE)
opts <- list()
for (arg in args) {
  if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) {
    stop("Expected --name=value argument, received: ", arg)
  }
  key <- sub("^--([^=]+)=.*$", "\\1", arg)
  opts[[key]] <- sub("^[^=]+=", "", arg)
}
cohort <- toupper(opts$cohort %||% "")
if (!cohort %in% c("AUMC", "MIMIC", "GENEVA")) {
  stop("Use --cohort=AUMC, --cohort=MIMIC or --cohort=GENEVA")
}
required <- c("aumc", "models")
if (cohort != "AUMC") required <- c(required, "external")
missing <- required[!vapply(required, function(x) {
  !is.null(opts[[x]]) && nzchar(opts[[x]])
}, logical(1))]
if (length(missing)) stop("Missing arguments: ", paste(missing, collapse=", "))
if (!dir.exists(opts$models)) stop("--models must be an existing directory")
for (key in setdiff(required, "models")) {
  if (!file.exists(opts[[key]])) stop("Not found: ", opts[[key]])
}
out_dir <- opts$out %||% file.path(opts$models, paste0("evaluation_",tolower(cohort)))
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
Sys.setenv(AKI_AUMC_RDS=normalizePath(opts$aumc),
           AKI_EXTERNAL_RDS=if (cohort == "AUMC") "" else normalizePath(opts$external),
           AKI_MODELS_DIR=normalizePath(opts$models),
           AKI_OUT_DIR=normalizePath(out_dir),
           AKI_MODEL_PATTERN=opts$pattern %||% "")
file_arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
if (!length(file_arg)) stop("Run this script using Rscript")
base_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
script <- switch(cohort,
                 AUMC="02_AUMC_model_evaluation.R",
                 MIMIC="02_MIMIC_external_assignment.R",
                 GENEVA="02_Genevadb_external_assignment.R")
cat("[EVAL] cohort:",cohort,"\n")
source(file.path(base_dir,"internal",script), chdir=FALSE)
