#!/usr/bin/env Rscript
# Landmark cause-specific Cox models for subsequent ICU death.
# Baseline clinical model versus clinical model plus frozen M2 probabilities.
# Paired patient-level bootstrap: refit both Cox models in each replicate.
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

argv <- commandArgs(trailingOnly=TRUE)
`%||%` <- function(x,y) if(is.null(x)) y else x
opt <- list()
for (a in argv) {
  if (!startsWith(a,"--") || !grepl("=",a,fixed=TRUE)) stop("Use --key=value: ",a)
  p <- strsplit(sub("^--","",a),"=",fixed=TRUE)[[1]]
  opt[[p[[1]]]] <- paste(p[-1],collapse="=")
}
root <- Sys.getenv("AKI_BASE_DIR")
default_path <- function(...) if (nzchar(root)) file.path(root,...) else ""
pred_path <- opt$pred %||% default_path("XGB_external_AUMC_train__MIMIC_G3_LANDMARK_RISKSET_v4","LANDMARK_pred_dt.csv")
clin_path <- opt$clinical %||% default_path("Prognostic_comparison_mortality_saved_outputs","comparison_variables_MIMIC.rds")
base_path <- opt$baseline %||% default_path("DF","MIMIC","mimic31_aki_admission_base_v3.rds")
scr_path <- opt$scr %||% default_path("DF","MIMIC","timeseries_mimic.rds")
aumc_clin_path <- opt$aumc_clinical %||% default_path("Prognostic_comparison_mortality_saved_outputs","comparison_variables_AUMC.rds")
aumc_base_path <- opt$aumc_baseline %||% default_path("DF","AUMC","aumc_bf.rds")
out_dir <- opt$out %||% default_path("Landmark_Cox_bootstrap")
n_boot <- as.integer(opt$boot %||% "1000")
seed <- as.integer(opt$seed %||% "20260921")
if (!is.finite(n_boot) || n_boot < 1) stop("--boot must be a positive integer")
if (!nzchar(out_dir)) stop("Provide --out or AKI_BASE_DIR")
for (p in c(pred_path,clin_path,base_path,scr_path,aumc_clin_path,aumc_base_path)) {
  if (!file.exists(p)) stop("Missing: ",p)
}
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

to_posix <- function(x) {
  if (inherits(x,"POSIXt")) return(as.POSIXct(x,tz="UTC"))
  if (is.numeric(x)) return(as.POSIXct(x,origin="1970-01-01",tz="UTC"))
  as.POSIXct(x,tz="UTC")
}
pred <- if (grepl("\\.rds$",pred_path,ignore.case=TRUE)) as.data.table(readRDS(pred_path)) else fread(pred_path)
required_pred <- c("id","model","day","pC1","pC2","pC3")
if (length(setdiff(required_pred,names(pred)))) stop("Predictions missing: ",paste(setdiff(required_pred,names(pred)),collapse=", "))
pred <- pred[model=="M2" & day %in% 0:2,
             .(stay_id=as.integer(id), landmark_h=as.integer(24*(day+1)),
               pC1=as.numeric(pC1),pC2=as.numeric(pC2),pC3=as.numeric(pC3))]
if (anyDuplicated(pred[,.(stay_id,landmark_h)])) stop("Duplicate M2 predictions per ICU stay and landmark")
if (any(!is.finite(pred$pC1+pred$pC2+pred$pC3)) ||
    any(abs(pred$pC1+pred$pC2+pred$pC3-1)>1e-6)) stop("Invalid M2 class probability sums")

cmp <- as.data.table(readRDS(clin_path))
clin_vars <- c("stay_id","egfr_baseline_mlmin","sofa_total_adm","vaso_rate_max_0_24","surgery_context")
if (length(setdiff(clin_vars,names(cmp)))) stop("Missing clinical variables: ",paste(setdiff(clin_vars,names(cmp)),collapse=", "))
clin <- unique(cmp[,..clin_vars],by="stay_id")
if (anyDuplicated(cmp$stay_id) && nrow(unique(cmp[,..clin_vars])) != nrow(clin)) stop("Inconsistent covariates within stay")
clin[,stay_id:=as.integer(stay_id)]

bf <- as.data.table(readRDS(base_path))
required_bf <- c("stay_id","intime","outtime","deathtime")
if (length(setdiff(required_bf,names(bf)))) stop("Missing event timestamps")
bf <- unique(bf[,..required_bf],by="stay_id")
bf[,stay_id:=as.integer(stay_id)]
for (v in c("intime","outtime","deathtime")) set(bf,j=v,value=to_posix(bf[[v]]))
bf[, death_icu := as.integer(!is.na(deathtime) & !is.na(intime) &
                              !is.na(outtime) & deathtime>=intime & deathtime<=outtime)]

dat <- merge(pred,clin,by="stay_id",all.x=TRUE)
dat <- merge(dat,bf,by="stay_id",all.x=TRUE)
dat[,surgery_context:=as.factor(surgery_context)]

# Same first- and last-observed SCr definitions as the H48 level/shape analysis.
ts <- as.data.table(readRDS(scr_path))
if (length(setdiff(c("stay_id","t_hours","var","value"),names(ts)))) {
  stop("SCr timeseries requires stay_id, t_hours, var and value")
}
ts[,stay_id:=as.integer(stay_id)]
scr48 <- ts[var=="scr" & is.finite(t_hours) & t_hours>=0 & t_hours<=48,
            .(stay_id,t_hours,scr=as.numeric(value))]
scr48 <- scr48[is.finite(scr)]
setorder(scr48,stay_id,t_hours)
scr_descriptor <- scr48[,.(scr_admission=first(scr),
                             scr_last_48h=last(scr)),by=stay_id]
scr_descriptor[,delta_scr_48h:=scr_last_48h-scr_admission]
build_landmark <- function(h) {
  d <- copy(dat[landmark_h==h])
  d[,landmark_time:=intime+h*3600]
  d <- d[!is.na(intime) & !is.na(outtime) & outtime>landmark_time &
           (is.na(deathtime) | deathtime>landmark_time)]
  d[,event_death_cs:=death_icu]
  d[,end_time:=outtime]
  d[event_death_cs==1L,end_time:=deathtime]
  d[,time_landmark:=as.numeric(difftime(end_time,landmark_time,units="days"))]
  d <- d[is.finite(time_landmark) & time_landmark>0]
  d <- d[complete.cases(d[,c("egfr_baseline_mlmin","sofa_total_adm",
                              "vaso_rate_max_0_24","surgery_context","pC1","pC2"),with=FALSE])]
  d
}
base_formula <- Surv(time_landmark,event_death_cs) ~
  egfr_baseline_mlmin + sofa_total_adm + vaso_rate_max_0_24 + surgery_context
plus_formula <- update(base_formula,. ~ . + pC1 + pC2)
simple_formula <- update(base_formula,. ~ . + scr_admission + delta_scr_48h)

fit_pair <- function(d, with_scr=FALSE) {
  # Identical rows for both models, no additional complete-case filtering.
  a <- coxph(base_formula,data=d,x=TRUE,y=TRUE,ties="efron")
  b <- coxph(plus_formula,data=d,x=TRUE,y=TRUE,ties="efron")
  ca <- as.numeric(concordance(a)$concordance)
  cb <- as.numeric(concordance(b)$concordance)
  res <- c(clinical=ca,clinical_M2=cb,delta=cb-ca)
  if (with_scr) {
    simple_fit <- coxph(simple_formula,data=d,x=TRUE,y=TRUE,ties="efron")
    cs <- as.numeric(concordance(simple_fit)$concordance)
    res <- c(res,clinical_SCr=cs,delta_SCr=cs-ca)
  }
  res
}
rows <- list()
for (h in c(24L,48L,72L)) {
  d <- build_landmark(h)
  if (h==48L) {
    d <- merge(d,scr_descriptor,by="stay_id",all.x=TRUE)
    d <- d[is.finite(scr_admission) & is.finite(delta_scr_48h)]
  }
  cat("H",h,": N=",nrow(d),", ICU deaths=",sum(d$event_death_cs),"\n",sep="")
  if (nrow(d)<50 || sum(d$event_death_cs)<2) stop("Insufficient events at H",h)
  observed <- fit_pair(d,with_scr=(h==48L))
  set.seed(seed+h)
  boots <- matrix(NA_real_,nrow=n_boot,ncol=length(observed),
                  dimnames=list(NULL,names(observed)))
  for (b in seq_len(n_boot)) {
    ix <- sample.int(nrow(d),nrow(d),replace=TRUE)
    boots[b,] <- tryCatch(fit_pair(d[ix],with_scr=(h==48L)),
                         error=function(e) rep(NA_real_,length(observed)))
  }
  valid <- complete.cases(boots)
  if (sum(valid)<ceiling(0.95*n_boot)) stop("Too many failed bootstrap refits at H",h,": ",sum(valid),"/",n_boot)
  ci <- apply(boots[valid,,drop=FALSE],2,quantile,probs=c(.025,.975),names=FALSE)
  rows[[as.character(h)]] <- data.table(
    landmark_h=h,n=nrow(d),icu_deaths=sum(d$event_death_cs),
    clinical=observed[["clinical"]],
    clinical_low=ci[1,"clinical"],clinical_high=ci[2,"clinical"],
    clinical_M2=observed[["clinical_M2"]],
    clinical_M2_low=ci[1,"clinical_M2"],clinical_M2_high=ci[2,"clinical_M2"],
    delta=observed[["delta"]],delta_low=ci[1,"delta"],delta_high=ci[2,"delta"],
    p_boot=min(1,2*min(mean(boots[valid,"delta"]<=0),mean(boots[valid,"delta"]>=0))),
    bootstrap_valid=sum(valid),
    clinical_SCr=if(h==48L) observed[["clinical_SCr"]] else NA_real_,
    clinical_SCr_low=if(h==48L) ci[1,"clinical_SCr"] else NA_real_,
    clinical_SCr_high=if(h==48L) ci[2,"clinical_SCr"] else NA_real_,
    delta_SCr=if(h==48L) observed[["delta_SCr"]] else NA_real_,
    delta_SCr_low=if(h==48L) ci[1,"delta_SCr"] else NA_real_,
    delta_SCr_high=if(h==48L) ci[2,"delta_SCr"] else NA_real_)
  fwrite(as.data.table(boots),file.path(out_dir,sprintf("bootstrap_H%d.csv",h)))
}
results <- rbindlist(rows)
print(results)
fwrite(results,file.path(out_dir,"Table_S19_landmark_clinical_vs_M2.csv"))

# Admission-time retrospective clinical benchmark: fit in AmsterdamUMCdb,
# evaluate the fixed model in MIMIC-IV, using ICU death/discharge only.
# This is intentionally distinct from the landmark Cox refits in MIMIC-IV.
make_admission_cohort <- function(cmp_file,baseline_file) {
  cc <- as.data.table(readRDS(cmp_file))
  cols <- c(clin_vars,"traj_class")
  missing <- setdiff(cols,names(cc))
  if (length(missing)) stop("Admission clinical dataset lacks: ",paste(missing,collapse=", "))
  cc <- unique(cc[,..cols],by="stay_id")
  cc[,stay_id:=as.integer(stay_id)]
  bb <- as.data.table(readRDS(baseline_file))
  if (length(setdiff(required_bf,names(bb)))) stop("Missing ICU timestamps")
  bb <- unique(bb[,..required_bf],by="stay_id")
  bb[,stay_id:=as.integer(stay_id)]
  for (v in c("intime","outtime","deathtime")) set(bb,j=v,value=to_posix(bb[[v]]))
  bb[,death_icu:=as.integer(!is.na(intime) & !is.na(outtime) &
                             !is.na(deathtime) & deathtime>=intime & deathtime<=outtime)]
  bb[,end_time:=outtime]
  bb[death_icu==1L,end_time:=deathtime]
  bb[,time_to_event_days:=as.numeric(difftime(end_time,intime,units="days"))]
  d <- merge(cc,bb[,.(stay_id,time_to_event_days,death_icu)],by="stay_id",all.x=TRUE)
  d <- d[is.finite(time_to_event_days) & time_to_event_days>0 &
           !is.na(death_icu) & !is.na(traj_class)]
  d <- d[complete.cases(d[,c("egfr_baseline_mlmin","sofa_total_adm",
                              "vaso_rate_max_0_24","surgery_context"),with=FALSE])]
  d[,surgery_context:=factor(surgery_context,levels=c(0,1),
                            labels=c("Non-surgical","Surgical"))]
  if (anyNA(d$surgery_context)) stop("Unexpected surgery_context coding")
  d
}
adm_A <- make_admission_cohort(aumc_clin_path,aumc_base_path)
adm_M <- make_admission_cohort(clin_path,base_path)
admission_formula <- Surv(time_to_event_days,death_icu) ~
  egfr_baseline_mlmin + sofa_total_adm + vaso_rate_max_0_24 + surgery_context
fit_admission <- coxph(admission_formula,data=adm_A,x=TRUE,y=TRUE,ties="efron")
lp_M <- as.numeric(predict(fit_admission,newdata=adm_M,type="lp"))
c_admission <- as.numeric(concordance(Surv(adm_M$time_to_event_days,adm_M$death_icu) ~ lp_M,
                                      reverse=TRUE)$concordance)
set.seed(seed+900L)
admission_boot <- replicate(n_boot,{
  ix <- sample.int(nrow(adm_M),nrow(adm_M),replace=TRUE)
  tryCatch(as.numeric(concordance(
    Surv(adm_M$time_to_event_days[ix],adm_M$death_icu[ix]) ~ lp_M[ix],
    reverse=TRUE)$concordance),error=function(e) NA_real_)
})
admission_boot <- admission_boot[is.finite(admission_boot)]
if (length(admission_boot)<ceiling(.95*n_boot)) stop("Too many invalid retrospective bootstrap samples")
retrospective_clinical <- data.table(
  population="MIMIC-IV external validation",
  outcome="ICU mortality",
  n=nrow(adm_M),icu_deaths=sum(adm_M$death_icu),
  cindex=c_admission,
  cindex_low=unname(quantile(admission_boot,.025)),
  cindex_high=unname(quantile(admission_boot,.975)),
  training_n=nrow(adm_A),
  training_icu_deaths=sum(adm_A$death_icu))
print(retrospective_clinical)
fwrite(retrospective_clinical,file.path(out_dir,"Table_S19_retrospective_clinical_reference.csv"))
writeLines(capture.output(sessionInfo()),file.path(out_dir,"sessionInfo.txt"))
