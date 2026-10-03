# Exports the winter-stratified CSD AUC and the AUC by forecast lead with and without surrounding fields as CSV tables.
# Inputs: 02_MART/ews_winter_strat.rds, intl_ew_pattern_store.rds.
# Outputs: 02_MART/ews_winter_strat.csv, pattern_by_lead_auc.csv.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
MART <- file.path(ROOT, "02_MART")

w <- as.data.table(readRDS(file.path(MART, "ews_winter_strat.rds")))
setnames(w, c("var_AUC", "ar1_AUC"), c("auc_variance_trend", "auc_ar1_trend"), skip_absent = TRUE)
print(w)
fwrite(w, file.path(MART, "ews_winter_strat.csv"), bom = TRUE)
cat(sprintf("
Winter-stratified range — variance %.3f–%.3f · AR1 %.3f–%.3f
",
            min(w$auc_variance_trend), max(w$auc_variance_trend),
            min(w$auc_ar1_trend), max(w$auc_ar1_trend)))
cat("Saved: 02_MART/ews_winter_strat.csv
")

auc41<-function(s,p){o<-is.finite(s);s<-s[o];p<-p[o];np<-as.numeric(sum(p));nn<-as.numeric(sum(!p))
  if(!np||!nn)return(NA_real_);r<-rank(s);(sum(r[p])-np*(np+1)/2)/(np*nn)}
ST<-readRDS(file.path(MART,"intl_ew_pattern_store.rds"))
LEAD<-rbindlist(lapply(names(ST),function(h){
  data.table(lead_min=as.integer(h), lead_h=sprintf("T+%d",as.integer(h)%/%60L),
             auc_rwis_only=auc41(ST[[h]]$b$pr,ST[[h]]$b$y==1L),
             auc_with_surrounding_fields=auc41(ST[[h]]$p$pr,ST[[h]]$p$y==1L))}))
LEAD[,difference:=auc_with_surrounding_fields-auc_rwis_only]
print(LEAD)
fwrite(LEAD,file.path(MART,"pattern_by_lead_auc.csv"),bom=TRUE)
cat(sprintf("Gain from surrounding fields, range %.3f–%.3f\n",min(LEAD$difference),max(LEAD$difference)))
cat("Saved: 02_MART/pattern_by_lead_auc.csv\n")
