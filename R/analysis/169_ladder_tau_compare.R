# Compares dose-only and full nowcast models when the dose decay time tau_d is set to 20, 60 or 1,440 min.
# Inputs: 02_MART/intl_ew_features.rds, intl_event_table.rds; 1-min RWIS files.
# Output: 02_MART/ladder_tau_compare.csv.
# Reads the full 1-min archive (slow).

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages({library(data.table);library(lightgbm)})
MART<-file.path(ROOT, "02_MART")
E<-new.env(); sys.source(code_file("12_intl_earlywarning.R"),envir=E)

set.seed(2026)
H<-180L
D<-readRDS(file.path(MART,"intl_ew_features.rds"))
allf<-c("cold_dose","dep_dose","ft","ts","rh","dts_1h")
D<-D[complete.cases(D[,..allf])]; setorder(D,stn,tm_min)
D[,y:=as.integer(!is.na(tno)&tno<=H)]
stns<-unique(D$stn); fold<-setNames(sample(rep(1:5,length.out=length(stns))),stns)
D[,fold:=fold[as.character(stn)]]
auc<-function(s,p){o<-is.finite(s);s<-s[o];p<-p[o];np<-as.numeric(sum(p));nn<-as.numeric(sum(!p))
  if(!np||!nn)return(NA);r<-rank(s);(sum(r[p])-np*(np+1)/2)/(np*nn)}
params<-list(objective="binary",learning_rate=0.05,num_leaves=31,min_data_in_leaf=100,
  feature_fraction=0.8,bagging_fraction=0.8,bagging_freq=1,verbosity=-1)
cv_auc<-function(X,fs){pr<-numeric(nrow(X))
  for(kf in 1:5){tr<-which(X$fold!=kf);te<-which(X$fold==kf)
    m<-lgb.train(params,lgb.Dataset(as.matrix(X[tr,..fs]),label=X$y[tr]),nrounds=300)
    p<-predict(m,as.matrix(X[te,..fs]));if(is.matrix(p))p<-p[,1];pr[te]<-p}
  auc(pr,X$y==1L)}

base_dyn<-cv_auc(D,c("ts","rh","dts_1h"))
cat(sprintf("Current+change {Ts,RH,ΔTs1h} AUC %.4f\n",base_dyn))
OUT<-list(data.table(tau_min=NA_real_,configuration="plus_dynamics",AUC=base_dyn,n=nrow(D)))
for(tau in c(1440,20,60)){
  t0<-Sys.time(); E$cfg$tau<-tau
  B<-E$build()[,.(stn,tm_min,cold_dose,dep_dose,ft)]
  X<-merge(D[,!c("cold_dose","dep_dose","ft")],B,by=c("stn","tm_min"),all.x=TRUE)
  setorder(X,stn,tm_min)
  stopifnot(nrow(X)==nrow(D), all(complete.cases(X[,..allf])))
  if(tau==1440) stopifnot(isTRUE(all.equal(X$cold_dose,D$cold_dose)), isTRUE(all.equal(X$dep_dose,D$dep_dose)), isTRUE(all.equal(X$ft,D$ft)))
  a_d<-cv_auc(X,c("cold_dose","dep_dose","ft")); a_f<-cv_auc(X,allf)
  cat(sprintf("τd %5.0f min: doses only %.4f · full %.4f · full−current/change %+.4f  (%.0f min)\n",tau,a_d,a_f,a_f-base_dyn,
              as.numeric(difftime(Sys.time(),t0,units="mins"))))
  OUT[[length(OUT)+1]]<-data.table(tau_min=tau,configuration=c("dose_only","full"),AUC=c(a_d,a_f),n=nrow(X))
}
R<-rbindlist(OUT); print(R)
fwrite(R,file.path(MART,"ladder_tau_compare.csv"),bom=TRUE)
cat("Saved: 02_MART/ladder_tau_compare.csv\n")
