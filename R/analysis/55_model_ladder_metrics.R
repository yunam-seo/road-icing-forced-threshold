# Computes ROC-AUC, PR-AUC, Brier skill and trigger metrics at equal POD for Models 0-2 of the information ladder.
# Input: 02_MART/intl_ew_features.rds.
# Output: 02_MART/model_ladder_metrics.csv.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages({library(data.table);library(lightgbm)})
set.seed(2026)
MART<-file.path(ROOT, "02_MART")
H<-180L; step<-30L; ppd<-24*60/step
source(code_file("_trigger_rule.R"))
D<-readRDS(file.path(MART,"intl_ew_features.rds"))
allf<-c("cold_dose","dep_dose","ft","ts","rh","dts_1h")
D<-D[complete.cases(D[,..allf])]; setorder(D,stn,tm_min)
D[,y:=as.integer(!is.na(tno)&tno<=H)]
stns<-unique(D$stn); fold<-setNames(sample(rep(1:5,length.out=length(stns))),stns)
D[,fold:=fold[as.character(stn)]]

auc<-function(s,p){o<-is.finite(s);s<-s[o];p<-p[o];np<-as.numeric(sum(p));nn<-as.numeric(sum(!p))
  if(!np||!nn)return(NA);r<-rank(s);(sum(r[p])-np*(np+1)/2)/(np*nn)}
ap<-function(s,y){o<-is.finite(s);s<-s[o];yy<-y[o];ordv<-order(-s)
  yy<-yy[ordv];tp<-cumsum(yy);prec<-tp/seq_along(yy);sum(prec[yy==1L])/sum(yy)}
bss<-function(pr,y){bs<-mean((pr-y)^2);bc<-mean((mean(y)-y)^2);1-bs/bc}

params<-list(objective="binary",learning_rate=0.05,num_leaves=31,min_data_in_leaf=100,
  feature_fraction=0.8,bagging_fraction=0.8,bagging_freq=1,verbosity=-1)
cv_pred<-function(fs){pr<-numeric(nrow(D))
  for(kf in 1:5){tr<-which(D$fold!=kf);te<-which(D$fold==kf)
    m<-lgb.train(params,lgb.Dataset(as.matrix(D[tr,..fs]),label=D$y[tr]),nrounds=300)
    p<-predict(m,as.matrix(D[te,..fs]));if(is.matrix(p))p<-p[,1];pr[te]<-p}; pr}

TG <- make_trigger(data.table(stn = D$stn, tm = D$tm_min, tno = D$tno), step)
trig_full <- function(score, target = 0.8) {
  th <- TG$thr_for_pod(score, H, 1L, target)
  if (is.na(th)) return(NULL)
  r <- TG$trig(score, H, 1L, th)
  fire <- score >= th
  TP <- sum(D$y == 1L & fire); FP <- sum(D$y == 0L & fire); FN <- sum(D$y == 1L & !fire)
  list(thr = th, pod = r$pod, far_day = r$far, lead = median(r$lead),
       far_pt = FP / (TP + FP), csi = TP / (TP + FN + FP))
}

cat("Evaluation points",nrow(D)," positive rate",round(mean(D$y),4)," H=",H,"min\n\n")
MODELS<-list(
  "Model0_threshold"=list(fs=NULL, prob=FALSE),
  "Model1_current"  =list(fs=c("ts","rh"), prob=TRUE),
  "Model2_recent"   =list(fs=allf, prob=TRUE))
res<-list()
for(nm in names(MODELS)){
  cfg<-MODELS[[nm]]
  pr<-if(is.null(cfg$fs)) -D$ts else cv_pred(cfg$fs)
  b<-trig_full(pr,0.8)
  row<-data.table(model=nm,
    ROC_AUC=round(auc(pr,D$y==1L),3), PR_AUC=round(ap(pr,D$y),3),
    BSS=if(cfg$prob) round(bss(pr,D$y),3) else NA_real_,
    POD=round(b$pod,2), FA_per_day=round(b$far_day,2),
    FAR_point=round(b$far_pt,3), CSI_point=round(b$csi,3),
    lead_med_min=round(b$lead))
  res[[nm]]<-row; print(row)
}
out<-rbindlist(res)
fwrite(out,file.path(MART,"model_ladder_metrics.csv"),bom=TRUE)
cat("\nSaved: 02_MART/model_ladder_metrics.csv\n")
