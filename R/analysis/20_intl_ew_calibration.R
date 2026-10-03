# Compares probability calibration options (positive-class weighting, unweighted, unweighted with isotonic calibration) by AUC, Brier skill and reliability.
# Input: 02_MART/intl_ew_features.rds.
# Output: console summary.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages({library(data.table);library(lightgbm)})
set.seed(2026)
MART<-file.path(ROOT, "02_MART")
feats<-c("cold_dose","dep_dose","ft","ts","rh","dts_1h"); H<-180L
D<-readRDS(file.path(MART,"intl_ew_features.rds")); D<-D[complete.cases(D[,..feats])]
D[,y:=as.integer(!is.na(tno)&tno<=H)]
stns<-unique(D$stn); fold<-setNames(sample(rep(1:5,length.out=length(stns))),stns)
D[,fold:=fold[as.character(stn)]]
auc<-function(s,p){o<-is.finite(s);s<-s[o];p<-p[o];np<-as.numeric(sum(p));nn<-as.numeric(sum(!p))
  if(!np||!nn)return(NA);r<-rank(s);(sum(r[p])-np*(np+1)/2)/(np*nn)}
bss<-function(y,p){1-mean((p-y)^2)/mean((mean(y)-y)^2)}
base<-list(objective="binary",learning_rate=0.05,num_leaves=31,min_data_in_leaf=100,
  feature_fraction=0.8,bagging_fraction=0.8,bagging_freq=1,verbosity=-1)

oof<-function(use_spw,calib){
  pr<-numeric(nrow(D))
  for(kf in 1:5){tr<-which(D$fold!=kf);te<-which(D$fold==kf)
    pars<-base; if(use_spw)pars<-c(pars,scale_pos_weight=round((1-mean(D$y[tr]))/mean(D$y[tr])))
    m<-lgb.train(pars,lgb.Dataset(as.matrix(D[tr,..feats]),label=D$y[tr]),nrounds=300)
    p<-predict(m,as.matrix(D[te,..feats]));if(is.matrix(p))p<-p[,1]
    if(calib){
      set.seed(kf); sub<-sample(tr,length(tr)*0.8); cal<-setdiff(tr,sub)
      mc<-lgb.train(pars,lgb.Dataset(as.matrix(D[sub,..feats]),label=D$y[sub]),nrounds=300)
      pc<-predict(mc,as.matrix(D[cal,..feats]));if(is.matrix(pc))pc<-pc[,1]
      ir<-isoreg(pc,D$y[cal])
      p<-approx(sort(pc),ir$yf,xout=p,rule=2)$y }
    pr[te]<-p}
  pr
}
report<-function(pr,lab){
  qb<-quantile(pr,0:10/10,na.rm=TRUE)
  bin<-if(anyDuplicated(qb)){
    as.integer(cut(rank(pr,ties.method="first"),
                   breaks=quantile(seq_along(pr),0:10/10),
                   include.lowest=TRUE,labels=FALSE))
  } else cut(pr,qb,include.lowest=TRUE,labels=FALSE)
  rel<-data.table(bin,y=D$y,p=pr)[,.(predicted=round(mean(p),3),observed=round(mean(y),3)),by=bin][order(bin)]
  cat(sprintf("\n[%s] AUC %.3f  Brier skill %.3f\n  Reliability (deciles 1/5/8/10): forecast %s / observed %s\n",
    lab,auc(pr,D$y==1L),bss(D$y,pr),
    paste(rel$predicted[c(1,5,8,10)],collapse=","),paste(rel$observed[c(1,5,8,10)],collapse=",")))
}
cat("Positive rate:",round(mean(D$y),4),"\n")
report(oof(TRUE ,FALSE),"i. pos-weight")
report(oof(FALSE,FALSE),"ii. natural (unweighted)")
report(oof(FALSE,TRUE ),"iii. natural + isotonic calibration")
