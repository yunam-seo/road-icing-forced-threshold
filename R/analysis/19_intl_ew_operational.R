# Evaluates the operational nowcast: AUC at several horizons, reliability and Brier skill, persistence triggers (k = 1, 2) and cost-loss value.
# Input: 02_MART/intl_ew_features.rds.
# Output: 05_REPORT/FIGURES/intl_ew_operational.png.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages({library(data.table);library(lightgbm)})
set.seed(2026)
MART<-file.path(ROOT, "02_MART"); FIG<-file.path(ROOT, "05_REPORT/FIGURES")
feats<-c("cold_dose","dep_dose","ft","ts","rh","dts_1h")
Hs<-c(60L,120L,180L,360L,720L); step<-30L; ppd<-24*60/step
D<-readRDS(file.path(MART,"intl_ew_features.rds"))
D<-D[complete.cases(D[,..feats])]; setorder(D,stn,tm_min)
stns<-unique(D$stn); fold<-setNames(sample(rep(1:5,length.out=length(stns))),stns)
D[,fold:=fold[as.character(stn)]]
auc<-function(s,p){o<-is.finite(s);s<-s[o];p<-p[o];np<-as.numeric(sum(p));nn<-as.numeric(sum(!p))
  if(!np||!nn)return(NA);r<-rank(s);(sum(r[p])-np*(np+1)/2)/(np*nn)}
params<-list(objective="binary",learning_rate=0.05,num_leaves=31,min_data_in_leaf=100,
  feature_fraction=0.8,bagging_fraction=0.8,bagging_freq=1,verbosity=-1)

oof_prob<-function(H){
  y<-as.integer(!is.na(D$tno)&D$tno<=H); pr<-numeric(nrow(D))
  for(kf in 1:5){tr<-which(D$fold!=kf);te<-which(D$fold==kf)
    m<-lgb.train(params,lgb.Dataset(as.matrix(D[tr,..feats]),label=y[tr]),nrounds=300)
    p<-predict(m,as.matrix(D[te,..feats]));if(is.matrix(p))p<-p[,1];pr[te]<-p}
  list(y=y,pr=pr)
}

cat("=== (A) Station-held-out CV AUC at multiple horizons ===\n")
store<-list()
for(H in Hs){o<-oof_prob(H);store[[as.character(H)]]<-o
  cat(sprintf("  H=%4dmin(%.0fh): positive rate %.4f  AUC %.3f\n",H,H/60,mean(o$y),auc(o$pr,o$y==1L)))}

brier_skill<-function(y,p){bs<-mean((p-y)^2);bcl<-mean((mean(y)-y)^2);1-bs/bcl}
trigger_eval<-function(H,k,thr){
  y<-store[[as.character(H)]]$y; p<-store[[as.character(H)]]$pr
  dt<-data.table(stn=D$stn,tm=D$tm_min,tno=D$tno,p=p,y=y); setorder(dt,stn,tm)
  dt[,fire:={a<-p>=thr
    if(k==1)a else {ok<-a; for(j in 1:(k-1)){prev<-shift(a,j)&(tm-shift(tm,j)==j*step); ok<-ok&(prev %in% TRUE)}; ok}},by=stn]
  onset_pts<-dt[!is.na(tno)&tno>0&tno<=H]
  onset_pts[,oid:=paste0(stn,"_",tm+tno)]
  lead<-onset_pts[fire==TRUE,.(lead=max(tno)),by=oid]
  allon<-unique(onset_pts$oid)
  pod<-length(unique(lead$oid))/length(allon)
  far_day<-mean(dt[ (is.na(tno)|tno>H) ]$fire)*ppd
  list(pod=pod,far_day=far_day,lead=lead$lead,n_onset=length(allon))
}

for(H in c(180L,360L)){
  o<-store[[as.character(H)]]
  cat(sprintf("\n=== (B/C/D) Operational H=%dh ===  AUC %.3f  Brier skill %.3f\n",
      H/60,auc(o$pr,o$y==1L),brier_skill(o$y,o$pr)))
  qb<-unique(quantile(o$pr,0:10/10)); bin<-cut(o$pr,qb,include.lowest=TRUE,labels=FALSE)
  rel<-data.table(bin,y=o$y,p=o$pr)[,.(predicted=round(mean(p),3),observed=round(mean(y),3),n=.N),by=bin][order(bin)]
  cat("  Reliability (forecast probability vs observed frequency) low/mid/high:\n");print(rel[c(1,round(.N/2),.N)])
  if(H==180L) rel180<-copy(rel)
  thrs<-quantile(o$pr,seq(.90,.999,.005))
  sw<-rbindlist(lapply(thrs,function(t){r<-trigger_eval(H,1L,t)
    data.table(thr=t,pod=r$pod,far=r$far_day,lead=median(r$lead))}))
  for(g in c(0.8,0.7)){op<-sw[which.min(abs(pod-g))]
    r<-trigger_eval(H,1L,op$thr)
    cat(sprintf("  [k=1] POD≈%.2f: achieved %.2f, false alarms %.2f/day, first-trigger lead median %.0f min (IQR %.0f-%.0f)\n",
        g,op$pod,op$far,median(r$lead),quantile(r$lead,.25),quantile(r$lead,.75)))}
  op<-sw[which.min(abs(pod-0.7))]; r2<-trigger_eval(H,2L,op$thr)
  cat(sprintf("  [k=2 persistence] same threshold: POD %.2f, false alarms %.2f/day, median lead %.0f min\n",
      r2$pod,r2$far_day,median(r2$lead)))
  if(H==180L) sw180<-sw
}

o<-store[["180"]]; y<-o$y; p<-o$pr; s<-mean(y)
thr<-quantile(p,seq(.5,.999,.002)); Hr<-Fr<-numeric(length(thr))
for(i in seq_along(thr)){w<-p>=thr[i];Hr[i]<-sum(w&y==1)/sum(y==1);Fr[i]<-sum(w&y==0)/sum(y==0)}
alphas<-seq(0.02,0.98,0.02)
Vmax<-sapply(alphas,function(a){
  v<-(pmin(a,s)-(a*(Hr*s+Fr*(1-s))+s*(1-Hr)))/(pmin(a,s)-s*a); max(v,na.rm=TRUE)})
cat(sprintf("\n(D) cost-loss maximum relative economic value %.3f (near cost/loss ratio α=%.2f)\n",
    max(Vmax),alphas[which.max(Vmax)]))

png(file.path(FIG,"intl_ew_operational.png"),width=2400,height=720,res=200)
pp<-par(mfrow=c(1,4),mar=c(4.4,4.6,2.6,1),cex.lab=1.2,cex.axis=1.05,cex.main=1.15)
av<-sapply(Hs,function(H)auc(store[[as.character(H)]]$pr,store[[as.character(H)]]$y==1L))
plot(seq_along(Hs),av,type="b",pch=19,col="navy",lwd=2,xaxt="n",
  xlab="Forecast lead time",ylab="Station-held-out CV AUC",
  main="(a) Skill by lead time",ylim=c(0.8,0.95))
axis(1,seq_along(Hs),paste0("T+",Hs/60))
plot(rel180$predicted,rel180$observed,type="b",pch=19,col="purple",lwd=2,xlab="Forecast probability",
  ylab="Observed frequency",main="(b) Reliability (T+3)",xlim=c(0,max(rel180$predicted)),ylim=c(0,max(rel180$predicted)))
abline(0,1,lty=2,col="grey50")
plot(sw180$far,sw180$pod,type="l",lwd=2,col="darkred",xlab="False alarms (per day per station)",
  ylab="Event detection (POD)",main="(c) Trigger POD vs false alarms (T+3)",ylim=c(0,1),xlim=c(0,10))
plot(alphas,Vmax,type="l",lwd=2,col="darkgreen",xlab="Cost/loss ratio \u03b1",
  ylab="Relative economic value",main="(d) Cost\u2013loss value (T+3)",ylim=c(0,0.5),xlim=c(0,0.5))
par(pp);dev.off()
cat("Saved: 05_REPORT/FIGURES/intl_ew_operational.png\n")
