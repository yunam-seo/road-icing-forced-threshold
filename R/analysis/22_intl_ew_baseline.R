# Compares nowcast models of increasing information (naive Ts, doses only, current state, recent change, full) with station-held-out CV.
# Input: 02_MART/intl_ew_features.rds.
# Outputs: 02_MART/model_ladder_full.csv; 05_REPORT/FIGURES/intl_ew_baseline.png.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages({library(data.table);library(lightgbm)})
set.seed(2026)
MART<-file.path(ROOT, "02_MART"); FIG<-file.path(ROOT, "05_REPORT/FIGURES")
H<-180L; step<-30L; ppd<-24*60/step
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
cv_auc<-function(fs){pr<-numeric(nrow(D))
  for(kf in 1:5){tr<-which(D$fold!=kf);te<-which(D$fold==kf)
    m<-lgb.train(params,lgb.Dataset(as.matrix(D[tr,..fs]),label=D$y[tr]),nrounds=300)
    p<-predict(m,as.matrix(D[te,..fs]));if(is.matrix(p))p<-p[,1];pr[te]<-p}; pr}

trig<-function(score,target=0.8,k=1L){
  dt<-data.table(stn=D$stn,tm=D$tm_min,tno=D$tno,s=score,y=D$y);setorder(dt,stn,tm)
  thrs<-quantile(score,seq(.80,.999,.002))
  best<-NULL
  for(t in thrs){
    dt[,fire:={a<-s>=t; if(k==1)a else {ok<-a;for(j in 1:(k-1)){p<-shift(a,j)&(tm-shift(tm,j)==j*step);ok<-ok&(p%in%TRUE)};ok}},by=stn]
    op<-dt[!is.na(tno)&tno>0&tno<=H];op[,oid:=paste0(stn,"_",tm+tno)]
    pod<-uniqueN(op[fire==TRUE]$oid)/uniqueN(op$oid)
    if(is.null(best)||abs(pod-target)<abs(best$pod-target)){
      far<-mean(dt[is.na(tno)|tno>H]$fire)*ppd
      lead<-op[fire==TRUE,.(l=max(tno)),by=oid]$l
      best<-list(thr=t,pod=pod,far=far,lead=median(lead))}
  }; best
}

cat("Positive rate:",round(mean(D$y),4)," H=",H,"min\n\n=== (1) Model AUC (station-held-out CV) ===\n")
sets<-list(
  "naive Ts (−Ts, no training)"=NULL,
  "instant {Ts,RH}"=c("ts","rh"),
  "+short-term dynamics {…,ΔTs1h}"=c("ts","rh","dts_1h"),
  "doses only {cold,dep,ft}"=c("cold_dose","dep_dose","ft"),
  "full {instant+dynamics+doses}"=allf)
aucs<-list()
for(nm in names(sets)){fs<-sets[[nm]]
  pr<-if(is.null(fs)) -D$ts else cv_auc(fs)
  aucs[[nm]]<-pr; cat(sprintf("  %-26s AUC %.3f\n",nm,auc(pr,D$y==1L)))}
cat(sprintf("\n  Added value: full−instant ΔAUC = %+.3f,  full−naive = %+.3f\n",
  auc(aucs[["full {instant+dynamics+doses}"]],D$y==1L)-auc(aucs[["instant {Ts,RH}"]],D$y==1L),
  auc(aucs[["full {instant+dynamics+doses}"]],D$y==1L)-auc(-D$ts,D$y==1L)))

cat("\n=== (2) Trigger performance: naive Ts vs full model (POD≈0.80) ===\n")
for(nm in c("naive Ts (−Ts, no training)","full {instant+dynamics+doses}")){
  b<-trig(aucs[[nm]],0.8)
  cat(sprintf("  %-26s POD %.2f, false alarms %.2f/day, median lead %.0f min\n",nm,b$pod,b$far,b$lead))}

LAD<-data.table(
  stage=c("naive_ts","dose_only","current","plus_dynamics","full"),
  description=c("surface-temperature threshold only","cumulative antecedent conditions only","current Ts and RH",
         "+ recent 1-h surface-temperature change","full (current + change + cumulative)"),
  AUC=round(sapply(c("naive Ts (−Ts, no training)","doses only {cold,dep,ft}","instant {Ts,RH}",
                     "+short-term dynamics {…,ΔTs1h}","full {instant+dynamics+doses}"),
                   function(nm) auc(aucs[[nm]],D$y==1L)),4))
fwrite(LAD,file.path(MART,"model_ladder_full.csv"),bom=TRUE)
cat("
Saved: 02_MART/model_ladder_full.csv
"); print(LAD)

png(file.path(FIG,"intl_ew_baseline.png"),width=560,height=400,res=110)
par(mar=c(4.4,4.4,2.4,1));cols<-c("grey55","#2980b9","#27ae60","#e67e22","#c0392b")
nb<-names(sets);plot(NA,xlim=c(0.5,5.5),ylim=c(0.6,0.95),xaxt="n",xlab="",ylab="Station-held-out CV AUC",
  main="Added value of cumulative doses (T+3)")
av<-sapply(nb,function(nm)auc(aucs[[nm]],D$y==1L))
points(1:5,av,pch=19,col=cols,cex=1.6);text(1:5,av,sprintf("%.3f",av),pos=3,cex=.8)
axis(1,1:5,c("Naive\nTs","Instant\nTs,RH","+ΔTs","Doses only","full"),padj=.5,cex.axis=.8)
abline(h=av[2],lty=3,col="grey");dev.off()
cat("\nSaved: 05_REPORT/FIGURES/intl_ew_baseline.png\n")
