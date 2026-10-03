# Builds the nowcast evaluation table (doses, Ts, RH, recent Ts change at each evaluation time) and evaluates onset prediction with station-held-out CV.
# Inputs: 02_MART/intl_event_table.rds; 1-min RWIS files.
# Outputs: 02_MART/intl_ew_features.rds, intl_earlywarning_eval.rds.
# Reads the full 1-min archive (slow).

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages({library(data.table);library(lightgbm)})
set.seed(2026)
MART<-file.path(ROOT, "02_MART")
source(code_file("07_intl_cumstress.R"),local=TRUE)
DDIR<-cfg$data_dir<-rwis_dir()
step_eval<-30L; lead_label<-180L

build<-function(nfile=NULL){
  ev<-readRDS(file.path(MART,"intl_event_table.rds"));setkey(ev,stn)
  files<-list.files(DDIR,pattern="\\.csv$",full.names=TRUE)
  if(!is.null(nfile))files<-head(files,nfile)
  sof<-as.integer(sub(".*rwis_(\\d+)_.*","\\1",basename(files)))
  out<-list();k<-0L
  for(i in seq_along(files)){
    sid<-sof[i];evs<-ev[stn==sid];if(!nrow(evs))next
    d<-read_station(files[i]);if(is.null(d)||nrow(d)<120L)next
    L<-leaky(d,cfg$tau,cfg$gap_max)
    M<-merge(d[,.(tm_min,ts,ta,rh)],L[,.(tm_min,cold_dose,dep_dose,ft,cls)],by="tm_min")
    setkey(M,tm_min)
    M[,ts_1h:=M[data.table(tm_min=M$tm_min-60L),on="tm_min",roll="nearest"]$ts]
    M[,dts_1h:=ts-ts_1h]
    E<-M[cls%in%c("DRY","WET","SNOW") & (tm_min%%step_eval==0L) & is.finite(cold_dose)]
    if(!nrow(E))next
    oo<-sort(unique(evs$onset_min))
    j<-findInterval(E$tm_min,oo)
    nxt<-ifelse(j<length(oo),oo[pmin(j+1L,length(oo))],NA_integer_)
    E[,tno:=nxt-tm_min]
    E[is.na(tno)|tno<0,tno:=NA_integer_]
    k<-k+1L
    out[[k]]<-E[,.(stn=sid,tm_min,cold_dose,dep_dose,ft,ts,rh,dts_1h,tno)]
    if(i%%20L==0L)message("  ",i,"/",length(files))
  }
  rbindlist(out,fill=TRUE)
}

if(sys.nframe()==0L){
  args<-commandArgs(trailingOnly=TRUE);test<-length(args)>=1L&&args[1]=="test"
  t0<-Sys.time();D<-build(if(test)8L else NULL)
  cat("Build:",round(difftime(Sys.time(),t0,units="secs"),1),"s  evaluation points",nrow(D),
      " stations",uniqueN(D$stn),"\n")
  D<-D[complete.cases(D[,.(cold_dose,dep_dose,ft,ts,rh,dts_1h)])]
  D[,y:=as.integer(!is.na(tno)&tno<=lead_label)]
  cat("Positive rate (",lead_label,"min):",round(mean(D$y),4)," positives",sum(D$y),"\n")

  feats<-c("cold_dose","dep_dose","ft","ts","rh","dts_1h")
  stns<-unique(D$stn);fold<-setNames(sample(rep(1:5,length.out=length(stns))),stns)
  D[,fold:=fold[as.character(stn)]]
  auc<-function(s,p){o<-is.finite(s);s<-s[o];p<-p[o];if(!sum(p)||!sum(!p))return(NA)
    r<-rank(s);{n1<-as.numeric(sum(p));n0<-as.numeric(sum(!p));(sum(r[p])-n1*(n1+1)/2)/(n1*n0)}}
  pr<-numeric(nrow(D))
  for(kf in 1:5){
    tr<-D[fold!=kf];te<-which(D$fold==kf)
    m<-lgb.train(list(objective="binary",learning_rate=0.05,num_leaves=31,
      min_data_in_leaf=100,feature_fraction=0.8,bagging_fraction=0.8,bagging_freq=1,
      scale_pos_weight=round((1-mean(tr$y))/mean(tr$y)),verbosity=-1),
      lgb.Dataset(as.matrix(tr[,..feats]),label=tr$y),nrounds=300)
    p<-predict(m,as.matrix(D[te,..feats]));if(is.matrix(p))p<-p[,1];pr[te]<-p
  }
  D[,score:=pr]
  cat(sprintf("\nStation-held-out 5-fold CV (LightGBM, label = onset within %dh):\n  AUC %.3f\n",
      lead_label,auc(D$score,D$y==1L)))
  thr<-quantile(D$score,probs=seq(.5,.999,.001))
  perf<-rbindlist(lapply(thr,function(t){pp<-D$score>=t
    data.table(thr=t,alarm=mean(pp),
      prec=sum(pp&D$y==1)/max(1,sum(pp)),rec=sum(pp&D$y==1)/sum(D$y==1))}))
  op<-perf[which.min(abs(perf$rec-0.7))]
  cat(sprintf("  [point level] operating point (recall≈0.70): precision %.3f, alarm rate %.3f\n",op$prec,op$alarm))
  D[,onset_id:=ifelse(!is.na(tno),paste0(stn,"_",tm_min+tno),NA_character_)]
  evD<-D[y==1L]
  det<-evD[,.(detected=any(score>=op$thr)),by=onset_id]
  ev_recall<-mean(det$detected)
  far_rate<-mean(D[y==0L]$score>=op$thr)
  far_per_day<-far_rate*(24*60/step_eval)
  hit<-evD[score>=op$thr]
  cat(sprintf("  [event level] detection rate %.3f (of %d onsets), false alarm rate %.4f (≈%.1f per station-day)\n",
      ev_recall,uniqueN(det$onset_id),far_rate,far_per_day))
  cat(sprintf("  Lead time of hits (min): median %.0f, IQR %.0f~%.0f\n",
      median(hit$tno),quantile(hit$tno,.25),quantile(hit$tno,.75)))
  cat("\nUnivariate AUC:\n")
  for(v in feats) cat(sprintf("  %s: %.3f\n",v,auc(D[[v]],D$y==1L)))
  if(!test){saveRDS(D[,.(stn,tm_min,score,y,tno)],
    file.path(MART,"intl_earlywarning_eval.rds"));cat("\nSaved: intl_earlywarning_eval.rds\n")
    saveRDS(D[,.(stn,tm_min,cold_dose,dep_dose,ft,ts,rh,dts_1h,tno)],
      file.path(MART,"intl_ew_features.rds"))
    cat("Saved: intl_ew_features.rds · rows",nrow(D),"\n")}
}
