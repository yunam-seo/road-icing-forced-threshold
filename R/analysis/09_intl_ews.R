# Computes critical-slowing-down indicators (rolling variance, lag-1 autocorrelation) of detrended surface temperature before onsets and at control times, summarised by Kendall tau.
# Inputs: 02_MART/intl_event_table.rds, intl_cluster_assign.rds; 1-min RWIS files.
# Output: 02_MART/intl_ews.rds.
# Reads the full 1-min archive (slow).

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
source(code_file("00_analysis_window.R"))
set.seed(2026)
MART<-file.path(ROOT, "02_MART")
FIG <-file.path(ROOT, "05_REPORT/FIGURES")
DDIR<-rwis_dir()
cfg<-list(W=720L, rw=360L, step=30L, bw=60, cov_thr=0.8,
          freeze_mon=c(11L,12L,1L,2L,3L), ctrl_per_event=2L, ctrl_buffer=1440L,
          ts_range=c(-40,40))

map_sfs_class<-function(s){o<-rep(NA_character_,length(s));o[s==1L]<-"DRY"
  o[s%in%c(2L,11:14)]<-"WET";o[s==3L]<-"SNOW";o[s%in%c(4L,15L)]<-"ICE";o}

read_ts<-function(f){
  d<-fread(f,select=c("datetime","TR1","SFS","mon"),na.strings=c("-999","-999.0",""))
  d<-d[mon%in%cfg$freeze_mon]; d<-filter_winter(d)
  if(!nrow(d))return(NULL)
  d[(TR1<cfg$ts_range[1])|(TR1>cfg$ts_range[2]),TR1:=NA_real_]
  d[,tm_min:=as.integer(as.POSIXct(strptime(sprintf("%012.0f",as.numeric(datetime)),
            "%Y%m%d%H%M",tz="UTC")))%/%60L]
  d<-d[!is.na(tm_min)]; setorder(d,tm_min)
  d[,cls:=map_sfs_class(SFS)]; d[,.(tm_min,ts=TR1,cls)]
}

gauss_trend<-function(v,bw){
  n<-length(v); idx<-seq_len(n); ok<-which(!is.na(v))
  if(length(ok)<3) return(rep(NA_real_,n))
  sapply(idx,function(i){w<-exp(-0.5*((ok-i)/bw)^2);sum(w*v[ok])/sum(w)})
}
ews_window<-function(v){
  if(mean(!is.na(v))<cfg$cov_thr) return(c(var_tau=NA,ar1_tau=NA,var_end=NA,ar1_end=NA))
  tr<-gauss_trend(v,cfg$bw); r<-v-tr
  st<-seq(1L,length(v)-cfg$rw,by=cfg$step)
  vv<-aa<-numeric(length(st))
  for(k in seq_along(st)){
    seg<-r[st[k]:(st[k]+cfg$rw-1L)]
    vv[k]<-var(seg,na.rm=TRUE)
    a<-seg[-length(seg)]; b<-seg[-1]
    good<-is.finite(a)&is.finite(b)
    aa[k]<-if(sum(good)>30) suppressWarnings(cor(a[good],b[good])) else NA
  }
  kt<-function(y){o<-is.finite(y);if(sum(o)<4)return(NA)
    suppressWarnings(cor(seq_along(y)[o],y[o],method="kendall"))}
  c(var_tau=kt(vv), ar1_tau=kt(aa),
    var_end=tail(vv[is.finite(vv)],1)[1], ar1_end=tail(aa[is.finite(aa)],1)[1])
}

build<-function(nfile=NULL){
  ev<-readRDS(file.path(MART,"intl_event_table.rds")); setkey(ev,stn)
  files<-list.files(DDIR,pattern="\\.csv$",full.names=TRUE)
  if(!is.null(nfile))files<-head(files,nfile)
  sof<-as.integer(sub(".*rwis_(\\d+)_.*","\\1",basename(files)))
  out<-list();k<-0L
  for(i in seq_along(files)){
    sid<-sof[i];evs<-ev[stn==sid];if(!nrow(evs))next
    d<-read_ts(files[i]);if(is.null(d)||nrow(d)<cfg$W)next
    setkey(d,tm_min); tmin<-d$tm_min[1]
    grab<-function(t_end){
      g<-data.table(tm_min=(t_end-cfg$W):t_end,key="tm_min")
      d[g,on="tm_min"]$ts
    }
    for(j in seq_len(nrow(evs))){
      e<-ews_window(grab(evs$onset_min[j]-1L))
      k<-k+1L;out[[k]]<-data.table(type="onset",stn=sid,
        eid=paste0(sid,"_",evs$onset_min[j]),t_end=evs$onset_min[j]-1L,t(e))
    }
    oo<-sort(unique(evs$onset_min))
    cand<-d[cls%in%c("DRY","WET","SNOW")&tm_min>tmin+cfg$W]
    if(nrow(cand)){
      jj<-findInterval(cand$tm_min,oo)
      dl<-ifelse(jj>=1,abs(cand$tm_min-oo[pmax(jj,1L)]),Inf)
      dh<-ifelse(jj<length(oo),abs(cand$tm_min-oo[pmin(jj+1L,length(oo))]),Inf)
      far<-cand[pmin(dl,dh)>=cfg$ctrl_buffer]
      nc<-min(nrow(far),cfg$ctrl_per_event*nrow(evs))
      if(nc>0)for(t_end in far[sample(.N,nc)]$tm_min){
        e<-ews_window(grab(t_end))
        k<-k+1L;out[[k]]<-data.table(type="control",stn=sid,eid=NA_character_,
          t_end=t_end,t(e))
      }
    }
    if(i%%20L==0L)message("  ",i,"/",length(files))
  }
  rbindlist(out,fill=TRUE)
}

if(sys.nframe()==0L){
  args<-commandArgs(trailingOnly=TRUE);test<-length(args)>=1L&&args[1]=="test"
  t0<-Sys.time();res<-build(if(test)5L else NULL)
  cat("Elapsed:",round(difftime(Sys.time(),t0,units="secs"),1),"s\n")
  asg<-readRDS(file.path(MART,"intl_cluster_assign.rds"))
  res<-merge(res,asg[,.(eid,cluster)],by="eid",all.x=TRUE)
  cat("Sample onset",sum(res$type=="onset")," control",sum(res$type=="control"),
      " (onsets with valid EWS",sum(res$type=="onset"&is.finite(res$var_tau)),")\n")
  auc<-function(s,p){o<-is.finite(s);s<-s[o];p<-p[o];r<-rank(s)
    (sum(r[p])-sum(p)*(sum(p)+1)/2)/(sum(p)*sum(!p))}
  pp<-res$type=="onset"
  cat("\n=== Median Kendall τ and discrimination AUC (onset vs control) ===\n")
  for(v in c("var_tau","ar1_tau")){
    cat(sprintf("  %s: onset %.3f / control %.3f  AUC %.3f\n",v,
      median(res[[v]][pp],na.rm=T),median(res[[v]][!pp],na.rm=T),auc(res[[v]],pp)))
  }
  cat("\n=== Onset Kendall τ by pathway (positive = CSD) ===\n")
  print(res[type=="onset"&!is.na(cluster),.(n=.N,
    var_tau=round(median(var_tau,na.rm=T),3),ar1_tau=round(median(ar1_tau,na.rm=T),3),
    pos_var=round(mean(var_tau>0,na.rm=T),2),pos_ar1=round(mean(ar1_tau>0,na.rm=T),2)),
    by=cluster][order(cluster)])
  if(!test){saveRDS(res,file.path(MART,"intl_ews.rds"));cat("\nSaved: intl_ews.rds\n")}
}
