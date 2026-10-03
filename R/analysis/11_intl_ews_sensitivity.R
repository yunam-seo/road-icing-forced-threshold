# Tests the CSD result across window length W and detrending bandwidth (9 combinations); all stations by default, an optional argument limits the number of stations.
# Inputs: 02_MART/intl_event_table.rds; 1-min RWIS files.
# Output: 02_MART/ews_sensitivity_9combo.csv.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
source(code_file("00_analysis_window.R"))
set.seed(2026)
MART<-file.path(ROOT, "02_MART")
DDIR<-rwis_dir()
ARGS<-commandArgs(trailingOnly=TRUE)
NSTN<-if(length(ARGS)>=1L) as.integer(ARGS[1]) else NA_integer_; Wmax<-1440L
Ws<-c(360L,720L,1440L); bws<-c(30,60,120); cov_thr<-0.8; ctrl_per<-2L; ctrl_buf<-1440L
fmon<-c(11L,12L,1L,2L,3L)
map_cls<-function(s){o<-rep(NA_character_,length(s));o[s==1L]<-"DRY"
  o[s%in%c(2L,11:14)]<-"WET";o[s==3L]<-"SNOW";o[s%in%c(4L,15L)]<-"ICE";o}
read_ts<-function(f){d<-fread(f,select=c("datetime","TR1","SFS","mon"),
    na.strings=c("-999","-999.0",""));d<-d[mon%in%fmon];d<-filter_winter(d)
    if(!nrow(d))return(NULL)
  d[(TR1< -40)|(TR1>40),TR1:=NA_real_]
  d[,tm_min:=as.integer(as.POSIXct(strptime(sprintf("%012.0f",as.numeric(datetime)),
    "%Y%m%d%H%M",tz="UTC")))%/%60L];d<-d[!is.na(tm_min)];setorder(d,tm_min)
  d[,cls:=map_cls(SFS)];setkey(d,tm_min);d[,.(tm_min,ts=TR1,cls)]}
gtrend<-function(v,bw){ok<-which(!is.na(v));if(length(ok)<3)return(rep(NA_real_,length(v)))
  sapply(seq_along(v),function(i){w<-exp(-0.5*((ok-i)/bw)^2);sum(w*v[ok])/sum(w)})}
ews<-function(v,W,bw){
  v<-tail(v,W+1L); if(mean(!is.na(v))<cov_thr)return(c(NA,NA))
  rw<-as.integer(W/2);step<-as.integer(W/24)
  r<-v-gtrend(v,bw); st<-seq(1L,length(v)-rw,by=step); vv<-aa<-numeric(length(st))
  for(k in seq_along(st)){seg<-r[st[k]:(st[k]+rw-1L)];vv[k]<-var(seg,na.rm=TRUE)
    a<-seg[-length(seg)];b<-seg[-1];g<-is.finite(a)&is.finite(b)
    aa[k]<-if(sum(g)>30)suppressWarnings(cor(a[g],b[g]))else NA}
  kt<-function(y){o<-is.finite(y);if(sum(o)<4)return(NA)
    suppressWarnings(cor(seq_along(y)[o],y[o],method="kendall"))}
  c(kt(vv),kt(aa))}

ev<-readRDS(file.path(MART,"intl_event_table.rds"));setkey(ev,stn)
files<-list.files(DDIR,pattern="\\.csv$",full.names=TRUE)
sof<-as.integer(sub(".*rwis_(\\d+)_.*","\\1",basename(files)))
sel<-which(sof%in%ev$stn); if(!is.na(NSTN)) sel<-sel[1:min(NSTN,length(sel))]
combos<-CJ(W=Ws,bw=bws)
acc<-vector("list",nrow(combos))
for(ci in seq_len(nrow(combos))) acc[[ci]]<-list(on=NULL,ct=NULL)

for(i in sel){
  sid<-sof[i];evs<-ev[stn==sid];if(!nrow(evs))next
  d<-read_ts(files[i]);if(is.null(d)||nrow(d)<Wmax)next
  tmin<-d$tm_min[1]
  grab<-function(te){d[data.table(tm_min=(te-Wmax):te,key="tm_min"),on="tm_min"]$ts}
  ons<-lapply(evs$onset_min-1L,grab)
  oo<-sort(unique(evs$onset_min));cand<-d[cls%in%c("DRY","WET","SNOW")&tm_min>tmin+Wmax]
  cts<-list()
  if(nrow(cand)){jj<-findInterval(cand$tm_min,oo)
    dl<-ifelse(jj>=1,abs(cand$tm_min-oo[pmax(jj,1L)]),Inf)
    dh<-ifelse(jj<length(oo),abs(cand$tm_min-oo[pmin(jj+1L,length(oo))]),Inf)
    far<-cand[pmin(dl,dh)>=ctrl_buf];nc<-min(nrow(far),ctrl_per*nrow(evs))
    if(nc>0)cts<-lapply(far[sample(.N,nc)]$tm_min,grab)}
  for(ci in seq_len(nrow(combos))){W<-combos$W[ci];bw<-combos$bw[ci]
    if(length(ons))acc[[ci]]$on<-rbind(acc[[ci]]$on,t(sapply(ons,ews,W=W,bw=bw)))
    if(length(cts))acc[[ci]]$ct<-rbind(acc[[ci]]$ct,t(sapply(cts,ews,W=W,bw=bw)))}
}
auc<-function(on,ct){s<-c(on,ct);p<-c(rep(TRUE,length(on)),rep(FALSE,length(ct)))
  o<-is.finite(s);s<-s[o];p<-p[o];if(!sum(p)||!sum(!p))return(NA);r<-rank(s)
  (sum(r[p])-sum(p)*(sum(p)+1)/2)/(sum(p)*sum(!p))}
cat(sprintf("Stations %d, combinations %d\n\n",length(sel),nrow(combos)))
cat(sprintf("%5s %4s | %8s %7s | %8s %7s\n","W(h)","bw","var_med","varAUC","ar1_med","ar1AUC"))
OUT<-rbindlist(lapply(seq_len(nrow(combos)),function(ci){
  o<-acc[[ci]]$on;c_<-acc[[ci]]$ct
  data.table(W_h=combos$W[ci]/60,bw=combos$bw[ci],
             var_med=round(median(o[,1],na.rm=TRUE),3),
             var_auc=round(auc(o[,1],c_[,1]),3),
             ar1_med=round(median(o[,2],na.rm=TRUE),3),
             ar1_auc=round(auc(o[,2],c_[,2]),3))}))
OUT[,n_stations:=length(sel)]
fwrite(OUT,file.path(MART,"ews_sensitivity_9combo.csv"),bom=TRUE)
cat(sprintf("
Saved: 02_MART/ews_sensitivity_9combo.csv (stations: %d)
",length(sel)))
cat(sprintf("Discrimination AUC range — variance %.3f-%.3f · AR1 %.3f-%.3f
",
  min(OUT$var_auc),max(OUT$var_auc),min(OUT$ar1_auc),max(OUT$ar1_auc)))

for(ci in seq_len(nrow(combos))){o<-acc[[ci]]$on;c_<-acc[[ci]]$ct
  cat(sprintf("%5.0f %4.0f | %8.3f %7.3f | %8.3f %7.3f\n",
    combos$W[ci]/60,combos$bw[ci],
    median(o[,1],na.rm=T),auc(o[,1],c_[,1]),
    median(o[,2],na.rm=T),auc(o[,2],c_[,2])))}
