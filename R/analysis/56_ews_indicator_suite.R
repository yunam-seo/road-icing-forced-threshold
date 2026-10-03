# Evaluates six CSD indicators with random and time-of-day-matched control windows.
# Inputs: 02_MART/intl_event_table.rds; 1-min RWIS files.
# Outputs: 02_MART/ews_indicator_suite.csv, ews_indicator_suite.rds.
# Reads the full 1-min archive (slow).

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
sys.source(code_file("09_intl_ews.R"), envir=environment(), toplevel.env=environment())
set.seed(2026)
MART<-file.path(ROOT, "02_MART")
NSTN<-Inf

skew<-function(x){x<-x[is.finite(x)];n<-length(x);if(n<8)return(NA_real_)
  m<-mean(x);s<-sd(x);if(!is.finite(s)||s==0)return(NA_real_);sum((x-m)^3)/(n*s^3)}
kurt<-function(x){x<-x[is.finite(x)];n<-length(x);if(n<8)return(NA_real_)
  m<-mean(x);s<-sd(x);if(!is.finite(s)||s==0)return(NA_real_);sum((x-m)^4)/(n*s^4)-3}
redness<-function(x){x<-x[is.finite(x)];n<-length(x);if(n<32)return(NA_real_)
  p<-Mod(fft(x-mean(x)))^2; p<-p[2:floor(n/2)]
  if(!length(p)||sum(p)<=0)return(NA_real_)
  sum(p[1:max(1,floor(length(p)*0.2))])/sum(p)}
act<-function(x){x<-x[is.finite(x)];n<-length(x);if(n<60)return(NA_real_)
  a<-suppressWarnings(acf(x,lag.max=min(60,n-2),plot=FALSE)$acf[-1])
  w<-which(a<exp(-1))[1]; if(is.na(w)) length(a)+1 else w}

kt<-function(y){o<-is.finite(y);if(sum(o)<4)return(NA_real_)
  suppressWarnings(cor(seq_along(y)[o],y[o],method="kendall"))}

suite_window<-function(v){
  NAout<-setNames(rep(NA_real_,6),c("var","ar1","skew","kurt","red","act"))
  if(mean(!is.na(v))<cfg$cov_thr) return(NAout)
  r<-v-gauss_trend(v,cfg$bw)
  st<-seq(1L,length(v)-cfg$rw,by=cfg$step)
  M<-matrix(NA_real_,length(st),6)
  for(k in seq_along(st)){
    seg<-r[st[k]:(st[k]+cfg$rw-1L)]
    M[k,1]<-var(seg,na.rm=TRUE)
    a<-seg[-length(seg)];b<-seg[-1];g<-is.finite(a)&is.finite(b)
    M[k,2]<-if(sum(g)>30) suppressWarnings(cor(a[g],b[g])) else NA_real_
    M[k,3]<-skew(seg); M[k,4]<-kurt(seg); M[k,5]<-redness(seg); M[k,6]<-act(seg)
  }
  setNames(apply(M,2,kt),names(NAout))
}

ev<-readRDS(file.path(MART,"intl_event_table.rds")); setkey(ev,stn)
files<-list.files(DDIR,pattern="\\.csv$",full.names=TRUE)
sof<-as.integer(sub(".*rwis_(\\d+)_.*","\\1",basename(files)))
if(is.finite(NSTN)) files<-head(files,NSTN)
hour_of<-function(tm) as.integer(format(as.POSIXct(tm*60,origin="1970-01-01",tz="UTC"),"%H"))

out<-list(); k<-0L
for(i in seq_along(files)){
  sid<-sof[i]; evs<-ev[stn==sid]; if(!nrow(evs))next
  d<-read_ts(files[i]); if(is.null(d)||nrow(d)<cfg$W)next
  setkey(d,tm_min); tmin<-d$tm_min[1]
  grab<-function(t_end){g<-data.table(tm_min=(t_end-cfg$W):t_end,key="tm_min"); d[g,on="tm_min"]$ts}
  for(j in seq_len(nrow(evs))){
    e<-suite_window(grab(evs$onset_min[j]-1L)); k<-k+1L
    out[[k]]<-data.table(type="onset",ctrl="-",stn=sid,t(e))
  }
  oo<-sort(unique(evs$onset_min))
  cand<-d[cls%in%c("DRY","WET","SNOW")&tm_min>tmin+cfg$W]
  if(nrow(cand)){
    jj<-findInterval(cand$tm_min,oo)
    dl<-ifelse(jj>=1,abs(cand$tm_min-oo[pmax(jj,1L)]),Inf)
    dh<-ifelse(jj<length(oo),abs(cand$tm_min-oo[pmin(jj+1L,length(oo))]),Inf)
    far<-cand[pmin(dl,dh)>=cfg$ctrl_buffer]
    if(nrow(far)){
      nc<-min(nrow(far),cfg$ctrl_per_event*nrow(evs))
      for(t_end in far[sample(.N,nc)]$tm_min){
        e<-suite_window(grab(t_end)); k<-k+1L
        out[[k]]<-data.table(type="control",ctrl="random",stn=sid,t(e))
      }
      far[,hh:=hour_of(tm_min)]
      oh<-hour_of(evs$onset_min)
      for(j in seq_along(oh)){
        want<-((oh[j]+c(-1,0,1))%%24)
        pool<-far[hh%in%want]
        if(nrow(pool)) for(t_end in pool[sample(.N,min(.N,cfg$ctrl_per_event))]$tm_min){
          e<-suite_window(grab(t_end)); k<-k+1L
          out[[k]]<-data.table(type="control",ctrl="hour_matched",stn=sid,t(e))
        }
      }
    }
  }
  if(i%%20L==0L)message("  ",i,"/",length(files))
}
R<-rbindlist(out,fill=TRUE)
saveRDS(R,file.path(MART,"ews_indicator_suite.rds"))

auc<-function(s,p){o<-is.finite(s);s<-s[o];p<-p[o]
  np<-as.numeric(sum(p));nn<-as.numeric(sum(!p)); if(!np||!nn)return(NA_real_)
  r<-rank(s);(sum(r[p])-np*(np+1)/2)/(np*nn)}
IND<-c(var="rolling variance",ar1="lag-1 autocorrelation",skew="skewness",kurt="kurtosis",red="spectral reddening",act="autocorrelation decay time")
res<-list()
for(cc in c("random","hour_matched")){
  sub<-R[type=="onset"|(type=="control"&ctrl==cc)]
  for(v in names(IND)){
    x<-sub[[v]]; p<-sub$type=="onset"
    res[[length(res)+1]]<-data.table(control=cc, indicator=IND[[v]], code=v,
      onset_median=round(median(x[p],na.rm=TRUE),3),
      control_median=round(median(x[!p],na.rm=TRUE),3),
      AUC=round(auc(x,p),3),
      frac_positive_tau_onset=round(mean(x[p]>0,na.rm=TRUE),3),
      n_onset=sum(is.finite(x[p])), n_control=sum(is.finite(x[!p])))
  }
}
S<-rbindlist(res); fwrite(S,file.path(MART,"ews_indicator_suite.csv"),bom=TRUE)
cat("\n=== 6 CSD indicators × 2 control designs: onset vs control discrimination AUC ===\n"); print(S)
cat("\nAUC range:",sprintf("%.3f–%.3f",min(S$AUC,na.rm=TRUE),max(S$AUC,na.rm=TRUE)),"\n")
cat("Saved: 02_MART/ews_indicator_suite.csv\n")
