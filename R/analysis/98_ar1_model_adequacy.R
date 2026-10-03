# Checks whether detrended surface-temperature residuals follow an AR(1) process by comparing AR(k) with AR(1)^k.
# Inputs: 02_MART/intl_event_table.rds; 1-min RWIS files.
# Output: 02_MART/ar1_model_adequacy.csv.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
source(code_file("00_analysis_window.R"))
set.seed(2026)
MART <- file.path(ROOT, "02_MART")
DDIR <- rwis_dir()
ARGS<-commandArgs(trailingOnly=TRUE)
NSTN <- if(length(ARGS)>=1L) as.integer(ARGS[1]) else NA_integer_; W <- 720L; RW <- 360L; BW <- 60
LAGS <- c(1L, 2L, 3L, 5L, 10L, 20L)
cov_thr <- 0.8; ctrl_per <- 1L; ctrl_buf <- 1440L; fmon <- c(11L,12L,1L,2L,3L)

map_cls <- function(s){o<-rep(NA_character_,length(s)); o[s==1L]<-"DRY"
  o[s%in%c(2L,11:14)]<-"WET"; o[s==3L]<-"SNOW"; o[s%in%c(4L,15L)]<-"ICE"; o}
read_ts <- function(f){
  d<-fread(f,select=c("datetime","TR1","SFS","mon"),na.strings=c("-999","-999.0",""))
  d<-d[mon%in%fmon]; d<-filter_winter(d); if(!nrow(d))return(NULL)
  d[(TR1< -40)|(TR1>40),TR1:=NA_real_]
  d[,tm_min:=as.integer(as.POSIXct(strptime(sprintf("%012.0f",as.numeric(datetime)),
    "%Y%m%d%H%M",tz="UTC")))%/%60L]
  d<-d[!is.na(tm_min)]; setorder(d,tm_min)
  d[,cls:=map_cls(SFS)]; setkey(d,tm_min); d[,.(tm_min,ts=TR1,cls)]}
gtrend <- function(v,bw){ok<-which(!is.na(v))
  if(length(ok)<3) return(rep(NA_real_,length(v)))
  sapply(seq_along(v),function(i){w<-exp(-0.5*((ok-i)/bw)^2); sum(w*v[ok])/sum(w)})}

acf_end <- function(v){
  v<-tail(v,W+1L); if(mean(!is.na(v))<cov_thr) return(rep(NA_real_,length(LAGS)))
  r<-v-gtrend(v,BW); seg<-tail(r,RW)
  sapply(LAGS,function(L){
    a<-seg[1:(length(seg)-L)]; b<-seg[(1+L):length(seg)]
    g<-is.finite(a)&is.finite(b)
    if(sum(g)>50) suppressWarnings(cor(a[g],b[g])) else NA_real_})}

ev<-readRDS(file.path(MART,"intl_event_table.rds")); setkey(ev,stn)
files<-list.files(DDIR,pattern="\\.csv$",full.names=TRUE)
sof<-as.integer(sub(".*rwis_(\\d+)_.*","\\1",basename(files)))
sel<-which(sof%in%ev$stn)
if(!is.na(NSTN)) sel <- sel[1:min(NSTN,length(sel))]

ON<-CT<-list()
for(ii in sel){
  d<-try(read_ts(files[ii]),silent=TRUE)
  if(inherits(d,"try-error")||is.null(d)||nrow(d)<5000) next
  sid<-sof[ii]; evs<-ev[.(sid)][!is.na(onset_min)]; if(!nrow(evs)) next
  grab<-function(t_end){g<-data.table(tm_min=(t_end-W):t_end,key="tm_min"); d[g,on="tm_min"]$ts}
  for(t0 in evs$onset_min-1L) ON[[length(ON)+1L]]<-acf_end(grab(t0))
  oo<-sort(unique(evs$onset_min))
  cand<-d[cls%in%c("DRY","WET","SNOW") & tm_min>d$tm_min[1]+W]
  if(nrow(cand)){
    jj<-findInterval(cand$tm_min,oo)
    dl<-ifelse(jj>=1,abs(cand$tm_min-oo[pmax(jj,1L)]),Inf)
    dh<-ifelse(jj<length(oo),abs(cand$tm_min-oo[pmin(jj+1L,length(oo))]),Inf)
    far<-cand[pmin(dl,dh)>=ctrl_buf]
    nc<-min(nrow(far),ctrl_per*nrow(evs))
    if(nc>0) for(t0 in far[sample(.N,nc)]$tm_min) CT[[length(CT)+1L]]<-acf_end(grab(t0))}
  cat(sprintf("  station %d  onsets %d\n",sid,nrow(evs)))
}
mat<-function(L){m<-do.call(rbind,L); m[complete.cases(m),,drop=FALSE]}
MO<-mat(ON); MC<-mat(CT)
cat(sprintf("\nOnset windows %d · control windows %d\n\n",nrow(MO),nrow(MC)))

summ<-function(M,lab){
  a1<-median(M[,1])
  obs<-apply(M,2,median)
  pred<-a1^LAGS
  cat(sprintf("== %s (AR1=%.4f) ==\n",lab,a1))
  cat(sprintf("%5s %10s %10s %9s\n","lag","observed AR(k)","AR(1)^k","observed/predicted"))
  for(i in seq_along(LAGS))
    cat(sprintf("%5d %10.4f %10.4f %9.2f\n",LAGS[i],obs[i],pred[i],obs[i]/pred[i]))
  cat("\n")
  data.table(sample_group=lab, lag=LAGS, observed=round(obs,4), predicted=round(pred,4),
             ratio=round(obs/pred,2),
             implied_inv_lambda_min=round(-LAGS/log(obs),1))}
out<-rbind(summ(MO,"onset"), summ(MC,"control"))
out[, n_stations := length(sel)]
fwrite(out, file.path(MART,"ar1_model_adequacy.csv"))

for(g in c("onset","control")){
  a<-out[sample_group==g]; r<-a[lag==1]$implied_inv_lambda_min / a[lag==20]$implied_inv_lambda_min
  cat(sprintf("%s: relaxation time is %.1f min at lag1 and %.1f min at lag20 -> %.1f-fold difference
",
      g, a[lag==1]$implied_inv_lambda_min, a[lag==20]$implied_inv_lambda_min, r))}
ro<-out[sample_group=="onset"&lag==1]$implied_inv_lambda_min/out[sample_group=="onset"&lag==20]$implied_inv_lambda_min
rc<-out[sample_group=="control"&lag==1]$implied_inv_lambda_min/out[sample_group=="control"&lag==20]$implied_inv_lambda_min
cat(if(max(ro,rc)>1.25)
 "
Result: autocorrelation does not decay exponentially with lag.
  A relaxation time derived from lag-1 alone therefore cannot be read as a single recovery time.
  The departure is larger for onsets, so the two groups differ in the shape of their autocorrelation,
  and comparing them at a single lag does not compare the same quantity.
"
 else "
Result: consistent with exponential decay.
")
cat("Saved: 02_MART/ar1_model_adequacy.csv\n")
