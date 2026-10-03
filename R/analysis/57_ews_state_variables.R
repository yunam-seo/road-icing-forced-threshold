# Repeats the CSD test with alternative state variables (deposition margin, air temperature, friction coefficient, water-film thickness).
# Inputs: 02_MART/intl_event_table.rds; 1-min RWIS files.
# Outputs: 02_MART/ews_state_variables.csv, ews_state_variables.rds.
# Reads the full 1-min archive (slow).

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
sys.source(code_file("09_intl_ews.R"), envir=environment(), toplevel.env=environment())
set.seed(2026)
MART<-file.path(ROOT, "02_MART")

aw<-17.62; bw_<-243.12; ai<-22.46; bi<-272.62
frost_point<-function(ta,rh){
  e<-(rh/100)*6.112*exp(aw*ta/(bw_+ta))
  l<-log(pmax(e,1e-6)/6.112)
  bi*l/(ai-l)
}
read_vars<-function(f){
  d<-fread(f,select=c("datetime","TR1","TA","HM","SFS","mon","IFT","FRC","WFT"),
           na.strings=c("-999","-999.0",""))
  d<-d[mon%in%cfg$freeze_mon]; d<-filter_winter(d); if(!nrow(d))return(NULL)
  d[(TR1< -40)|(TR1>40),TR1:=NA_real_]; d[(TA< -40)|(TA>40),TA:=NA_real_]
  d[(HM<0)|(HM>100),HM:=NA_real_]
  d[(FRC<0)|(FRC>1),FRC:=NA_real_]
  d[(WFT<0)|(WFT>20),WFT:=NA_real_]
  d[,tm_min:=as.integer(as.POSIXct(strptime(sprintf("%012.0f",as.numeric(datetime)),
      "%Y%m%d%H%M",tz="UTC")))%/%60L]
  d<-d[!is.na(tm_min)]; setorder(d,tm_min)
  d[,cls:=map_sfs_class(SFS)]
  d[,tf:=frost_point(TA,HM)]
  d[,.(tm_min,cls,ts=TR1,ta=TA,m_dep=TR1-tf,frc=FRC,wft=WFT,ift=IFT)]
}
tau2<-function(v){
  if(mean(!is.na(v))<cfg$cov_thr) return(c(var_tau=NA_real_,ar1_tau=NA_real_))
  r<-v-gauss_trend(v,cfg$bw)
  st<-seq(1L,length(v)-cfg$rw,by=cfg$step); vv<-aa<-numeric(length(st))
  for(k in seq_along(st)){
    seg<-r[st[k]:(st[k]+cfg$rw-1L)]
    vv[k]<-var(seg,na.rm=TRUE)
    a<-seg[-length(seg)];b<-seg[-1];g<-is.finite(a)&is.finite(b)
    aa[k]<-if(sum(g)>30) suppressWarnings(cor(a[g],b[g])) else NA_real_}
  kt<-function(y){o<-is.finite(y);if(sum(o)<4)return(NA_real_)
    suppressWarnings(cor(seq_along(y)[o],y[o],method="kendall"))}
  c(var_tau=kt(vv), ar1_tau=kt(aa))
}

ev<-readRDS(file.path(MART,"intl_event_table.rds")); setkey(ev,stn)
files<-list.files(DDIR,pattern="[.]csv$",full.names=TRUE)
sof<-as.integer(sub(".*rwis_(\\d+)_.*","\\1",basename(files)))
VARS<-c("ts","m_dep","ta","frc","wft"); out<-list(); k<-0L; ift_sd<-c()

for(i in seq_along(files)){
  sid<-sof[i]; evs<-ev[stn==sid]; if(!nrow(evs))next
  d<-read_vars(files[i]); if(is.null(d)||nrow(d)<cfg$W)next
  setkey(d,tm_min); tmin<-d$tm_min[1]
  grab<-function(t_end,v){g<-data.table(tm_min=(t_end-cfg$W):t_end,key="tm_min"); d[g,on="tm_min"][[v]]}
  for(j in seq_len(nrow(evs))){
    te<-evs$onset_min[j]-1L; k<-k+1L
    row<-data.table(type="onset",stn=sid)
    for(v in VARS){e<-tau2(grab(te,v)); row[[paste0(v,"_var")]]<-e[[1]]; row[[paste0(v,"_ar1")]]<-e[[2]]}
    x<-grab(te,"ift"); ift_sd<-c(ift_sd, sd(x,na.rm=TRUE))
    out[[k]]<-row
  }
  oo<-sort(unique(evs$onset_min))
  cand<-d[cls%in%c("DRY","WET","SNOW")&tm_min>tmin+cfg$W]
  if(nrow(cand)){
    jj<-findInterval(cand$tm_min,oo)
    dl<-ifelse(jj>=1,abs(cand$tm_min-oo[pmax(jj,1L)]),Inf)
    dh<-ifelse(jj<length(oo),abs(cand$tm_min-oo[pmin(jj+1L,length(oo))]),Inf)
    far<-cand[pmin(dl,dh)>=cfg$ctrl_buffer]
    nc<-min(nrow(far),cfg$ctrl_per_event*nrow(evs))
    if(nc>0)for(te in far[sample(.N,nc)]$tm_min){
      k<-k+1L; row<-data.table(type="control",stn=sid)
      for(v in VARS){e<-tau2(grab(te,v)); row[[paste0(v,"_var")]]<-e[[1]]; row[[paste0(v,"_ar1")]]<-e[[2]]}
      out[[k]]<-row
    }
  }
  if(i%%20L==0L)message("  ",i,"/",length(files))
}
R<-rbindlist(out,fill=TRUE); saveRDS(R,file.path(MART,"ews_state_variables.rds"))

auc<-function(s,p){o<-is.finite(s);s<-s[o];p<-p[o]
  np<-as.numeric(sum(p));nn<-as.numeric(sum(!p)); if(!np||!nn)return(NA_real_)
  r<-rank(s);(sum(r[p])-np*(np+1)/2)/(np*nn)}
LAB<-c(ts="surface temperature Ts (baseline)", m_dep="deposition margin M_dep = Ts − Tf", ta="air temperature Ta",
       frc="friction coefficient FRC", wft="water-film thickness WFT")
res<-list()
for(v in VARS) for(ind in c("var","ar1")){
  cl<-paste0(v,"_",ind); x<-R[[cl]]; p<-R$type=="onset"
  res[[length(res)+1]]<-data.table(state_variable=LAB[[v]], indicator=ifelse(ind=="var","rolling variance","lag-1 autocorrelation"),
    AUC=round(auc(x,p),3), frac_positive_tau_onset=round(mean(x[p]>0,na.rm=TRUE),3),
    n_valid_windows_onset=sum(is.finite(x[p])), n_valid_windows_control=sum(is.finite(x[!p])))
}
S<-rbindlist(res); fwrite(S,file.path(MART,"ews_state_variables.csv"),bom=TRUE)
cat("\n=== 5 state variables × 2 indicators: onset vs control discrimination AUC ===\n"); print(S)
cat(sprintf("\nAUC range: %.3f–%.3f\n",min(S$AUC,na.rm=TRUE),max(S$AUC,na.rm=TRUE)))
cat(sprintf("Note: median IFT standard deviation in the window before onset = %.3f (near 0 → unusable as a state variable)\n",
            median(ift_sd,na.rm=TRUE)))
cat("Saved: 02_MART/ews_state_variables.csv\n")
