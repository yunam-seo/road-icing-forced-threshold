# Measures how the CSD discrimination AUC of friction coefficient and surface temperature changes as the window end is moved before onset.
# Inputs: 02_MART/intl_event_table.rds; 1-min RWIS files.
# Outputs: 02_MART/frc_variance_lead.csv, frc_variance_lead.rds.
# Reads the full 1-min archive (slow).

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
sys.source(code_file("09_intl_ews.R"), envir=environment(), toplevel.env=environment())
set.seed(2026)
MART<-file.path(ROOT, "02_MART")
LEADS<-c(0L,30L,60L,120L,180L,360L,720L)

read_vars<-function(f){
  d<-fread(f,select=c("datetime","TR1","SFS","mon","FRC"),na.strings=c("-999","-999.0",""))
  d<-d[mon%in%cfg$freeze_mon]; d<-filter_winter(d); if(!nrow(d))return(NULL)
  d[(TR1< -40)|(TR1>40),TR1:=NA_real_]
  d[(FRC<0)|(FRC>1),FRC:=NA_real_]
  d[,tm_min:=as.integer(as.POSIXct(strptime(sprintf("%012.0f",as.numeric(datetime)),
      "%Y%m%d%H%M",tz="UTC")))%/%60L]
  d<-d[!is.na(tm_min)]; setorder(d,tm_min)
  d[,cls:=map_sfs_class(SFS)]
  d[,.(tm_min,cls,ts=TR1,frc=FRC)]
}
tau2<-function(v){
  if(mean(!is.na(v))<cfg$cov_thr) return(c(NA_real_,NA_real_))
  r<-v-gauss_trend(v,cfg$bw)
  st<-seq(1L,length(v)-cfg$rw,by=cfg$step); vv<-aa<-numeric(length(st))
  for(k in seq_along(st)){
    seg<-r[st[k]:(st[k]+cfg$rw-1L)]
    vv[k]<-var(seg,na.rm=TRUE)
    a<-seg[-length(seg)];b<-seg[-1];g<-is.finite(a)&is.finite(b)
    aa[k]<-if(sum(g)>30) suppressWarnings(cor(a[g],b[g])) else NA_real_}
  kt<-function(y){o<-is.finite(y);if(sum(o)<4)return(NA_real_)
    suppressWarnings(cor(seq_along(y)[o],y[o],method="kendall"))}
  c(kt(vv),kt(aa))
}

ev<-readRDS(file.path(MART,"intl_event_table.rds")); setkey(ev,stn)
files<-list.files(DDIR,pattern="[.]csv$",full.names=TRUE)
sof<-as.integer(sub(".*rwis_(\\d+)_.*","\\1",basename(files)))
out<-list(); k<-0L
for(i in seq_along(files)){
  sid<-sof[i]; evs<-ev[stn==sid]; if(!nrow(evs))next
  d<-read_vars(files[i]); if(is.null(d)||nrow(d)<cfg$W)next
  setkey(d,tm_min); tmin<-d$tm_min[1]
  grab<-function(t_end,v){g<-data.table(tm_min=(t_end-cfg$W):t_end,key="tm_min"); d[g,on="tm_min"][[v]]}
  addrow<-function(type,t0){
    row<-data.table(type=type,stn=sid)
    for(L in LEADS){ te<-t0-L
      for(v in c("frc","ts")){ e<-tau2(grab(te,v))
        row[[sprintf("%s_var_L%d",v,L)]]<-e[1]; row[[sprintf("%s_ar1_L%d",v,L)]]<-e[2] } }
    row }
  for(j in seq_len(nrow(evs))){ k<-k+1L; out[[k]]<-addrow("onset",evs$onset_min[j]-1L) }
  oo<-sort(unique(evs$onset_min))
  cand<-d[cls%in%c("DRY","WET","SNOW")&tm_min>tmin+cfg$W+max(LEADS)]
  if(nrow(cand)){
    jj<-findInterval(cand$tm_min,oo)
    dl<-ifelse(jj>=1,abs(cand$tm_min-oo[pmax(jj,1L)]),Inf)
    dh<-ifelse(jj<length(oo),abs(cand$tm_min-oo[pmin(jj+1L,length(oo))]),Inf)
    far<-cand[pmin(dl,dh)>=cfg$ctrl_buffer]
    nc<-min(nrow(far),cfg$ctrl_per_event*nrow(evs))
    if(nc>0)for(te in far[sample(.N,nc)]$tm_min){ k<-k+1L; out[[k]]<-addrow("control",te) }
  }
  if(i%%20L==0L)message("  ",i,"/",length(files))
}
R<-rbindlist(out,fill=TRUE); saveRDS(R,file.path(MART,"frc_variance_lead.rds"))

auc<-function(s,p){o<-is.finite(s);s<-s[o];p<-p[o]
  np<-as.numeric(sum(p));nn<-as.numeric(sum(!p)); if(!np||!nn)return(NA_real_)
  r<-rank(s);(sum(r[p])-np*(np+1)/2)/(np*nn)}
p<-R$type=="onset"; res<-list()
for(L in LEADS) for(v in c("frc","ts")) for(ind in c("var","ar1")){
  x<-R[[sprintf("%s_%s_L%d",v,ind,L)]]
  res[[length(res)+1]]<-data.table(
    variable=ifelse(v=="frc","friction coefficient FRC","surface temperature Ts"),
    indicator=ifelse(ind=="var","rolling variance","lag-1 autocorrelation"),
    lead_min=L, AUC=round(auc(x,p),3),
    n_valid_windows_onset=sum(is.finite(x[p])))
}
S<-rbindlist(res); fwrite(S,file.path(MART,"frc_variance_lead.csv"),bom=TRUE)
cat("\n=== Discrimination AUC with the window end shifted L min before onset ===\n")
print(dcast(S,variable+indicator~lead_min,value.var="AUC"))
cat("\nSaved: 02_MART/frc_variance_lead.csv\n")
