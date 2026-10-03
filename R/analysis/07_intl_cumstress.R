# Computes leaky-integrated cold dose, deposition dose and freeze-thaw accumulation at icing onsets and control times.
# Inputs: 02_MART/intl_event_table.rds, intl_cluster_assign.rds; 1-min RWIS files.
# Output: 02_MART/intl_cumstress.rds.
# Reads the full 1-min archive (slow).

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
source(code_file("00_analysis_window.R"))
set.seed(2026)
MART <- file.path(ROOT, "02_MART")
DDIR <- rwis_dir()
cfg <- list(tau=1440, gap_max=15L, freeze_mon=c(11L,12L,1L,2L,3L),
            ctrl_per_event=5L, ctrl_buffer=720L,
            ts_range=c(-40,40), ta_range=c(-40,40), hm_range=c(0,100))

map_sfs_class <- function(s){o<-rep(NA_character_,length(s));o[s==1L]<-"DRY"
  o[s%in%c(2L,11:14)]<-"WET";o[s==3L]<-"SNOW";o[s%in%c(4L,15L)]<-"ICE";o}
frost_point <- function(Ta,RH){aw<-17.62;bw<-243.12;ai<-22.46;bi<-272.62
  e<-(pmax(RH,1e-3)/100)*6.112*exp(aw*Ta/(bw+Ta));l<-log(e/6.112);bi*l/(ai-l)}

read_station <- function(f){
  d<-fread(f,select=c("datetime","STN_ID","TR1","TA","HM","SFS","mon"),
           na.strings=c("-999","-999.0",""))
  d<-d[mon%in%cfg$freeze_mon]; d<-filter_winter(d)
  if(!nrow(d))return(NULL)
  d[(TR1<cfg$ts_range[1])|(TR1>cfg$ts_range[2]),TR1:=NA_real_]
  d[(TA<cfg$ta_range[1])|(TA>cfg$ta_range[2]),TA:=NA_real_]
  d[(HM<cfg$hm_range[1])|(HM>cfg$hm_range[2]),HM:=NA_real_]
  d[,tm_min:=as.integer(as.POSIXct(strptime(sprintf("%012.0f",as.numeric(datetime)),
            "%Y%m%d%H%M",tz="UTC")))%/%60L]
  d<-d[!is.na(tm_min)]; setorder(d,tm_min)
  d[,cls:=map_sfs_class(SFS)]
  d[,.(tm_min,ts=TR1,ta=TA,rh=HM,cls)]
}

leaky <- function(d,tau,gap_max){
  n<-nrow(d); cold<-dep<-ft<-numeric(n)
  Sc<-Sd<-Sf<-0; prev_t<-d$tm_min[1]; prev_sign<-NA
  fp<-frost_point(d$ta,d$rh)
  xc<-pmax(0,-d$ts); xd<-pmax(0,-(d$ts-fp))
  for(i in seq_len(n)){
    dt<-d$tm_min[i]-prev_t; dec<-exp(-dt/tau)
    Sc<-Sc*dec; Sd<-Sd*dec; Sf<-Sf*dec
    if(!is.na(dt)&&dt<=gap_max){
      if(!is.na(xc[i])) Sc<-Sc+xc[i]*dt
      if(!is.na(xd[i])) Sd<-Sd+xd[i]*dt
      s<-if(!is.na(d$ts[i])) sign(d$ts[i]) else NA
      if(!is.na(s)&&s!=0){
        if(!is.na(prev_sign)&&s!=prev_sign) Sf<-Sf+1
        prev_sign<-s
      }
    } else {
      prev_sign<-NA
    }
    cold[i]<-Sc; dep[i]<-Sd; ft[i]<-Sf; prev_t<-d$tm_min[i]
  }
  data.table(tm_min=d$tm_min,cold_dose=cold,dep_dose=dep,ft=ft,cls=d$cls)
}

build <- function(nfile=NULL){
  ev<-readRDS(file.path(MART,"intl_event_table.rds")); setkey(ev,stn)
  files<-list.files(DDIR,pattern="\\.csv$",full.names=TRUE)
  if(!is.null(nfile)) files<-head(files,nfile)
  sof<-as.integer(sub(".*rwis_(\\d+)_.*","\\1",basename(files)))
  out<-list(); k<-0L
  for(i in seq_along(files)){
    sid<-sof[i]; evs<-ev[stn==sid]; if(!nrow(evs))next
    d<-read_station(files[i]); if(is.null(d)||nrow(d)<10L)next
    L<-leaky(d,cfg$tau,cfg$gap_max); setkey(L,tm_min)
    on<-L[data.table(tm_min=evs$onset_min,key="tm_min"),roll="nearest",on="tm_min"]
    on[,`:=`(type="onset",stn=sid,onset_min=evs$onset_min)]
    cand<-L[cls%in%c("DRY","WET","SNOW")]
    oo<-sort(unique(evs$onset_min))
    j<-findInterval(cand$tm_min,oo)
    dlo<-ifelse(j>=1L,abs(cand$tm_min-oo[pmax(j,1L)]),Inf)
    dhi<-ifelse(j<length(oo),abs(cand$tm_min-oo[pmin(j+1L,length(oo))]),Inf)
    far<-cand[pmin(dlo,dhi)>=cfg$ctrl_buffer]
    nc<-min(nrow(far),cfg$ctrl_per_event*nrow(evs))
    if(nc>0){
      cc<-far[sample(.N,nc)]; cc[,`:=`(type="control",stn=sid,onset_min=NA_integer_)]
      k<-k+1L; out[[k]]<-rbind(on,cc,fill=TRUE)
    } else { k<-k+1L; out[[k]]<-on }
    if(i%%20L==0L) message("  ",i,"/",length(files))
  }
  rbindlist(out,fill=TRUE)
}

if(sys.nframe()==0L){
  args<-commandArgs(trailingOnly=TRUE)
  test<-length(args)>=1L && args[1]=="test"
  t0<-Sys.time()
  res<-build(if(test)3L else NULL)
  cat("Build time:",round(difftime(Sys.time(),t0,units="secs"),1),"s\n")
  asg<-readRDS(file.path(MART,"intl_cluster_assign.rds"))
  res[type=="onset",eid:=paste0(stn,"_",onset_min)]
  res<-merge(res,asg[,.(eid,cluster)],by="eid",all.x=TRUE)

  auc<-function(score,pos){ r<-rank(score); (sum(r[pos])-sum(pos)*(sum(pos)+1)/2)/
      (sum(pos)*sum(!pos)) }
  pos<-res$type=="onset"
  cat("Sample: onset",sum(pos)," control",sum(!pos),"\n")
  for(v in c("cold_dose","dep_dose","ft"))
    cat(sprintf("  AUC %s: %.3f  (onset median %.0f vs control %.0f)\n",v,
        auc(res[[v]],pos), median(res[[v]][pos],na.rm=T),
        median(res[[v]][!pos],na.rm=T)))
  cat("\nOnset doses by pathway (median):\n")
  print(res[type=="onset"&!is.na(cluster),
    .(n=.N,cold=round(median(cold_dose)),dep=round(median(dep_dose)),
      ft=round(median(ft),2)),by=cluster][order(cluster)])
  if(!test){ saveRDS(res,file.path(MART,"intl_cumstress.rds"))
    cat("\nSaved: intl_cumstress.rds\n") }
}
