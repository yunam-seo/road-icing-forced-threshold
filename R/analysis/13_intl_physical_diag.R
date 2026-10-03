# Computes meteorological diagnostics of the 24 h before each onset (wind, pressure tendency, precipitation, onset hour, night Ts-Ta, Ts range).
# Inputs: 02_MART/intl_event_table.rds, intl_cluster_assign.rds; 1-min RWIS files.
# Output: 02_MART/intl_physical_diag.csv.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
source(code_file("00_analysis_window.R"))
MART<-file.path(ROOT, "02_MART")
DDIR<-rwis_dir()
W<-1440L; fmon<-c(11L,12L,1L,2L,3L)

read_phys<-function(f){
  d<-fread(f,select=c("datetime","TR1","TA","WS","PA","RN05_YN","SD","TFZ","hour","mon"),
           na.strings=c("-999","-999.0",""))
  d<-d[mon%in%fmon];d<-filter_winter(d)
  if(!nrow(d))return(NULL)
  d[(TR1< -40)|(TR1>40),TR1:=NA_real_]; d[(TA< -40)|(TA>40),TA:=NA_real_]
  d[(WS<0)|(WS>60),WS:=NA_real_]; d[(PA<800)|(PA>1100),PA:=NA_real_]
  d[RN05_YN<0,RN05_YN:=NA_real_]
  d[,tm_min:=as.integer(as.POSIXct(strptime(sprintf("%012.0f",as.numeric(datetime)),
    "%Y%m%d%H%M",tz="UTC")))%/%60L];d<-d[!is.na(tm_min)];setorder(d,tm_min);setkey(d,tm_min)
  d[,.(tm_min,ts=TR1,ta=TA,ws=WS,pa=PA,rn=as.integer(RN05_YN>0),sd=SD,tfz=TFZ,hour)]
}

ev<-readRDS(file.path(MART,"intl_event_table.rds"));setkey(ev,stn)
files<-list.files(DDIR,pattern="\\.csv$",full.names=TRUE)
sof<-as.integer(sub(".*rwis_(\\d+)_.*","\\1",basename(files)))
args<-commandArgs(trailingOnly=TRUE);test<-length(args)>=1L&&args[1]=="test"
if(test)files<-head(files,8L)

out<-list();k<-0L
for(i in seq_along(files)){
  sid<-sof[i];evs<-ev[stn==sid];if(!nrow(evs))next
  d<-read_phys(files[i]);if(is.null(d)||nrow(d)<W)next
  for(j in seq_len(nrow(evs))){
    o<-evs$onset_min[j]
    win<-d[tm_min>=o-W & tm_min<o]
    if(nrow(win)<W*0.5)next
    at<-function(t)d[data.table(tm_min=t),on="tm_min",roll="nearest"]
    pa_o<-at(o)$pa; pa_6<-at(o-360L)$pa
    hr<-at(o)$hour
    night<-win[hour>=22L|hour<=7L]
    k<-k+1L
    out[[k]]<-data.table(eid=paste0(sid,"_",o),stn=sid,
      onset_hour=hr,
      ws_mean=mean(win$ws,na.rm=TRUE), ws_onset=at(o)$ws,
      pa_tend6=pa_o-pa_6,
      precip_any=as.integer(any(win$rn>0,na.rm=TRUE)),
      precip_frac=mean(win$rn>0,na.rm=TRUE),
      noct_ts_ta=mean(night$ts-night$ta,na.rm=TRUE),
      ts_amp=max(win$ts,na.rm=TRUE)-min(win$ts,na.rm=TRUE),
      tfz_ok=as.integer(any(is.finite(win$tfz))))
  }
  if(i%%20L==0L)message("  ",i,"/",length(files))
}
ph<-rbindlist(out,fill=TRUE)
asg<-readRDS(file.path(MART,"intl_cluster_assign.rds"))
ph<-merge(ph,asg[,.(eid,cluster)],by="eid")

mn<-function(x)round(mean(x,na.rm=TRUE),2);md<-function(x)round(median(x,na.rm=TRUE),1)
cat("\n=== Physical diagnostics by pathway (mean) ===\n")
print(ph[,.(n=.N,
  wind_speed=mn(ws_mean), pressure_tendency_6h=mn(pa_tend6),
  precipitation_any=mn(precip_any), precipitation_fraction=mn(precip_frac),
  night_ts_minus_ta=mn(noct_ts_ta), ts_amplitude=mn(ts_amp)),by=cluster][order(cluster)])
cat("\nOnset local hour (median/quartiles by cluster):\n")
print(ph[,.(median_hour=md(onset_hour),Q1=round(quantile(onset_hour,.25,na.rm=T)),
  Q3=round(quantile(onset_hour,.75,na.rm=T))),by=cluster][order(cluster)])
cat("\nOnset hour histogram (cluster × 6-h bins):\n")
ph[,hbin:=cut(onset_hour,c(-1,5,11,17,23),labels=c("0-5","6-11","12-17","18-23"))]
print(table(cluster=ph$cluster,hour_bin=ph$hbin))
cat("\nShare of events with TFZ (freezing point) available:",round(mean(ph$tfz_ok,na.rm=T),3),"\n")

fwrite(ph,file.path(MART,"intl_physical_diag.csv"))
cat("Saved: intl_physical_diag.csv\n")
