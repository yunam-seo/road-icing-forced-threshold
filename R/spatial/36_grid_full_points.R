# Builds the event-control time list for the grid analysis (controls in freezing months, at least 12 h from any onset).
# Inputs: 02_MART/intl_event_table.rds, road_grid500_match.csv.
# Output: 02_MART/grid_full_points.rds.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
set.seed(2026)
MART<-file.path(ROOT, "02_MART")
DDIR<-rwis_dir()
ctrl_ratio<-3L; buffer<-720L; freeze_mon<-c(11L,12L,1L,2L,3L)

mt<-fread(file.path(MART,"road_grid500_match.csv"))
ev<-as.data.table(readRDS(file.path(MART,"intl_event_table.rds")))
cat("All events:",nrow(ev)," stations:",uniqueN(ev$stn)," years:",
    paste(names(table(ev$year)),table(ev$year),sep="="),"\n")

map_cls<-function(s){o<-rep(NA_character_,length(s));o[s==1L]<-"DRY"
  o[s%in%c(2L,11:14)]<-"WET";o[s==3L]<-"SNOW";o[s%in%c(4L,15L)]<-"ICE";o}
fs<-list.files(DDIR,pattern="csv$",full.names=TRUE)
sof<-as.integer(sub(".*rwis_([0-9]+)_.*","\\1",basename(fs)))

out<-list()
for(sid in unique(ev$stn)){
  i<-which(sof==sid); if(!length(i))next
  d<-fread(fs[i], select=c("datetime","SFS","mon"), na.strings=c("-999",""))
  d<-d[mon %in% freeze_mon]
  d[,tm:=as.integer(as.POSIXct(strptime(sprintf("%012.0f",as.numeric(datetime)),
        "%Y%m%d%H%M",tz="UTC")))%/%60L]
  d<-d[!is.na(tm)]; d[,cls:=map_cls(SFS)]
  oo<-sort(ev[stn==sid]$onset_min)
  cand<-d[cls%in%c("DRY","WET","SNOW")]
  if(nrow(cand)&&length(oo)){
    nd<-sapply(cand$tm,function(t)min(abs(t-oo)))
    far<-cand[nd>=buffer]
    nc<-min(nrow(far),ctrl_ratio*length(oo))
    if(nc>0){cc<-far[sample(.N,nc)]
      out[[length(out)+1L]]<-data.table(stn=sid,tm_min=cc$tm,label=0L)}
  }
}
cc<-rbindlist(out)
pts<-rbind(data.table(stn=ev$stn,tm_min=ev$onset_min,label=1L), cc)
pts<-merge(pts, mt[,.(stn=STN_ID,xii,yii)], by="stn")
pts[,tm5:=(tm_min%/%5L)*5L]
cat("Cases:",sum(pts$label),"  controls:",sum(pts$label==0L),"  total:",nrow(pts),"\n")
pts[,yr:=as.integer(format(as.POSIXct(tm_min*60,origin="1970-01-01",tz="UTC"),"%Y"))]
cat("Cases by year:\n"); print(table(pts[label==1L]$yr))
saveRDS(pts, file.path(MART,"grid_full_points.rds"))
cat("Saved: 02_MART/grid_full_points.rds\n")
