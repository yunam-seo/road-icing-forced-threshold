# Extracts gridded air temperature and dew point around each event and control point, with RWIS state after QC.
# Inputs: 02_MART/grid_full_points.rds; 500 m NetCDF fields; 1-min RWIS files.
# Output: 02_MART/grid_full_raw.rds.
# Reads the full grid archive (slow).

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages({library(data.table);library(terra)})
set.seed(2026)
MART<-file.path(ROOT, "02_MART")
AGB<-grid_fields()
DDIR<-rwis_dir()
gvars<-c(ta="ta", td="td"); lags<-c(0L,30L); Rcell<-6L; Dup<-10L
ts_range<-c(-40,40); ta_range<-c(-40,40); hm_range<-c(0,100); spike_max<-10

pts<-readRDS(file.path(MART,"grid_full_points.rds"))
on<-pts[label==1L]; non<-on[,.(non=.N),by=stn]
ct<-merge(pts[label==0L], non, by="stn")[, .SD[sample(.N, min(.N, non[1]))], by=stn]
ct[,non:=NULL]
P<-rbind(on, ct, fill=TRUE); setkey(P,stn)
args<-commandArgs(trailingOnly=TRUE)
if(length(args)>=1L && args[1]=="test") P<-P[sample(.N, 40)]
cat("Cases",sum(P$label),"  controls",sum(P$label==0L)," total",nrow(P),"\n")

off<-CJ(dy=-Rcell:Rcell, dx=-Rcell:Rcell)[sqrt(dy^2+dx^2)<=Rcell]
cat("Circular neighbourhood cells:",nrow(off),"\n")

qc_station<-function(d){
  setorder(d,tm)
  d[(ts<ts_range[1])|(ts>ts_range[2]), ts:=NA_real_]
  d[(ta<ta_range[1])|(ta>ta_range[2]), ta:=NA_real_]
  d[(rh<hm_range[1])|(rh>hm_range[2]), rh:=NA_real_]
  d[,dp:=ts-data.table::shift(ts)]; d[,dn:=data.table::shift(ts,type="lead")-ts]
  d[!is.na(dp)&!is.na(dn)&abs(dp)>spike_max&abs(dn)>spike_max&sign(dp)!=sign(dn), ts:=NA_real_]
  d[,c("dp","dn"):=NULL]; d }
fs<-list.files(DDIR,pattern="csv$",full.names=TRUE)
sof<-as.integer(sub(".*rwis_([0-9]+)_.*","\\1",basename(fs)))
rwis<-list()
for(sid in unique(P$stn)){ i<-which(sof==sid); if(!length(i))next
  d<-fread(fs[i],select=c("datetime","TR1","TA","HM","WS","WD"),na.strings=c("-999",""))
  d[,tm:=as.integer(as.POSIXct(strptime(sprintf("%012.0f",as.numeric(datetime)),
        "%Y%m%d%H%M",tz="UTC")))%/%60L]
  d<-d[!is.na(tm), .(tm,ts=TR1,ta=TA,rh=HM,ws=WS,wd=WD)]
  rwis[[as.character(sid)]]<-qc_station(d); setkey(rwis[[as.character(sid)]],tm) }

P[,tm5:=(tm_min%/%5L)*5L]
need<-CJ(eid=seq_len(nrow(P)), lag=lags)[, .(row=eid, lag, ts=P$tm5[eid]-lag)]
ts_list<-sort(unique(need$ts))
cat("Unique times:",length(ts_list)," × ",length(gvars),"variables =",length(ts_list)*length(gvars),"files\n")
tspath<-function(t,stem){o<-as.POSIXct(t*60,origin="1970-01-01",tz="UTC")
  file.path(AGB,format(o,"%Y/%m/%d"),paste0("sfc_grid_",stem,"_",format(o,"%Y%m%d%H%M"),".nc"))}

SC<-vector("list", length(ts_list)*length(gvars)); k<-0L; t0<-Sys.time(); nmiss<-0L
rc<-data.table(row=seq_len(nrow(P)), xii=P$xii, yii=P$yii, stn=P$stn)
for(ti in seq_along(ts_list)){ t<-ts_list[ti]
  rows_t<-need[ts==t]; if(!nrow(rows_t))next
  rr<-unique(rows_t$row); cell<-rc[row %in% rr]
  for(vn in names(gvars)){ f<-tspath(t,gvars[[vn]])
    if(!file.exists(f)){nmiss<-nmiss+1L; next}
    m<-as.matrix(rast(f),wide=TRUE); m[m< -9000]<-NA
    vc<-vreg<-gx<-gy<-lap<-numeric(nrow(cell))
    for(j in seq_len(nrow(cell))){ x<-cell$xii[j]; y<-cell$yii[j]
      vc[j]<-m[y,x]
      vreg[j]<-mean(m[cbind(y+off$dy, x+off$dx)],na.rm=TRUE)
      gx[j]<-(m[y,x+1]-m[y,x-1])/2; gy[j]<-(m[y+1,x]-m[y-1,x])/2
      lap[j]<-m[y+1,x]+m[y-1,x]+m[y,x+1]+m[y,x-1]-4*m[y,x] }
    k<-k+1L; SC[[k]]<-data.table(row=cell$row, lag=rows_t$lag[match(cell$row,rows_t$row)],
       gv=vn, t=t, xii=cell$xii, yii=cell$yii, vc,vreg,gx,gy,lap) }
  if(ti%%500L==0L) message("  ",ti,"/",length(ts_list)," (",
     round(difftime(Sys.time(),t0,units="mins"),1),"min)") }
S<-rbindlist(SC)
cat("Missing grid files:",nmiss," / extracted scalar rows:",nrow(S),"\n")
saveRDS(list(P=P, S=S, rwis=rwis, cfg=list(Rcell=Rcell,Dup=Dup,lags=lags,qc=TRUE)),
        file.path(MART,"grid_full_raw.rds"))
cat("Saved: 02_MART/grid_full_raw.rds  elapsed:",
    round(difftime(Sys.time(),t0,units="mins"),1),"min\n")
