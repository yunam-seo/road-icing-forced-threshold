# Extracts surrounding-field components (station cell, 10 km disc mean, upwind and downwind cells) from the 500 m objectively analysed fields at every nowcast evaluation time.
# Inputs: 02_MART/intl_ew_features.rds, road_grid500_match.csv; 500 m NetCDF fields.
# Output: 02_MART/grid_ew_scale_raw.rds.
# Reads the full grid archive (slow); keeps checkpoints and a wind cache in 02_MART/_ckpt and resumes from them.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages({library(data.table);library(terra);library(parallel)})
MART<-file.path(ROOT, "02_MART")
AGB<-grid_fields()
DDIR<-rwis_dir()
Rreg<-20L; CELLKM<-0.5; Dkm<-c(5L,10L,20L); Dc<-as.integer(Dkm/CELLKM)
args<-commandArgs(trailingOnly=TRUE); TEST<-length(args)>=1L && args[1]=="test"
NW<-6L

ew<-as.data.table(readRDS(file.path(MART,"intl_ew_features.rds")))[,.(stn,tm_min)]
mt<-fread(file.path(MART,"road_grid500_match.csv"))
P<-merge(ew, mt[,.(stn=STN_ID,xii,yii)], by="stn", sort=FALSE)
stopifnot(nrow(P)==nrow(ew))
if(TEST){ set.seed(1); P<-P[sample(.N,300)] }
P[,tm5:=(tm_min%/%5L)*5L]
N<-nrow(P); cat("Evaluation points:",N," stations:",uniqueN(P$stn),"\n")

CKROOT<-Sys.getenv("RICE_CKPT", file.path(MART,"_ckpt"))
WCACHE<-if(TEST) "" else file.path(CKROOT,"scale_wind_cache.rds")
t_w<-Sys.time()
if(nzchar(WCACHE) && file.exists(WCACHE)){ wc<-readRDS(WCACHE); stopifnot(identical(wc$stn,P$stn), identical(wc$tm_min,P$tm_min))
  P[,`:=`(wd0=wc$wd0, ws0=wc$ws0, wind_gap_min=wc$wind_gap_min)]; cat("Wind direction: using cache\n")
} else {
fs<-list.files(DDIR,pattern="csv$",full.names=TRUE)
sof<-as.integer(sub(".*rwis_([0-9]+)_.*","\\1",basename(fs)))
wl<-list()
for(sid in unique(P$stn)){ i<-which(sof==sid); if(!length(i))next
  d<-fread(fs[i],select=c("datetime","WS","WD"),na.strings=c("-999",""))
  d[,tm:=as.integer(as.POSIXct(strptime(sprintf("%012.0f",as.numeric(datetime)),
        "%Y%m%d%H%M",tz="UTC")))%/%60L]
  d<-d[!is.na(tm)&!is.na(WD),.(stn=sid,tm,wtm=tm,ws=WS,wd=WD)]
  wl[[as.character(sid)]]<-d; rm(d) }
windall<-rbindlist(wl); rm(wl); gc(verbose=FALSE); setkey(windall,stn,tm)
wj<-windall[data.table(stn=P$stn,tm=P$tm_min), roll="nearest", on=.(stn,tm)]
P[,`:=`(wd0=wj$wd, ws0=wj$ws, wind_gap_min=abs(wj$wtm-P$tm_min))]
if(nzchar(WCACHE)){ dir.create(dirname(WCACHE),showWarnings=FALSE,recursive=TRUE); saveRDS(P[,.(stn,tm_min,wd0,ws0,wind_gap_min)],WCACHE) } }
cat("Attaching wind direction:",round(difftime(Sys.time(),t_w,units="mins"),1),"min\n")
P[,fx:=sin(wd0*pi/180)]; P[,fy:=cos(wd0*pi/180)]
cat("Wind direction missing:",sum(is.na(P$wd0))," | nearest wind observation >30 min away:",sum(P$wind_gap_min>30,na.rm=TRUE),"\n")

P[,row:=.I]
need0<-P[,.(row,t=tm5)]; need30<-P[,.(row,t=tm5-30L)]
ts_list<-sort(unique(c(need0$t,need30$t)))
cat("Unique times:",length(ts_list),"\n")
CK<-file.path(CKROOT, if(TEST) "scale_test" else "scale")
dir.create(CK,showWarnings=FALSE,recursive=TRUE); if(nzchar(WCACHE)) dir.create(dirname(WCACHE),showWarnings=FALSE,recursive=TRUE)
CH<-split(ts_list, ceiling(seq_along(ts_list)/250))
IDX0<-split(need0$row, need0$t); IDX30<-split(need30$row, need30$t)
PX<-P[,.(xii,yii,fx,fy)]
for(o in c("windall","wj","ew","wc")) if(exists(o)) rm(list=o); gc(verbose=FALSE)
cat("Chunks:",length(CH),"(250 times per chunk)  chunks already done:",sum(file.exists(file.path(CK,sprintf("ch%04d.rds",seq_along(CH))))),"\n")

work<-function(ci){
  suppressPackageStartupMessages({library(data.table);library(terra)}); terraOptions(memfrac=0.1)
  out<-file.path(CK,sprintf("ch%04d.rds",ci)); if(file.exists(out)) return(0L)
  off<-CJ(dy=-Rreg:Rreg,dx=-Rreg:Rreg)[sqrt(dy^2+dx^2)<=Rreg]
  tspath<-function(t,stem){o<-as.POSIXct(t*60,origin="1970-01-01",tz="UTC")
    file.path(AGB,format(o,"%Y/%m/%d"),paste0("sfc_grid_",stem,"_",format(o,"%Y%m%d%H%M"),".nc"))}
  samp<-function(m,y,x){ ok<-!is.na(y)&!is.na(x)&y>=1&x>=1&y<=nrow(m)&x<=ncol(m)
    v<-rep(NA_real_,length(y)); v[ok]<-m[cbind(y[ok],x[ok])]; v }
  R0<-list(); R30<-list(); nmiss<-0L
  for(t in CH[[ci]]){
    kt<-as.character(t); r0<-IDX0[[kt]]; r30<-IDX30[[kt]]
    for(vn in c("ta","td")){ f<-tspath(t,vn)
      if(!file.exists(f)){ nmiss<-nmiss+1L; next }
      m<-as.matrix(rast(f),wide=TRUE); m[m< -9000]<-NA
      if(length(r0)){ x<-PX$xii[r0]; y<-PX$yii[r0]; fx<-PX$fx[r0]; fy<-PX$fy[r0]
        yy<-pmin(pmax(outer(y,off$dy,"+"),1L),nrow(m)); xx<-pmin(pmax(outer(x,off$dx,"+"),1L),ncol(m))
        vr<-matrix(m[cbind(as.vector(yy),as.vector(xx))],nrow=length(r0))
        d<-data.table(row=r0,var=vn,vc=samp(m,y,x),vreg=rowMeans(vr,na.rm=TRUE))
        for(k in seq_along(Dkm)){
          d[[paste0("up",Dkm[k])]]<-samp(m,round(y+Dc[k]*fy),round(x+Dc[k]*fx))
          d[[paste0("dn",Dkm[k])]]<-samp(m,round(y-Dc[k]*fy),round(x-Dc[k]*fx)) }
        R0[[length(R0)+1L]]<-d }
      if(length(r30)) R30[[length(R30)+1L]]<-data.table(row=r30,var=vn,vc30=samp(m,PX$yii[r30],PX$xii[r30]))
    } }
  saveRDS(list(R0=rbindlist(R0),R30=rbindlist(R30),nmiss=nmiss),out); 1L }

t0<-Sys.time()
cl<-makeCluster(NW)
clusterExport(cl,c("CH","CK","PX","IDX0","IDX30","Rreg","Dkm","Dc","AGB"))
todo<-which(!file.exists(file.path(CK,sprintf("ch%04d.rds",seq_along(CH)))))
done<-0L
for(b in split(todo, ceiling(seq_along(todo)/NW))){
  parLapply(cl,b,work); done<-done+length(b)
  message(sprintf("  chunk %d/%d  %.1f min", done, length(todo), as.numeric(difftime(Sys.time(),t0,units="mins")))) }
stopCluster(cl)

L<-lapply(file.path(CK,sprintf("ch%04d.rds",seq_along(CH))),readRDS)
R0<-rbindlist(lapply(L,`[[`,"R0")); R30<-rbindlist(lapply(L,`[[`,"R30"))
nmiss<-sum(sapply(L,`[[`,"nmiss"))
w0<-dcast(R0, row~var, value.var=c("vc","vreg",paste0("up",Dkm),paste0("dn",Dkm)))
w30<-dcast(R30, row~var, value.var="vc30")
setnames(w0, names(w0), sub("^(vc|vreg|up5|up10|up20|dn5|dn10|dn20)_(ta|td)$","\\2_\\1",names(w0)))
setnames(w30, c("ta","td"), c("ta_vc30","td_vc30"))
setnames(w0, c("ta_vc","td_vc","ta_vreg","td_vreg"), c("ta_vc0","td_vc0","ta_vreg0","td_vreg0"))
O<-Reduce(function(a,b)merge(a,b,by="row",all.x=TRUE),
          list(P[,.(row,stn,tm_min,wd0,ws0,wind_gap_min)], w0, w30))
setorder(O,row); O[,row:=NULL]
stopifnot(nrow(O)==N)
cat("Missing grid files:",nmiss,"\n")
cat("Missing rate (%): ta_vc0",round(100*mean(is.na(O$ta_vc0)),3)," ta_vc30",round(100*mean(is.na(O$ta_vc30)),3),
    " ta_vreg0",round(100*mean(is.na(O$ta_vreg0)),3)," td_up10",round(100*mean(is.na(O$td_up10)),3),"\n")
outf<-file.path(MART, if(TEST) "grid_ew_scale_raw_TEST.rds" else "grid_ew_scale_raw.rds")
saveRDS(O,outf)
cat("Saved:",outf," elapsed:",round(difftime(Sys.time(),t0,units="mins"),1),"min\n")
