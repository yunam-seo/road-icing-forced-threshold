# Extracts the same surrounding-field components for the event-control sample at the evaluation time and 1 h and 3 h earlier.
# Inputs: 02_MART/grid_full_raw.rds, grid_upwind_raw.rds; 500 m NetCDF fields.
# Output: 02_MART/grid_past_scale_raw.rds.
# Reads the full grid archive (slow).

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages({library(data.table);library(terra);library(parallel)})
MART<-file.path(ROOT, "02_MART")
AGB<-grid_fields()
Rreg<-20L; CELLKM<-0.5; Dkm<-c(5L,10L,20L); Dc<-as.integer(Dkm/CELLKM)
LAGS<-c(0L,30L,60L,90L,180L,210L)
args<-commandArgs(trailingOnly=TRUE); TEST<-length(args)>=1L && args[1]=="test"
NW<-6L

x<-readRDS(file.path(MART,"grid_full_raw.rds")); P<-as.data.table(x$P); rwis<-x$rwis; rm(x); gc()
P[,row:=.I]
U<-as.data.table(readRDS(file.path(MART,"grid_upwind_raw.rds")))
stopifnot(nrow(U)==nrow(P), all(U$stn==P$stn), all(U$tm_min==P$tm_min), all(U$label==P$label))
cat("Sample: events",sum(P$label),"control",sum(P$label==0L),"— rows match the output of script 47\n")
if(TEST){ set.seed(1); P<-P[sample(.N,40)] }

Q<-rbindlist(lapply(LAGS,function(L) P[,.(row,stn,label,tm_min,xii,yii,lag=L,tq=tm_min-L)]))
Q[,t5:=(tq%/%5L)*5L]
windall<-rbindlist(lapply(names(rwis),function(s) rwis[[s]][!is.na(wd),.(stn=as.integer(s),tm,wtm=tm,wd,ws)]))
setkey(windall,stn,tm)
wv<-windall[data.table(stn=Q$stn,tm=Q$tq), roll="nearest", on=.(stn,tm)]
Q[,`:=`(wd=wv$wd,ws=wv$ws,wind_gap_min=abs(wv$wtm-Q$tq))]
Q[,fx:=sin(wd*pi/180)]; Q[,fy:=cos(wd*pi/180)]; Q[,q:=.I]
rm(rwis,windall,wv,U); gc(verbose=FALSE)
cat("Requests:",nrow(Q)," wind direction missing:",sum(is.na(Q$wd))," nearest wind observation >30 min away:",sum(Q$wind_gap_min>30,na.rm=TRUE),"\n")

ts_list<-sort(unique(Q$t5)); cat("Unique times:",length(ts_list),"\n")
CKROOT<-Sys.getenv("RICE_CKPT", file.path(MART,"_ckpt"))
CK<-file.path(CKROOT, if(TEST) "past_test" else "past")
dir.create(CK,showWarnings=FALSE,recursive=TRUE)
CH<-split(ts_list, ceiling(seq_along(ts_list)/250))

work<-function(ci){
  suppressPackageStartupMessages({library(data.table);library(terra)}); terraOptions(memfrac=0.1)
  out<-file.path(CK,sprintf("ch%04d.rds",ci)); if(file.exists(out)) return(0L)
  off<-CJ(dy=-Rreg:Rreg,dx=-Rreg:Rreg)[sqrt(dy^2+dx^2)<=Rreg]
  tspath<-function(t,stem){o<-as.POSIXct(t*60,origin="1970-01-01",tz="UTC")
    file.path(AGB,format(o,"%Y/%m/%d"),paste0("sfc_grid_",stem,"_",format(o,"%Y%m%d%H%M"),".nc"))}
  samp<-function(m,y,x){ ok<-!is.na(y)&!is.na(x)&y>=1&x>=1&y<=nrow(m)&x<=ncol(m)
    v<-rep(NA_real_,length(y)); v[ok]<-m[cbind(y[ok],x[ok])]; v }
  R<-list(); nmiss<-0L
  for(t in CH[[ci]]){ qq<-Q[t5==t]
    for(vn in c("ta","td")){ f<-tspath(t,vn)
      if(!file.exists(f)){ nmiss<-nmiss+1L; next }
      m<-as.matrix(rast(f),wide=TRUE); m[m< -9000]<-NA
      y<-qq$yii; x<-qq$xii
      yy<-pmin(pmax(outer(y,off$dy,"+"),1L),nrow(m)); xx<-pmin(pmax(outer(x,off$dx,"+"),1L),ncol(m))
      vr<-matrix(m[cbind(as.vector(yy),as.vector(xx))],nrow=nrow(qq))
      d<-data.table(q=qq$q,var=vn,vc=samp(m,y,x),vreg=rowMeans(vr,na.rm=TRUE))
      for(k in seq_along(Dkm)){
        d[[paste0("up",Dkm[k])]]<-samp(m,round(y+Dc[k]*qq$fy),round(x+Dc[k]*qq$fx))
        d[[paste0("dn",Dkm[k])]]<-samp(m,round(y-Dc[k]*qq$fy),round(x-Dc[k]*qq$fx)) }
      R[[length(R)+1L]]<-d } }
  saveRDS(list(R=rbindlist(R),nmiss=nmiss),out); 1L }

t0<-Sys.time()
cl<-makeCluster(NW); clusterExport(cl,c("CH","CK","Q","Rreg","Dkm","Dc","AGB"))
todo<-which(!file.exists(file.path(CK,sprintf("ch%04d.rds",seq_along(CH))))); done<-0L
for(b in split(todo, ceiling(seq_along(todo)/NW))){
  parLapply(cl,b,work); done<-done+length(b)
  message(sprintf("  chunk %d/%d  %.1f min", done, length(todo), as.numeric(difftime(Sys.time(),t0,units="mins")))) }
stopCluster(cl)

L<-lapply(file.path(CK,sprintf("ch%04d.rds",seq_along(CH))),readRDS)
R<-rbindlist(lapply(L,`[[`,"R")); nmiss<-sum(sapply(L,`[[`,"nmiss"))
w<-dcast(R, q~var, value.var=c("vc","vreg",paste0("up",Dkm),paste0("dn",Dkm)))
setnames(w, names(w), sub("^(vc|vreg|up5|up10|up20|dn5|dn10|dn20)_(ta|td)$","\\2_\\1",names(w)))
O<-merge(Q[,.(q,row,stn,label,tm_min,lag,wd,ws,wind_gap_min)], w, by="q", all.x=TRUE)
setorder(O,lag,row); O[,q:=NULL]
stopifnot(nrow(O)==nrow(Q))
cat("Missing grid files:",nmiss,"\n")
print(O[,.(n=.N, missing_ta_vc=sum(is.na(ta_vc)), missing_td_up10=sum(is.na(td_up10))),by=lag])
outf<-file.path(MART, if(TEST) "grid_past_scale_raw_TEST.rds" else "grid_past_scale_raw.rds")
saveRDS(O,outf)
cat("Saved:",outf," elapsed:",round(difftime(Sys.time(),t0,units="mins"),1),"min\n")
