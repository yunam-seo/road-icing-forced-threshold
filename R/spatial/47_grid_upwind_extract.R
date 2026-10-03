# Extracts upwind and downwind air temperature and dew point at 5, 10 and 20 km and the 10 km ring mean for each event and control point.
# Inputs: 02_MART/grid_full_points.rds; 500 m NetCDF fields.
# Output: 02_MART/grid_upwind_raw.rds.
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
gvars<-c(ta="ta", td="td"); Rreg<-20L; CELLKM<-0.5
Dkm<-c(5L,10L,20L); Dc<-as.integer(Dkm/CELLKM)
ts_range<-c(-40,40); ta_range<-c(-40,40); hm_range<-c(0,100); spike_max<-10

pts<-readRDS(file.path(MART,"grid_full_points.rds"))
on<-pts[label==1L]; non<-on[,.(non=.N),by=stn]
ct<-merge(pts[label==0L], non, by="stn")[, .SD[sample(.N, min(.N, non[1]))], by=stn]; ct[,non:=NULL]
P<-rbind(on, ct, fill=TRUE); setkey(P,stn)
args<-commandArgs(trailingOnly=TRUE); if(length(args)>=1L && args[1]=="test") P<-P[sample(.N,40)]
cat("Cases",sum(P$label),"  controls",sum(P$label==0L)," total",nrow(P),"\n")

qc_station<-function(d){ setorder(d,tm)
  d[(ts<ts_range[1])|(ts>ts_range[2]), ts:=NA_real_]; d[(ta<ta_range[1])|(ta>ta_range[2]), ta:=NA_real_]
  d[(rh<hm_range[1])|(rh>hm_range[2]), rh:=NA_real_]
  d[,dp:=ts-data.table::shift(ts)]; d[,dn:=data.table::shift(ts,type="lead")-ts]
  d[!is.na(dp)&!is.na(dn)&abs(dp)>spike_max&abs(dn)>spike_max&sign(dp)!=sign(dn), ts:=NA_real_]
  d[,c("dp","dn"):=NULL]; d }
fs<-list.files(DDIR,pattern="csv$",full.names=TRUE); sof<-as.integer(sub(".*rwis_([0-9]+)_.*","\\1",basename(fs)))
rwis<-list()
for(sid in unique(P$stn)){ i<-which(sof==sid); if(!length(i))next
  d<-fread(fs[i],select=c("datetime","TR1","TA","HM","WS","WD"),na.strings=c("-999",""))
  d[,tm:=as.integer(as.POSIXct(strptime(sprintf("%012.0f",as.numeric(datetime)),"%Y%m%d%H%M",tz="UTC")))%/%60L]
  d<-d[!is.na(tm), .(tm,ts=TR1,ta=TA,rh=HM,ws=WS,wd=WD)]; rwis[[as.character(sid)]]<-qc_station(d); setkey(rwis[[as.character(sid)]],tm) }
P[,tm5:=(tm_min%/%5L)*5L]
wj<-rbindlist(lapply(seq_len(nrow(P)),function(i){w<-rwis[[as.character(P$stn[i])]]
  if(is.null(w))return(data.table(wd=NA_real_,ws=NA_real_)); r<-w[.(P$tm_min[i]),roll="nearest",on="tm"]; data.table(wd=r$wd,ws=r$ws)}))
P[,`:=`(wd0=wj$wd, ws0=wj$ws)]
P[,fx:=sin(wd0*pi/180)]; P[,fy:=cos(wd0*pi/180)]

off<-CJ(dy=-Rreg:Rreg,dx=-Rreg:Rreg)[sqrt(dy^2+dx^2)<=Rreg]
ts_list<-sort(unique(P$tm5))
cat("Unique times:",length(ts_list)," × 2 variables\n")
tspath<-function(t,stem){o<-as.POSIXct(t*60,origin="1970-01-01",tz="UTC")
  file.path(AGB,format(o,"%Y/%m/%d"),paste0("sfc_grid_",stem,"_",format(o,"%Y%m%d%H%M"),".nc"))}
NR<-nrow(P); B<-list()
mk<-function(nm) B[[nm]]<<-rep(NA_real_,NR)
for(vn in gvars) for(s in c("vc","vreg", paste0("up",Dkm), paste0("dn",Dkm))) mk(paste(vn,s,sep="_"))
samp<-function(m,y,x){ if(y<1||x<1||y>nrow(m)||x>ncol(m)) NA_real_ else m[y,x] }

t0<-Sys.time(); nmiss<-0L
for(ti in seq_along(ts_list)){ t<-ts_list[ti]; rr<-which(P$tm5==t); if(!length(rr))next
  for(vn in names(gvars)){ f<-tspath(t,gvars[[vn]]); if(!file.exists(f)){nmiss<-nmiss+1L; next}
    m<-as.matrix(rast(f),wide=TRUE); m[m< -9000]<-NA
    for(i in rr){ x<-P$xii[i]; y<-P$yii[i]; fx<-P$fx[i]; fy<-P$fy[i]
      B[[paste0(vn,"_vc")]][i]<-samp(m,y,x)
      B[[paste0(vn,"_vreg")]][i]<-mean(m[cbind(pmin(pmax(y+off$dy,1),nrow(m)),pmin(pmax(x+off$dx,1),ncol(m)))],na.rm=TRUE)
      if(!is.na(fx)) for(k in seq_along(Dkm)){ dc<-Dc[k]
        B[[paste0(vn,"_up",Dkm[k])]][i]<-samp(m, round(y+dc*fy), round(x+dc*fx))
        B[[paste0(vn,"_dn",Dkm[k])]][i]<-samp(m, round(y-dc*fy), round(x-dc*fx)) }
    } }
  if(ti%%500L==0L) message("  ",ti,"/",length(ts_list)," (",round(difftime(Sys.time(),t0,units="mins"),1),"min, missing",nmiss,")") }
F<-as.data.table(B); F[,`:=`(row=seq_len(NR), stn=P$stn, label=P$label, tm_min=P$tm_min, rts=NA)]
rv<-rbindlist(lapply(seq_len(NR),function(i){w<-rwis[[as.character(P$stn[i])]]; if(is.null(w))return(data.table(rts=NA_real_,rta=NA,rrh=NA,rws=NA))
  r<-w[.(P$tm_min[i]),roll="nearest",on="tm"]; data.table(rts=r$ts,rta=r$ta,rrh=r$rh,rws=r$ws)}))
F[,`:=`(rts=rv$rts,rta=rv$rta,rrh=rv$rrh,rws=rv$rws)]
cat("Missing grid files:",nmiss,"\n")
saveRDS(F, file.path(MART,"grid_upwind_raw.rds"))
cat("Saved: 02_MART/grid_upwind_raw.rds  elapsed:",round(difftime(Sys.time(),t0,units="mins"),1),"min\n")
