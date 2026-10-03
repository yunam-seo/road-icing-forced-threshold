# Measures spatial correlation of the objectively analysed air temperature and dew point as a function of distance from road points.
# Inputs: 02_MART/road_grid500_match.csv, intl_event_table.rds; 500 m NetCDF fields.
# Outputs: 02_MART/spatial_corr_scale.csv, spatial_corr_scale.rds; 05_REPORT/FIGURES/intl_spatial_corr.png.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages({library(data.table);library(terra)})
set.seed(2026)
MART<-file.path(ROOT, "02_MART"); FIG<-file.path(ROOT, "05_REPORT/FIGURES")
AGB<-grid_fields()
mt<-fread(file.path(MART,"road_grid500_match.csv"))
ev<-as.data.table(readRDS(file.path(MART,"intl_event_table.rds")))
ev[,tm5:=(onset_min%/%5L)*5L]
ts_samp<-sort(sample(unique(ev$tm5), 250))
cat("Sampled times:",length(ts_samp)," road points:",nrow(mt),"\n")
tspath<-function(t,stem){o<-as.POSIXct(t*60,origin="1970-01-01",tz="UTC")
  file.path(AGB,format(o,"%Y/%m/%d"),paste0("sfc_grid_",stem,"_",format(o,"%Y%m%d%H%M"),".nc"))}
dirs<-list(c(1,0),c(-1,0),c(0,1),c(0,-1),c(1,1),c(1,-1),c(-1,1),c(-1,-1))
dist_cells<-1:12; CELLKM<-0.5
nmiss<-0L; pairs<-list()
for(t in ts_samp){ for(vn in c("ta","td")){ f<-tspath(t,vn)
  if(!file.exists(f)){nmiss<-nmiss+1L; next}
  m<-as.matrix(rast(f),wide=TRUE); m[m< -9000]<-NA
  for(j in seq_len(nrow(mt))){ x<-mt$xii[j]; y<-mt$yii[j]; c0<-m[y,x]; if(is.na(c0))next
    for(dd in dist_cells){ vals<-sapply(dirs,function(D){yy<-y+D[1]*dd; xx<-x+D[2]*dd
      if(yy<1||xx<1||yy>nrow(m)||xx>ncol(m)) NA else m[yy,xx]})
      pairs[[length(pairs)+1L]]<-data.table(var=vn,d=dd,c0=c0,v=mean(vals,na.rm=TRUE)) } } } }
P<-rbindlist(pairs)[is.finite(c0)&is.finite(v)]
cat("Missing grids:",nmiss," pairs:",nrow(P),"\n")
cor_d<-P[,.(r=cor(c0,v), km=d[1]*CELLKM, n=.N),by=.(var,d)][order(var,d)]
cat("\n=== Spatial correlation r(d) by distance ===\n"); print(cor_d[,.(var,km,r=round(r,3),n)])
declen<-function(km,r,thr){ i<-which(r<thr)[1]; if(is.na(i)||i==1) return(NA)
  km[i-1]+(km[i]-km[i-1])*(r[i-1]-thr)/(r[i-1]-r[i]) }
for(vv in c("ta","td")){ s<-cor_d[var==vv]
  cat(sprintf("%s: r(3km)=%.3f · L(r=0.5)=%.2fkm · L(r=1/e)=%.2fkm\n",
    vv, s[km==3]$r, declen(s$km,s$r,0.5), declen(s$km,s$r,exp(-1)))) }
saveRDS(cor_d,file.path(MART,"spatial_corr_scale.rds"))
fwrite(cor_d,file.path(MART,"spatial_corr_scale.csv"),bom=TRUE)

png(file.path(FIG,"intl_spatial_corr.png"),width=1100,height=620,res=140)
par(mar=c(4.6,4.6,3,1)); s_ta<-cor_d[var=="ta"]; s_td<-cor_d[var=="td"]
plot(s_ta$km,s_ta$r,type="b",pch=19,col="#c0392b",lwd=2.2,ylim=c(0,1),xlab="Distance from road point (km)",
  ylab="Spatial correlation r", main="Spatial correlation of the objectively analysed fields (basis for neighbourhood radius R)")
lines(s_td$km,s_td$r,type="b",pch=17,col="#2e86c1",lwd=2.2)
abline(v=3,lty=2,col="#7f7f7f"); abline(h=c(0.5,exp(-1)),lty=3,col="#aaa")
text(3,0.05,"R=3km (used)",pos=4,cex=0.85,col="#555")
legend("topright",c("air temperature ta","dew point td"),col=c("#c0392b","#2e86c1"),pch=c(19,17),lwd=2.2,bty="n")
dev.off(); cat("\nFigure: 05_REPORT/FIGURES/intl_spatial_corr.png\n")
