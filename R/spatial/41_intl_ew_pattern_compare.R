# Compares the operational nowcast with and without surrounding-field features at several horizons (station-held-out CV).
# Inputs: 02_MART/intl_ew_features.rds, grid_ew_scale_raw.rds.
# Outputs: 02_MART/intl_ew_pattern_store.rds; 05_REPORT/FIGURES/intl_ew_pattern.png.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages({library(data.table);library(lightgbm)})
set.seed(2026)
MART<-file.path(ROOT, "02_MART"); FIG<-file.path(ROOT, "05_REPORT/FIGURES")
source(code_file("_grid_features.R"))
base<-c("cold_dose","dep_dose","ft","ts","rh","dts_1h")
pat<-FEAT_SPATIAL_FINAL
Hs<-c(60L,120L,180L,360L,720L); step<-30L; ppd<-24*60/step

D<-as.data.table(readRDS(file.path(MART,"intl_ew_features.rds")))
G<-as.data.table(readRDS(file.path(MART,"grid_ew_scale_raw.rds")))
stopifnot(!anyDuplicated(G[,.(stn,tm_min)]))
setnames(G,c("ta_vc0","td_vc0","ta_vc30","td_vc30","ta_vreg0","td_vreg0"),
         c("ta_vc","td_vc","ta_vc_prev","td_vc_prev","ta_vreg","td_vreg"))
Fp<-build_final_features(G, keys=c("stn","tm_min"))
D<-merge(D, Fp, by=c("stn","tm_min"), all.x=TRUE)
cat("Merged:",nrow(D)," missing features (ta_anom10):",sum(is.na(D$ta_anom10)),
    " missing upwind (dd_up10):",sum(is.na(D$dd_up10))," (rows masked by wind gap >30 min:",sum(Fp$wind_masked),")\n")
D<-D[complete.cases(D[,..base])]; setorder(D,stn,tm_min)
stns<-unique(D$stn); fold<-setNames(sample(rep(1:5,length.out=length(stns))),stns)
D[,fold:=fold[as.character(stn)]]
auc<-function(s,p){o<-is.finite(s);s<-s[o];p<-p[o];np<-as.numeric(sum(p));nn<-as.numeric(sum(!p))
  if(!np||!nn)return(NA);r<-rank(s);(sum(r[p])-np*(np+1)/2)/(np*nn)}
params<-list(objective="binary",learning_rate=0.05,num_leaves=31,min_data_in_leaf=100,
  feature_fraction=0.8,bagging_fraction=0.8,bagging_freq=1,verbosity=-1)
oof<-function(H,fs){fs<-intersect(fs,names(D)); y<-as.integer(!is.na(D$tno)&D$tno<=H)
  pr<-rep(NA_real_,nrow(D))
  for(kf in 1:5){tr<-which(D$fold!=kf);te<-which(D$fold==kf)
    m<-lgb.train(params,lgb.Dataset(as.matrix(D[tr,..fs]),label=y[tr]),nrounds=300)
    p<-predict(m,as.matrix(D[te,..fs]));if(is.matrix(p))p<-p[,1];pr[te]<-p}
  list(y=y,pr=pr)}

cat("\n=== (A) Station-held-out CV AUC at multiple horizons: RWIS vs +surrounding fields ===\n")
ST<-list()
for(H in Hs){ ob<-oof(H,base); op<-oof(H,c(base,pat)); ST[[as.character(H)]]<-list(b=ob,p=op)
  cat(sprintf("  T+%-3.0fh  RWIS %.4f → +surrounding fields %.4f  (Δ%+.4f)\n",
      H/60, auc(ob$pr,ob$y==1L), auc(op$pr,op$y==1L), auc(op$pr,op$y==1L)-auc(ob$pr,ob$y==1L))) }
attr(ST,"meta")<-list(features=pat, base=base, date=as.character(Sys.time()), n=nrow(D),
                      n_wind_masked=sum(Fp$wind_masked))
saveRDS(ST,file.path(MART,"intl_ew_pattern_store.rds"))
cat("Saved: 02_MART/intl_ew_pattern_store.rds\n")

png(file.path(FIG,"intl_ew_pattern.png"),width=1600,height=1000,res=200)
par(mar=c(4.6,4.8,3,1),cex.lab=1.25,cex.axis=1.1,cex.main=1.25)
ab<-sapply(Hs,function(H)auc(ST[[as.character(H)]]$b$pr,ST[[as.character(H)]]$b$y==1L))
ap<-sapply(Hs,function(H)auc(ST[[as.character(H)]]$p$pr,ST[[as.character(H)]]$p$y==1L))
plot(seq_along(Hs),ap,type="b",pch=19,col="#c0392b",lwd=2.5,xaxt="n",ylim=c(0.8,0.98),
  xlab="Forecast lead time",ylab="Station-held-out CV AUC",
  main="Operational early warning: added value of surrounding fields")
lines(seq_along(Hs),ab,type="b",pch=17,col="#2c3e50",lwd=2.5,lty=2)
axis(1,seq_along(Hs),paste0("T+",Hs/60,"h"))
legend("topright",c("RWIS + surrounding fields","RWIS only"),col=c("#c0392b","#2c3e50"),
  pch=c(19,17),lwd=2.5,lty=c(1,2),bty="n")
dev.off(); cat("\nSaved: 05_REPORT/FIGURES/intl_ew_pattern.png (diagnostic figure)\n")
