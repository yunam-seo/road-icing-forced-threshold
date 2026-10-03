# Tests how much prior surface state explains onset conditions and analyses the principal components of onset conditions (Hartigan's dip test).
# Inputs: 02_MART/intl_event_table.rds, intl_lagged_table.rds, intl_cumstress.rds, intl_physical_diag.csv, intl_earlywarning_eval.rds.
# Outputs: 02_MART/intl_continuum.csv; 05_REPORT/FIGURES/intl_equifinality.png.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages({library(data.table);library(lightgbm)})
set.seed(2026)
MART<-file.path(ROOT, "02_MART"); FIG<-file.path(ROOT, "05_REPORT/FIGURES")

ev<-readRDS(file.path(MART,"intl_event_table.rds")); ev[,eid:=paste0(stn,"_",onset_min)]
ph<-fread(file.path(MART,"intl_physical_diag.csv"))
cs<-readRDS(file.path(MART,"intl_cumstress.rds"))[type=="onset",.(eid,cold_dose)]
lg<-readRDS(file.path(MART,"intl_lagged_table.rds"))
D<-merge(ph,cs,by="eid"); D<-merge(D,ev[,.(eid,pre_state,ts0,rh0,m_liq)],by="eid")

cat("=== (1) Equifinality: onset conditions by prior state (convergence?) ===\n")
print(D[,.(n=.N, ts0=round(median(ts0,na.rm=T),2), rh0=round(median(rh0,na.rm=T),1),
  Mliq=round(median(m_liq,na.rm=T),2), amplitude=round(median(ts_amp,na.rm=T),1)),
  by=pre_state][order(-n)])
eta2<-function(yv,g){g<-factor(g);gm<-tapply(yv,g,mean,na.rm=TRUE)
  ss_b<-sum(table(g)*(gm-mean(yv,na.rm=TRUE))^2,na.rm=TRUE)
  ss_t<-sum((yv-mean(yv,na.rm=TRUE))^2,na.rm=TRUE);ss_b/ss_t}
for(v in c("ts0","rh0","m_liq","ts_amp"))
  cat(sprintf("  η²(prior state→%s)=%.3f\n",v,eta2(D[[v]],D$pre_state)))

w<-dcast(melt(lg,id.vars=c("eid","lag"),measure.vars=c("ts","m_dep","rh")),
         eid~variable+lag,value.var="value")
w<-merge(w,ev[,.(eid,pre_state)],by="eid");w<-w[complete.cases(w)]
w[,y:=as.integer(factor(pre_state,levels=c("DRY","WET","SNOW")))-1L]
Xc<-as.matrix(w[,setdiff(names(w),c("eid","pre_state","y")),with=FALSE])
tr<-sample(nrow(w),round(.7*nrow(w)))
m<-lgb.train(list(objective="multiclass",num_class=3,learning_rate=0.05,num_leaves=15,
  min_data_in_leaf=20,verbosity=-1),lgb.Dataset(Xc[tr,],label=w$y[tr]),nrounds=200)
pp<-predict(m,Xc[-tr,]);if(is.null(dim(pp)))pp<-matrix(pp,ncol=3,byrow=TRUE)
pred<-max.col(pp)-1L;yte<-w$y[-tr]
acc<-mean(pred==yte);maj<-max(table(yte))/length(yte)
bacc<-mean(sapply(0:2,function(k)if(sum(yte==k))mean(pred[yte==k]==k)else NA),na.rm=TRUE)
cat(sprintf("\n  Trajectory→prior state classification accuracy %.3f (majority baseline %.3f, balanced accuracy %.3f)\n",acc,maj,bacc))

cat("\n=== (2) Continuous spectrum (PCA) ===\n")
fe<-c("ts_amp","onset_hour","precip_frac","noct_ts_ta","cold_dose","rh0","m_liq")
P<-D[complete.cases(D[,..fe])]; Z<-scale(P[,..fe]); Z[,"cold_dose"]<-scale(log1p(P$cold_dose))
pc<-prcomp(Z)
cat("  PC1 explained variance:",round(pc$sdev[1]^2/sum(pc$sdev^2),3),
    " PC2:",round(pc$sdev[2]^2/sum(pc$sdev^2),3),"\n")
cat("  PC1 loadings (sign = direction):\n");print(round(sort(pc$rotation[,1]),2))
P[,PC1:=pc$x[,1]];P[,PC2:=pc$x[,2]]
dip_p<-NA
if(requireNamespace("diptest",quietly=TRUE)){
  dt<-diptest::dip.test(P$PC1);dip_p<-dt$p.value
  cat(sprintf("  Hartigan dip test p=%.3f (>0.05 = unimodal = continuum, not discrete clusters)\n",dip_p))
} else cat("  (diptest not installed — check unimodality on the histogram)\n")

ee<-readRDS(file.path(MART,"intl_earlywarning_eval.rds"))
ee<-ee[!is.na(tno)&tno>0&tno<=180]
ee[,oid:=paste0(stn,"_",tm_min+tno)]
wn<-ee[,.(maxsc=max(score)),by=oid]
P<-merge(P,wn,by.x="eid",by.y="oid",all.x=TRUE)
P[,pc1bin:=cut(PC1,quantile(PC1,0:4/4),include.lowest=TRUE,labels=c("Q1 slow-wet","Q2","Q3","Q4 fast-dry"))]
cat("\n=== (3) Characteristics and warning difficulty by PC1 quartile ===\n")
print(P[!is.na(pc1bin),.(n=.N,amplitude=round(median(ts_amp),1),precipitation=round(mean(precip_any),2),
  median_onset_hour=round(median(onset_hour),0),dose=round(median(cold_dose)),
  detection_score=round(median(maxsc,na.rm=T),3)),by=pc1bin][order(pc1bin)])

png(file.path(FIG,"intl_equifinality.png"),width=2300,height=760,res=200)
op<-par(mfrow=c(1,3),mar=c(4.4,4.6,2.6,1),cex.lab=1.25,cex.axis=1.1,cex.main=1.2)
cols3<-c(DRY="#e67e22",WET="#2980b9",SNOW="#7f8c8d")
plot(NA,xlim=c(-8,4),ylim=c(0,0.5),xlab="Surface temperature at onset (\u00b0C)",ylab="Density",
  main="(a) Equifinality: onset condition by prior state")
for(s in names(cols3)){d<-density(D[pre_state==s]$ts0,na.rm=TRUE);lines(d,col=cols3[s],lwd=2)}
abline(v=0,lty=3);legend("topleft",names(cols3),col=cols3,lwd=2.5,bty="n",cex=1.0)
hist(P$PC1,breaks=40,col="#b8d8e8",border="white",xlab="PC1 (slow\u00b7wet\u00b7dawn \u2192 fast\u00b7dry)",
  main=sprintf("(b) Continuous spectrum%s",if(!is.na(dip_p))sprintf(" (dip p=%.2f)",dip_p)else""))
P[,p1d:=cut(PC1,quantile(PC1,0:10/10),include.lowest=TRUE,labels=FALSE)]
g<-P[,.(pc=median(PC1),precip=mean(precip_any,na.rm=T),amp=median(ts_amp,na.rm=T),
        hour=median(onset_hour,na.rm=T)),by=p1d][order(p1d)]
plot(g$pc,g$precip,type="b",pch=19,col="#2980b9",lwd=2,ylim=c(0,1.32),yaxt="n",
  xlab="PC1 (slow\u00b7wet \u2192 fast\u00b7dry)",ylab="Fraction / normalised",main="(c) Physical gradient along PC1")
axis(2,at=seq(0,1,0.2))
lines(g$pc,(g$amp-min(g$amp))/(max(g$amp)-min(g$amp)),type="b",pch=17,col="#c0392b",lwd=2)
lines(g$pc,g$hour/24,type="b",pch=15,col="#27ae60",lwd=2)
legend("topleft",inset=c(0.01,0.01),c("Precipitation co-occurrence","Surface temp. amplitude (norm.)","Onset hour / 24"),
  col=c("#2980b9","#c0392b","#27ae60"),pch=c(19,17,15),lwd=2,bty="n",cex=.9)
par(op);dev.off()
fwrite(P[,.(eid,PC1,PC2,ts_amp,onset_hour,precip_frac,cold_dose,pre_state,maxsc)],
  file.path(MART,"intl_continuum.csv"))
cat("\nSaved: 05_REPORT/FIGURES/intl_equifinality.png, intl_continuum.csv\n")
