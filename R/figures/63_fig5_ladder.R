# Draws Fig. 5: ROC-AUC, PR-AUC, Brier skill and trigger metrics of the information ladder (Models 0-2).
# Inputs: 02_MART/model_ladder_metrics.csv, intl_ew_features.rds.
# Outputs: 05_REPORT/FIGURES/fig5_ladder.png; 02_MART/fig5_base_rate.csv.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
MART<-file.path(ROOT, "02_MART")
FIG<-file.path(ROOT, "05_REPORT/FIGURES")
source(code_file("_figure_style.R"))
m<-fread(file.path(MART,"model_ladder_metrics.csv"))
stopifnot(nrow(m)==3)
labs<-c("Model 0\nTs threshold","Model 1\ncurrent {Ts, RH}","Model 2\n+ recent change\n+ cumulative")
cols<-c("grey55","#2980b9","#c0392b")
EW<-readRDS(file.path(MART,"intl_ew_features.rds")); allf<-c("cold_dose","dep_dose","ft","ts","rh","dts_1h")
EW<-EW[complete.cases(EW[,..allf])]; BASE<-mean(!is.na(EW$tno)&EW$tno<=180L); rm(EW)
cat(sprintf("Climatological base rate (T+3 h positive rate) %.5f
",BASE))
n_eval<-nrow(readRDS(file.path(MART,"intl_ew_features.rds")))
BR<-data.table(item="T+3 h positive rate (climatological base rate)",value=BASE,n_evaluation_points=n_eval)
fwrite(BR,file.path(MART,"fig5_base_rate.csv"),bom=TRUE)

figure_device(file.path(FIG,"fig5_ladder.png"),width=2200,height=1700,res=200,pointsize=13); panel_letters(letters[1:4])
op<-par(mfrow=c(2,2),mar=c(5.2,4.6,2.8,1),cex.lab=1.10,cex.axis=1.0,cex.main=1.15)

ladder<-function(v,main,ylab,ylim,digits=3,hline=NULL,htxt=NULL){
  plot(NA,xlim=c(.6,3.4),ylim=ylim,xaxt="n",xlab="",ylab=ylab,main=main)
  if(!is.null(hline)){abline(h=hline,lty=2,lwd=2,col="grey40")
    text(3.35,hline,htxt,adj=c(1,-0.5),col="grey40",cex=.95)}
  lines(1:3,v,lwd=2,col="grey60")
  points(1:3,v,pch=19,cex=2.1,col=cols)
  text(1:3,v,formatC(v,format="f",digits=digits),pos=3,offset=.8,cex=1.1)
  axis(1,1:3,labs,padj=.5,cex.axis=.95)
}
ladder(m$ROC_AUC,"","ROC-AUC",c(0.68,0.95),3,0.5,NULL)
ladder(m$PR_AUC,"","PR-AUC (average precision)",
       c(0,0.105),3,BASE,sprintf("climatological base rate %.4f",BASE))

plot(NA,xlim=c(.6,3.4),ylim=c(0,0.05),xaxt="n",xlab="",ylab="Brier skill score",
     main="")
rect(1:3-.28,0,1:3+.28,c(NA,m$BSS[2],m$BSS[3]),col=adjustcolor(cols,.75),border=cols,lwd=2)
text(2:3,m$BSS[2:3],formatC(m$BSS[2:3],format="f",digits=3),pos=3,offset=.5,cex=1.1)
text(1,0.004,"n/a\n(no probability)",cex=1.0,col="grey40")
axis(1,1:3,labs,padj=.5,cex.axis=.95)

plot(NA,xlim=c(.55,3.45),ylim=c(0,max(m$FA_per_day)*1.42),xaxt="n",xlab="",
     ylab="False alarms per station-day",main="")
rect(1:3-.28,0,1:3+.28,m$FA_per_day,col=adjustcolor(cols,.75),border=cols,lwd=2)
text(1:3,m$FA_per_day,sprintf("%.2f",m$FA_per_day),pos=3,offset=.4,cex=1.1)
pod_lab<-sprintf("POD %.2f",m$POD)
ann<-sprintf("%s\nCSI %.3f\nlead %d min",pod_lab,m$CSI_point,m$lead_med_min)
text(1:3,m$FA_per_day+max(m$FA_per_day)*0.10,ann,pos=3,cex=.95)
axis(1,1:3,labs,padj=.5,cex.axis=.95)
par(op);dev.off()
cat("Saved: 05_REPORT/FIGURES/fig5_ladder.png\n")
