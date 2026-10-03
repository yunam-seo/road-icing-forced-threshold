# Draws Fig. 6: AUC by lead with and without surrounding fields, false alarms at equal POD, and event-control AUC under three holdout schemes.
# Inputs: 02_MART/pattern_by_lead_auc.csv, trigger_unified.csv, grid_holdout_auc.csv.
# Output: 05_REPORT/FIGURES/fig6_falsealarm.png.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
MART<-file.path(ROOT, "02_MART")
FIG<-file.path(ROOT, "05_REPORT/FIGURES")
source(code_file("_figure_style.R"))
C0<-"grey50"; C1<-"#c0392b"; C2<-"#2980b9"; C3<-"#16a085"

LEAD<-fread(file.path(MART,"pattern_by_lead_auc.csv"),encoding="UTF-8")[order(lead_min)]
TR<-fread(file.path(MART,"trigger_unified.csv"),encoding="UTF-8")[group=="Fig. 6b (with and without surrounding fields)"]
HO<-fread(file.path(MART,"grid_holdout_auc.csv"),encoding="UTF-8")[model %in% c("RWIS","concurrent surrounding fields","surrounding fields 1 h earlier","surrounding fields 3 h earlier")]
stopifnot(nrow(LEAD)==5L, nrow(TR)==4L, nrow(HO)==12L)
tr<-function(k){ r<-TR[configuration==k]; stopifnot(nrow(r)==1L); r }

figure_device(file.path(FIG,"fig6_falsealarm.png"),width=2600,height=850,res=200,pointsize=19); panel_letters(letters[1:3])
op<-par(mfrow=c(1,3),mar=c(4.6,4.6,2.8,1),cex.lab=1.10,cex.axis=1.0,cex.main=1.15)

lead<-LEAD$lead_min/60; rwis<-LEAD$auc_rwis_only; spat<-LEAD$auc_with_surrounding_fields
yl<-range(c(rwis,spat)); yl<-c(yl[1]-0.02, yl[2]+0.025)
plot(NA,xlim=c(0.8,13.5),ylim=yl,log="x",xaxt="n",
     xlab="Forecast horizon (h, log scale)",ylab="Station-held-out CV AUC",
     main="")
axis(1,lead,labels=FALSE)
text(lead,par("usr")[3]-0.006,paste0("T+",lead),adj=c(0.5,1),xpd=TRUE)
lines(lead,rwis,type="b",pch=19,cex=1.5,lwd=2.2,col=C0)
lines(lead,spat,type="b",pch=17,cex=1.5,lwd=2.2,col=C1)
legend("bottomleft",c("+ surrounding fields","RWIS only"),col=c(C1,C0),
       pch=c(17,19),lwd=2.2,pt.cex=1.5,bty="n",cex=1.05)
g<-spat-rwis
text(2.6,yl[2]-0.005,sprintf("\u0394AUC +%.3f to +%.3f",min(g),max(g)),col=C1,cex=1.0)

k<-c("T+6 RWIS","T+6 RWIS + surrounding fields","T+3 RWIS","T+3 RWIS + surrounding fields")
fa<-sapply(k,function(z) tr(z)$false_alarms_per_station_day)
xx<-c(1,1.85,3.3,4.15)
plot(NA,xlim=c(0.4,4.8),ylim=c(0,max(fa)*1.60),xaxt="n",xlab="",
     ylab="False alarms per station-day",main="")
rect(xx-.32,0,xx+.32,fa,col=adjustcolor(c(C0,C1,C0,C1),.75),border=c(C0,C1,C0,C1),lwd=2)
text(xx,fa,sprintf("%.2f",fa),pos=3,cex=1.05)
axis(1,c(1.42,3.72),c("T+6 h","T+3 h"))
legend("topright",c("RWIS only","+ surrounding fields"),fill=adjustcolor(c(C0,C1),.75),border=c(C0,C1),bty="n",cex=.95)
for(j in c(1,3)){ xm<-(xx[j]+xx[j+1])/2; top<-fa[j]*0.94; bot<-fa[j+1]*1.14
  arrows(xm,top,xm,bot,length=.08,lwd=2,col=C1)
  text(xm+0.08,(top+bot)/2,sprintf("\u2212%.0f%%",100*(fa[j]-fa[j+1])/fa[j]),adj=c(0,0.5),col=C1,cex=0.95,font=2) }

sp<-c("stn","day","winter"); cv<-c("Station-\nheld-out","Date-\nheld-out","Winter-\nheld-out")
get<-function(m) sapply(sp,function(s){ v<-HO[holdout==s & model==m]$AUC; stopifnot(length(v)==1L); v })
b0<-get("RWIS"); b1<-get("concurrent surrounding fields"); p1<-get("surrounding fields 1 h earlier"); p3<-get("surrounding fields 3 h earlier")
yl<-range(c(b0,b1,p1,p3)); yl<-c(yl[1]-0.008, yl[2]+0.024)
xx<-1:3
plot(NA,xlim=c(0.35,3.75),ylim=yl,xaxt="n",xlab="",
     ylab="Discrimination AUC",main="")
segments(xx,b0,xx,b1,lwd=2,col="grey70")
points(xx,b0,pch=19,cex=1.7,col=C0)
points(xx,b1,pch=17,cex=1.7,col=C1)
points(xx+0.18,p1,pch=5,cex=1.6,lwd=2,col=C2)
points(xx+0.30,p3,pch=0,cex=1.5,lwd=2,col=C3)
text(xx,b0,sprintf("%.3f",b0),pos=1,cex=1.0)
text(xx,b1,sprintf("%.3f",b1),pos=3,cex=1.0)
axis(1,xx,cv,padj=.5,cex.axis=1.0)
legend("top",c("RWIS only","+ fields, same time","+ fields, 1 h earlier","+ fields, 3 h earlier"),
       col=c(C0,C1,C2,C3),pch=c(19,17,5,0),pt.cex=1.3,bty="n",cex=0.9,ncol=2)
par(op);dev.off()
cat("Saved: 05_REPORT/FIGURES/fig6_falsealarm.png\n")
