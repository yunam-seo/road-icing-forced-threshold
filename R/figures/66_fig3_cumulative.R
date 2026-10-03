# Draws Fig. 3: cold dose, deposition dose, freeze-thaw accumulation, the tau_d sweep and the information ladder.
# Inputs: 02_MART/intl_cumstress.rds, tau_sweep_lagged.csv, model_ladder_full.csv.
# Output: 05_REPORT/FIGURES/fig3_cumulative.png.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
MART<-file.path(ROOT, "02_MART")
FIG<-file.path(ROOT, "05_REPORT/FIGURES")
source(code_file("_figure_style.R"))
cs<-readRDS(file.path(MART,"intl_cumstress.rds"))
cs[,grp:=factor(type,c("control","onset"),c("Control","At onset"))]
sw<-fread(file.path(MART,"tau_sweep_lagged.csv"),encoding="UTF-8")[order(evaluation_time,tau_min)]
stopifnot(setequal(unique(sw$evaluation_time),c("t0-0","t0-30","t0-60")), nrow(sw)>=30L)

figure_device(file.path(FIG,"fig3_cumulative.png"),width=2600,height=1600,res=200,pointsize=19); panel_letters(letters[1:5])
layout(matrix(c(1,2,3,4,5,5),2,3,byrow=TRUE))
par(mar=c(4.2,4.6,2.8,1),cex.lab=1.10,cex.axis=1.0,cex.main=1.15)
cc<-c("grey80","tomato")

boxplot(cold_dose~grp,cs,outline=FALSE,col=cc,xlab="",
        ylab=expression("Cold dose ("*degree*"C"%.%"min)"),
        main="")
boxplot(log10(1+dep_dose)~grp,cs,outline=FALSE,col=cc,xlab="",
        ylab=expression(log[10]*"(1 + deposition dose)"),
        main="")
boxplot(ft~grp,cs,outline=FALSE,col=cc,xlab="",
        ylab="Freeze\u2013thaw accumulation",main="")

yl<-range(sw$AUC); yl<-c(yl[1]-0.008,yl[2]+0.012)
CC<-c("t0-0"="navy","t0-30"="#16a085","t0-60"="#c0392b"); PCH<-c("t0-0"=19,"t0-30"=17,"t0-60"=15)
plot(NA,xlim=range(sw$tau_min),ylim=yl,log="x",xaxt="n",
     xlab=expression("Decay time constant "*tau[d]*" (log scale)"),
     ylab="Onset discrimination AUC",
     main="")
for(k in names(CC)){ g<-sw[evaluation_time==k]
  lines(g$tau_min,g$AUC,type="b",pch=PCH[k],col=CC[k],lwd=2,cex=1.1)
  im<-which.max(g$AUC); points(g$tau_min[im],g$AUC[im],pch=21,bg="gold",cex=1.8,lwd=1.2) }
tk<-c(1,10,60,360,4320); axis(1,at=tk,labels=FALSE)
text(tk,par("usr")[3]-0.004,c("1 min","10 min","1 h","6 h","72 h"),adj=c(0.5,1),xpd=TRUE,cex=0.95)
legend("bottomleft",c("at onset","30 min before","60 min before"),col=CC,pch=PCH,lwd=2,bty="n",cex=0.95)

LAD<-fread(file.path(MART,"model_ladder_full.csv"),encoding="UTF-8")
setkey(LAD,"stage")
pick<-function(k){v<-LAD[k]$AUC; stopifnot(length(v)==1L,!is.na(v)); v}
lad<-data.table(
  lab=c("Ts threshold","Antecedent\nconditions only",
        "Current\n{Ts, RH}","+ recent\nchange","Full\n(+ cumulative)"),
  auc=c(pick("naive_ts"),pick("dose_only"),pick("current"),
        pick("plus_dynamics"),pick("full")))
plot(NA,xlim=c(.6,6.4),ylim=c(0.63,0.95),xaxt="n",xlab="",
     ylab="Station-held-out CV AUC (T+3 h)",main="")
colsE<-c("grey55","#8e44ad","#2980b9","#27ae60","#c0392b")
yc<-lad$auc[3]; yd<-lad$auc[4]; yf<-lad$auc[5]
segments(3,yc,5,yc,lty=3,col="grey50")
lines(1:5,lad$auc,lwd=2,col="grey70")
points(1:5,lad$auc,pch=19,cex=2.0,col=colsE)
text(1:5,lad$auc,sprintf("%.3f",lad$auc),pos=3,offset=.7,cex=1.05)
axis(1,1:5,lad$lab,padj=.5,cex.axis=.92)
arrows(5.30,yc,5.30,yf,length=.06,code=3,lwd=1.8,col="#c0392b")
text(5.30,yc-0.010,sprintf("recent + cumulative\n%+.3f",yf-yc),adj=c(0.5,1),cex=.82,col="#c0392b")
arrows(5.75,yd,5.75,yf,length=.06,code=3,lwd=1.8,col="#27ae60")
text(5.75,yf+0.010,sprintf("cumulative\n%+.3f",yf-yd),adj=c(0.5,0),cex=.82,col="#27ae60")
dev.off()
cat("Saved: 05_REPORT/FIGURES/fig3_cumulative.png\n")
