# Draws Fig. 4: Kendall tau of rolling variance and lag-1 autocorrelation, sensitivity across settings, and AUC with CIs for both networks.
# Inputs: 02_MART/intl_ews.rds, ews_sensitivity_9combo.csv, table2_auc_ci.csv.
# Output: 05_REPORT/FIGURES/fig4_csd.png.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
MART<-file.path(ROOT, "02_MART")
FIG<-file.path(ROOT, "05_REPORT/FIGURES")
source(code_file("_figure_style.R"))

combo<-fread(file.path(MART,"ews_sensitivity_9combo.csv"),encoding="UTF-8")
stopifnot(nrow(combo)==9L,!any(is.na(combo$var_auc)),!any(is.na(combo$ar1_auc)))

ews<-readRDS(file.path(MART,"intl_ews.rds"))
ews<-ews[type%in%c("onset","control")]
ci<-fread(file.path(MART,"table2_auc_ci.csv"),encoding="UTF-8")

COL_ON<-"#c0392b"; COL_CT<-"grey45"; COL_V<-"#2980b9"; COL_A<-"#e67e22"

figure_device(file.path(FIG,"fig4_csd.png"),width=2200,height=1700,res=200,pointsize=13); panel_letters(letters[1:4])
op<-par(mfrow=c(2,2),mar=c(4.4,4.6,2.8,1),cex.lab=1.10,cex.axis=1.0,cex.main=1.15)

dens_ref<-function(x){x<-x[is.finite(x)]; bw<-bw.nrd0(x)
  d<-density(c(x,-2-x,2-x),bw=bw,from=-1,to=1,n=512); list(x=c(-1,d$x,1),y=c(0,3*d$y,0))}
dens_panel<-function(v,main,xlab){
  don<-dens_ref(ews[type=="onset"][[v]])
  dct<-dens_ref(ews[type=="control"][[v]])
  plot(NA,xlim=c(-1,1),ylim=c(0,max(don$y,dct$y)*1.35),main=main,
       xlab=xlab,ylab="Density")
  polygon(dct,col=adjustcolor(COL_CT,.30),border=COL_CT,lwd=2)
  polygon(don,col=adjustcolor(COL_ON,.25),border=COL_ON,lwd=2)
  segments(0,0,0,max(don$y,dct$y)*1.08,lty=3,col="grey30")
  mon<-median(ews[type=="onset"][[v]],na.rm=TRUE)
  mct<-median(ews[type=="control"][[v]],na.rm=TRUE)
  segments(mon,0,mon,par("usr")[4]*.06,col=COL_ON,lwd=3)
  segments(mct,0,mct,par("usr")[4]*.06,col=COL_CT,lwd=3)
  legend("top",c(sprintf("Pre-onset (n = %s)",format(ews[type=="onset",sum(is.finite(get(v)))],big.mark=",")),
                     sprintf("Control (n = %s)",format(ews[type=="control",sum(is.finite(get(v)))],big.mark=","))),
         fill=adjustcolor(c(COL_ON,COL_CT),.3),border=c(COL_ON,COL_CT),bty="n",cex=1.0,ncol=2)
}
dens_panel("var_tau","",expression("Kendall "*tau[K]*" (rolling variance)"))
dens_panel("ar1_tau","",expression("Kendall "*tau[K]*" (AR1)"))

plot(NA,xlim=c(0.6,9.4),ylim=c(0.38,0.66),xaxt="n",xlab="Window W (h)  \u00d7  detrending bandwidth (min)",
     ylab="Discrimination AUC",main="")
abline(h=0.5,lty=2,lwd=2,col="grey30")
abline(v=c(3.5,6.5),lty=3,col="grey80")
points(1:9,combo$var_auc,pch=19,cex=1.5,col=COL_V,type="b",lwd=2)
points(1:9,combo$ar1_auc,pch=17,cex=1.5,col=COL_A,type="b",lwd=2,lty=2)
axis(1,1:9,rep("",9))
text(1:9,par("usr")[3]-0.012,paste0(combo$W_h,"\u00d7",combo$bw),srt=45,adj=1,xpd=TRUE,cex=.95)
text(0.62,0.496,"AUC = 0.5",adj=c(0,1),col="grey30")
legend("topright",c("Rolling variance","AR1"),col=c(COL_V,COL_A),
       pch=c(19,17),lwd=2,lty=c(1,2),bty="n",cex=1.05)

lab<-c("Variance","AR1","Variance","AR1")
ord<-data.table(network=c("primary","primary","independent","independent"),
                indicator=c("rolling variance tau","lag-1 autocorrelation tau","rolling variance tau","lag-1 autocorrelation tau"))
ci<-ci[ord,on=c("network","indicator")]
stopifnot(nrow(ci)==4L,!any(is.na(ci$AUC)))
xx<-c(1,2,3.4,4.4); colp<-c(COL_V,COL_A,COL_V,COL_A)
plot(NA,xlim=c(0.4,5.0),ylim=c(0.42,0.58),xaxt="n",xlab="",
     ylab="Discrimination AUC (95% CI)",main="")
abline(h=0.5,lty=2,lwd=2,col="grey30")
arrows(xx,ci$ci_low,xx,ci$ci_high,angle=90,code=3,length=.06,lwd=2.2,col=colp)
points(xx,ci$AUC,pch=c(19,17,19,17),cex=1.9,col=colp)
text(xx,ci$AUC,sprintf("%.3f",ci$AUC),pos=4,offset=.6,cex=1.0)
axis(1,xx,lab,cex.axis=.95)
abline(v=2.7,lty=3,col="grey70")
mtext("Primary (121 stations)",side=1,line=2.6,at=1.5,cex=.85,font=2)
mtext("Independent (127 stations)",side=1,line=2.6,at=3.9,cex=.85,font=2)
par(op);dev.off()
cat("Saved: 05_REPORT/FIGURES/fig4_csd.png\n")
