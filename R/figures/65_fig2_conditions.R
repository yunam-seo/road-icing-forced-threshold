# Draws Fig. 2: onset conditions by prior state, cooling trajectories and the principal-component continuum.
# Inputs: 02_MART/intl_event_table.rds, intl_physical_diag.csv, intl_lagged_table.rds, intl_continuum.csv.
# Outputs: 05_REPORT/FIGURES/fig2_conditions.png; 02_MART/fig2c_trajectory_median.csv.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages({library(data.table);library(diptest)})
MART<-file.path(ROOT, "02_MART")
FIG<-file.path(ROOT, "05_REPORT/FIGURES")
source(code_file("_figure_style.R"))
ev<-readRDS(file.path(MART,"intl_event_table.rds")); ev[,eid:=paste0(stn,"_",onset_min)]
ph<-fread(file.path(MART,"intl_physical_diag.csv"))
D<-merge(ph,ev[,.(eid,pre_state,ts0,rh0,m_liq)],by="eid")
lg<-readRDS(file.path(MART,"intl_lagged_table.rds"))
P<-fread(file.path(MART,"intl_continuum.csv"))

eta2<-function(yv,g){g<-factor(g);gm<-tapply(yv,g,mean,na.rm=TRUE)
  ss_b<-sum(table(g)*(gm-mean(yv,na.rm=TRUE))^2,na.rm=TRUE)
  ss_t<-sum((yv-mean(yv,na.rm=TRUE))^2,na.rm=TRUE);ss_b/ss_t}
evars<-c(ts0="T[s] at onset",rh0="RH at onset",ts_amp="24-h temperature range")
e2<-sapply(names(evars),function(v)eta2(D[[v]],D$pre_state))
cat("η²:",paste(sprintf("%s=%.3f",names(e2),e2),collapse=" · "),"\n")
dip_p<-dip.test(P$PC1)$p.value; cat(sprintf("dip p=%.3f (unimodality not rejected)\n",dip_p))

cols3<-c(DRY="#e67e22",WET="#2980b9",SNOW="#7f8c8d")
figure_device(file.path(FIG,"fig2_conditions.png"),width=2600,height=1600,res=200,pointsize=19); panel_letters(letters[1:5])
layout(matrix(c(1,2,3,4,5,5),2,3,byrow=TRUE),widths=c(1,1,1))
par(mar=c(4.4,4.6,2.8,1),cex.lab=1.10,cex.axis=1.0,cex.main=1.15)

plot(NA,xlim=c(-8,4),ylim=c(0,0.5),xlab=expression("Road surface temperature at onset ("*degree*"C)"),
     ylab="Density",main="")
for(s in names(cols3)){d<-density(D[pre_state==s]$ts0,na.rm=TRUE);lines(d,col=cols3[s],lwd=2.4)}
abline(v=0,lty=3,col="grey40")
legend("topleft",names(cols3),col=cols3,lwd=2.4,bty="n",cex=1.05)

bp<-barplot(e2,names.arg=rep("",length(e2)),col="#b8d8e8",border="#2980b9",ylim=c(0,max(0.1,max(e2)*1.5)),
  ylab=expression("Proportion of variation explained, "*eta^2),
  main="",axes=FALSE)
axis(2,gap.axis=-1)   # draw every tick label; R versions differ in which overlapping labels they drop
text(bp,par("usr")[3]-0.04*diff(par("usr")[3:4]),c(expression(T[s0]),expression(RH[0]),"24-h range"),
     srt=30,adj=1,xpd=TRUE,cex=1.0)
text(bp,e2,sprintf("%.3f",e2),pos=3,cex=1.0)

tr<-lg[lag<=1440,.(md=median(ts,na.rm=TRUE),q1=quantile(ts,.25,na.rm=TRUE),
                   q3=quantile(ts,.75,na.rm=TRUE)),by=.(pre_state,lag)]
TRW<-dcast(tr,lag~pre_state,value.var="md"); TRW[,median_range_across_states:=pmax(DRY,SNOW,WET)-pmin(DRY,SNOW,WET)]
fwrite(TRW[order(-lag)],file.path(MART,"fig2c_trajectory_median.csv"),bom=TRUE)
plot(NA,xlim=c(-24,0),ylim=c(-6,10),xlab="Hours before onset",
     ylab=expression("Road surface temperature ("*degree*"C)"),main="")
abline(h=0,lty=3,col="grey40")
for(s in names(cols3)){g<-tr[pre_state==s][order(lag)]
  polygon(c(-g$lag/60,rev(-g$lag/60)),c(g$q1,rev(g$q3)),col=adjustcolor(cols3[s],.15),border=NA)
  lines(-g$lag/60,g$md,col=cols3[s],lwd=2.4)}
legend("top",names(cols3),col=cols3,lwd=2.4,bty="n",cex=1.0,horiz=TRUE,inset=0.005,seg.len=1.3,x.intersp=0.5,text.width=c(3.2,3.2,3.2))

hist(P$PC1,breaks=40,col="#b8d8e8",border="white",
  xlab=expression("PC1 (wet"%.%"near 0 "*degree*"C "%->%" dry"%.%"cold)"),
  main="")
box()

evx<-readRDS(file.path(MART,"intl_event_table.rds")); evx[,eid:=paste0(stn,"_",onset_min)]
P<-merge(P,evx[,.(eid,rh0,m_liq)],by="eid")
P[,pc_bin:=cut(PC1,quantile(PC1,seq(0,1,.1)),include.lowest=TRUE,labels=FALSE)]
g<-P[,.(pc=median(PC1),precip=mean(precip_frac,na.rm=TRUE),rh=mean(rh0,na.rm=TRUE),
        tso=mean(m_liq,na.rm=TRUE),cold=mean(log1p(cold_dose),na.rm=TRUE)),by=pc_bin][order(pc)]
norm01<-function(x)(x-min(x))/(max(x)-min(x))
plot(NA,xlim=range(g$pc),ylim=c(0,1.25),yaxt="n",xlab=expression("PC1 (wet"%.%"near 0 "*degree*"C "%->%" dry"%.%"cold)"),
     ylab="Normalized value",main="")
axis(2,seq(0,1,.25))
lines(g$pc,norm01(g$precip),type="b",pch=19,col="#2980b9",lwd=2)
lines(g$pc,norm01(g$rh),type="b",pch=17,col="#16a085",lwd=2)
lines(g$pc,norm01(g$tso),type="b",pch=15,col="#e67e22",lwd=2)
lines(g$pc,norm01(g$cold),type="b",pch=18,col="#c0392b",lwd=2)
lab_e<-c("Precipitation minutes (24 h)","RH at onset","Ts at onset","Cold dose (log)")
legend("top",lab_e,
       col=c("#2980b9","#16a085","#e67e22","#c0392b"),pch=c(19,17,15,18),lwd=2,
       bty="n",ncol=4,cex=0.88,x.intersp=0.55,seg.len=1.4,
       text.width=strwidth(lab_e,cex=0.88)+strwidth("MMM",cex=0.88),xpd=NA)
dev.off()
cat("Saved: 05_REPORT/FIGURES/fig2_conditions.png\n")
