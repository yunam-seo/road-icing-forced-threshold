# Draws Fig. 1: station map of both networks, icing onsets per winter and the analysis design.
# Inputs: meta/STN_AWS3_MTU_OBS.csv, stn_exclusion_intl.csv; 02_MART/intl_event_table.rds, icetarget_event_table.rds, data_availability.csv.
# Output: 05_REPORT/FIGURES/fig1_overview.png.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages({library(data.table);library(maps)})
META<-META_DIR
MART<-file.path(ROOT, "02_MART")
FIG<-file.path(ROOT, "05_REPORT/FIGURES")
source(code_file("_figure_style.R"))

stn<-fread(meta_file("STN_AWS3_MTU_OBS.csv"),select=c("STN_ID","LAT","LON"))
stn<-unique(stn,by="STN_ID")
DDIR<-rwis_dir()
allstn<-as.integer(sub(".*rwis_(\\d+)_.*","\\1",list.files(DDIR,pattern="\\.csv$")))
excl<-fread(meta_file("stn_exclusion_intl.csv"))$STN_ID
s121<-setdiff(allstn,excl)
s127<-unique(readRDS(file.path(MART,"icetarget_event_table.rds"))$stn)
stopifnot(length(intersect(s121,s127))==0)
cat("Primary",length(s121),"stations · independent",length(s127),"stations · overlap 0\n")
p121<-stn[STN_ID%in%s121]; p127<-stn[STN_ID%in%s127]
cat("Coordinates matched:",nrow(p121),"/",length(s121)," · ",nrow(p127),"/",length(s127),"\n")

wint<-function(rds){
  ev<-readRDS(rds)
  tt<-as.POSIXct(ev$onset_min*60,origin="1970-01-01",tz="UTC")
  y<-as.integer(format(tt,"%Y")); m<-as.integer(format(tt,"%m"))
  w<-ifelse(m>=11L,y,y-1L); table(factor(w,levels=2023:2025))
}
w121<-wint(file.path(MART,"intl_event_table.rds"))
w127<-wint(file.path(MART,"icetarget_event_table.rds"))
cat("By winter, 121:",w121," total",sum(w121),"\nBy winter, 127:",w127," total",sum(w127),"\n")

C1<-"#2980b9"; C2<-"#c0392b"
figure_device(file.path(FIG,"fig1_overview.png"),width=2400,height=1250,res=200,pointsize=14); panel_letters(letters[1:3])
layout(matrix(c(1,1,2,3),2,2),widths=c(1.0,1.45),heights=c(0.92,1.08))
par(mar=c(3.6,3.8,2.6,0.8),cex.lab=1.10,cex.axis=1.0,cex.main=1.15)

map("world",regions="South Korea",fill=TRUE,col="grey93",border="grey55",
    xlim=c(125.8,129.7),ylim=c(33.9,38.4),mar=c(4.2,3.8,2.6,0.8))
map.axes(cex.axis=1.0)
points(p121$LON,p121$LAT,pch=21,bg=adjustcolor(C1,.75),col="white",cex=1.25,lwd=.6)
points(p127$LON,p127$LAT,pch=24,bg=adjustcolor(C2,.75),col="white",cex=1.1,lwd=.6)
mtext("Longitude (\u00b0E)",side=1,line=2.3,cex=.85)
mtext("Latitude (\u00b0N)",side=2,line=2.5,cex=.85)
legend("bottomright",c(sprintf("Primary (%d stations)",length(s121)),
                    sprintf("Independent (%d stations)",length(s127))),
       pch=c(21,24),pt.bg=adjustcolor(c(C1,C2),.75),col="white",
       pt.cex=c(1.4,1.25),bty="o",bg="white",box.col="grey70",cex=1.0)

par(mar=c(3.8,4.4,2.6,0.8))
bp<-barplot(rbind(as.numeric(w121),as.numeric(w127)),beside=TRUE,
  names.arg=c("2023/24","2024/25","2025/26"),col=adjustcolor(c(C1,C2),.8),
  border=c(C1,C2),ylim=c(0,max(w127,w121)*1.3),ylab="Icing-onset events",
  main="")
text(bp,rbind(as.numeric(w121),as.numeric(w127)),
     format(rbind(as.numeric(w121),as.numeric(w127)),big.mark=",",trim=TRUE),pos=3,cex=.95)
legend("topleft",c(sprintf("Primary (total %s)",format(sum(w121),big.mark=",")),
                   sprintf("Independent (total %s)",format(sum(w127),big.mark=","))),
       fill=adjustcolor(c(C1,C2),.8),border=c(C1,C2),bty="n",cex=1.0)

par(mar=c(0.3,0.3,2.6,0.3))
plot(NA,xlim=c(-0.05,10.05),ylim=c(0,10),axes=FALSE,xlab="",ylab="",main="",xaxs="i",yaxs="i")
box_<-function(x0,x1,y0,y1,txt,col="white",border="grey35",cex=.92,fcol="black",lwd=1.4){
  rect(x0,y0,x1,y1,col=col,border=border,lwd=lwd)
  text((x0+x1)/2,(y0+y1)/2,txt,cex=cex,col=fcol)}
arr<-function(x0,y0,x1,y1,col="grey30",lty=1) arrows(x0,y0,x1,y1,length=.075,angle=22,lwd=1.5,col=col,lty=lty)
YA<-c(7.95,9.95); YB<-c(4.75,7.10); YC<-c(0.25,3.05)
DA<-fread(file.path(MART,"data_availability.csv"),encoding="UTF-8")
tot_min<-as.numeric(DA[item=="total observed minutes"]$value); stopifnot(length(tot_min)==1)
box_(0.10,7.30,YA[1],YA[2],sprintf("Primary network: %d RWIS stations, 1-min observations\n%s winters · %.1f million minutes · %s icing onsets",
     length(s121),c("one","two","three","four")[length(w121)],tot_min/1e6,format(sum(w121),big.mark=",")),
     col=adjustcolor(C1,.12),border=C1)
box_(7.55,9.90,YA[1],YA[2],sprintf("Independent network\n%d stations\n%s icing onsets",length(s127),format(sum(w127),big.mark=",")),
     col=adjustcolor(C2,.12),border=C2,cex=.86)
XB<-list(c(0.10,3.20),c(3.40,6.50),c(6.70,9.90))
box_(XB[[1]][1],XB[[1]][2],YB[1],YB[2],"Conditions at\nicing onset\n(continuity vs. pathways)")
box_(XB[[2]][1],XB[[2]][2],YB[1],YB[2],"Cumulative antecedent\nconditions\n(leaky integrator)")
box_(XB[[3]][1],XB[[3]][2],YB[1],YB[2],"Critical slowing\ndown (CSD) indicators\n(variance, AR1)")
for(x in c(1.65,4.95,7.00)) arr(x,YA[1],x,YB[2])
arr(8.72,YA[1],8.72,YB[2],col=C2,lty=2)
text(8.84,(YA[1]+YB[2])/2,"replication",adj=0,cex=.80,col=C2,font=3)
box_(0.10,4.85,YC[1],YC[2],"",col="grey95")
text(2.475,2.30,"Road-icing prediction model\ncurrent, recent and cumulative conditions",cex=.92)
text(2.475,0.95,"also applied to the independent network\nwithout retraining",cex=.80,col=C2,font=3)
box_(5.35,9.90,YC[1],YC[2],"+ surrounding meteorological fields\n(500 m objectively analyzed fields)",col="grey95")
arr(1.65,YB[1],1.65,YC[2]); arr(4.30,YB[1],4.30,YC[2])
arr(4.85,1.675,5.35,1.675)
par(mfrow=c(1,1));dev.off()
cat("Saved: 05_REPORT/FIGURES/fig1_overview.png\n")
