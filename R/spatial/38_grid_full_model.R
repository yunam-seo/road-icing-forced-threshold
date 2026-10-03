# Assembles grid pattern features and compares station-held-out AUC of RWIS-only and pattern-augmented models.
# Input: 02_MART/grid_full_raw.rds.
# Outputs: 02_MART/grid_full_features.csv; 05_REPORT/FIGURES/intl_grid_pattern.png.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages({library(data.table);library(lightgbm)})
set.seed(2026)
MART<-file.path(ROOT, "02_MART"); FIG<-file.path(ROOT, "05_REPORT/FIGURES")
x<-readRDS(file.path(MART,"grid_full_raw.rds")); P<-as.data.table(x$P); S<-as.data.table(x$S); rwis<-x$rwis
P[,row:=.I]

w<-dcast(S, row ~ gv + lag, value.var=c("vc","vreg","gx","gy","lap"))
need<-c("vc_ta_0","vc_ta_30","vc_td_0","vc_td_30","vreg_ta_0","vreg_td_0",
        "gx_ta_0","gy_ta_0","gx_td_0","gy_td_0","lap_ta_0","lap_td_0",
        "vreg_ta_30","vreg_td_30")
for(c0 in need) if(!c0 %in% names(w)) w[[c0]]<-NA_real_

rv<-rbindlist(lapply(seq_len(nrow(P)),function(i){
  d<-rwis[[as.character(P$stn[i])]]; if(is.null(d))return(data.table(row=P$row[i]))
  r<-d[.(P$tm_min[i]),roll="nearest",on="tm"]
  data.table(row=P$row[i], rts=r$ts, rta=r$ta, rrh=r$rh, rws=r$ws, rwd=r$wd)}),fill=TRUE)
D<-Reduce(function(a,b)merge(a,b,by="row",all.x=TRUE), list(P[,.(row,stn,label)], w, rv))

D[,`:=`(
  ta_anom=vc_ta_0-vreg_ta_0, td_anom=vc_td_0-vreg_td_0,
  ta_grad=sqrt(gx_ta_0^2+gy_ta_0^2), td_grad=sqrt(gx_td_0^2+gy_td_0^2),
  ta_lap=lap_ta_0, td_lap=lap_td_0,
  dd=vc_ta_0-vc_td_0, dd_anom=(vc_ta_0-vc_td_0)-(vreg_ta_0-vreg_td_0),
  d_ta=vc_ta_0-vc_ta_30, d_td=vc_td_0-vc_td_30,
  d_dd=(vc_ta_0-vc_td_0)-(vc_ta_30-vc_td_30),
  d_ta_anom=(vc_ta_0-vreg_ta_0)-(vc_ta_30-vreg_ta_30),
  d_td_anom=(vc_td_0-vreg_td_0)-(vc_td_30-vreg_td_30))]
D[,u:=-rws*sin(rwd*pi/180)]; D[,v:=-rws*cos(rwd*pi/180)]
D[,adv_ta:=u*gx_ta_0+v*gy_ta_0]; D[,adv_td:=u*gx_td_0+v*gy_td_0]

feat_rwis<-c("rts","rta","rrh","rws")
feat_pt<-c("vc_ta_0","vc_td_0")
feat_pat<-c("ta_anom","td_anom","ta_grad","td_grad","ta_lap","td_lap",
            "dd","dd_anom","d_ta","d_td","d_dd","d_ta_anom","d_td_anom","adv_ta","adv_td")
auc<-function(s,p){o<-is.finite(s);s<-s[o];p<-p[o];np<-as.numeric(sum(p));nn<-as.numeric(sum(!p))
  if(!np||!nn)return(NA);r<-rank(s);(sum(r[p])-np*(np+1)/2)/(np*nn)}
stns<-unique(D$stn); fold<-setNames(sample(rep(1:5,length.out=length(stns))),stns)
D[,fold:=fold[as.character(stn)]]
pars<-list(objective="binary",learning_rate=0.05,num_leaves=15,min_data_in_leaf=20,
  feature_fraction=0.9,bagging_fraction=0.9,bagging_freq=1,verbosity=-1)
cv<-function(fs){fs<-intersect(fs,names(D))
  Dc<-D[complete.cases(D[,..fs])]; pr<-rep(NA_real_,nrow(Dc))
  for(kf in 1:5){tr<-which(Dc$fold!=kf);te<-which(Dc$fold==kf)
    if(!length(te)||!sum(Dc$label[tr]))next
    m<-lgb.train(pars,lgb.Dataset(as.matrix(Dc[tr,..fs]),label=Dc$label[tr]),nrounds=200)
    p<-predict(m,as.matrix(Dc[te,..fs]));if(is.matrix(p))p<-p[,1];pr[te]<-p}
  list(auc=auc(pr,Dc$label==1L), n=nrow(Dc))}

cat("Sample:",nrow(D)," cases:",sum(D$label),"\n=== onset vs control station-held-out 5-fold CV AUC ===\n")
res<-list()
for(nm in list(c("RWIS-only",feat_rwis), c("RWIS+point values",c(feat_rwis,feat_pt)),
               c("RWIS+patterns",c(feat_rwis,feat_pat)), c("full",c(feat_rwis,feat_pt,feat_pat)),
               c("patterns only",feat_pat))){
  r<-cv(nm[-1]); res[[nm[1]]]<-r
  cat(sprintf("  %-14s AUC %.3f  (n=%d, features %d)\n",nm[1],r$auc,r$n,length(intersect(nm[-1],names(D))))) }
cat(sprintf("\nAdded value: full-RWIS = %+.3f,  +patterns-RWIS = %+.3f\n",
    res[["full"]]$auc-res[["RWIS-only"]]$auc, res[["RWIS+patterns"]]$auc-res[["RWIS-only"]]$auc))
cat("\nUnivariate pattern AUC (top 6):\n")
ua<-sort(sapply(feat_pat,function(v)auc(D[[v]],D$label==1L)),decreasing=TRUE)
for(nm in names(head(ua,6))) cat(sprintf("  %-12s %.3f\n",nm,ua[nm]))

png(file.path(FIG,"intl_grid_pattern.png"),width=1500,height=900,res=160)
par(mar=c(7,5,3,1.5))
labs<-c("RWIS-only","RWIS+point values","RWIS+patterns","full","patterns only")
av<-sapply(labs,function(l)res[[l]]$auc)
bp<-barplot(av,names.arg=NA,ylim=c(0.5,max(av,na.rm=TRUE)+0.03),col="#3b6fb0",
  ylab="Station-held-out 5-fold CV AUC",main="Added value of surrounding fields for icing nowcasting (all winters)")
text(bp,av+0.008,sprintf("%.3f",av),cex=1.0,font=2)
text(bp,par("usr")[3]-0.012,labs,srt=30,adj=1,xpd=TRUE,cex=0.95)
abline(h=av["RWIS-only"],lty=2,col="#c0392b")
dev.off(); cat("\nFigure: 05_REPORT/FIGURES/intl_grid_pattern.png\n")
fwrite(D, file.path(MART,"grid_full_features.csv"))
cat("Saved: 02_MART/grid_full_features.csv\n")
