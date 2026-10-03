# Clusters pre-onset trajectories ({ts, m_dep, rh} at 14 lags) with k-means and compares clusters with the prior surface state.
# Inputs: 02_MART/intl_lagged_table.rds, intl_event_table.rds.
# Outputs: 02_MART/intl_cluster_assign.rds, intl_cluster_profiles.csv; 05_REPORT/FIGURES/intl_pathway_clusters.png.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages({
  library(data.table); library(cluster)
})
set.seed(2026)
MART <- file.path(ROOT, "02_MART")
FIG  <- file.path(ROOT, "05_REPORT/FIGURES")
dir.create(FIG, showWarnings=FALSE, recursive=TRUE)

cfg <- list(vars=c("ts","m_dep","rh"), cov_thr=0.7, ks=2:6)

lg <- readRDS(file.path(MART,"intl_lagged_table.rds"))
lags <- sort(unique(lg$lag))
setorder(lg, eid, lag)

fill_event <- function(v, x) {
  if (all(is.na(v))) return(v)
  approx(x[!is.na(v)], v[!is.na(v)], xout=x, rule=2)$y
}
cov <- lg[, .(cov=mean(!is.na(ts))), by=eid]
keep <- cov[cov>=cfg$cov_thr, eid]
cat("Events:", uniqueN(lg$eid), " coverage>=", cfg$cov_thr, ":", length(keep), "\n")
d <- lg[eid %in% keep]
for (v in cfg$vars) d[, (v):=fill_event(get(v), lag), by=eid]

w <- dcast(melt(d, id.vars=c("eid","lag"), measure.vars=cfg$vars),
           eid ~ variable + lag, value.var="value")
feat <- as.matrix(w[, -1]); rownames(feat) <- w$eid
feat <- feat[complete.cases(feat), , drop=FALSE]
Z <- scale(feat)
cat("Clustering input:", nrow(Z), "events ×", ncol(Z), "features\n")

sil <- sapply(cfg$ks, function(k){
  km <- kmeans(Z, centers=k, nstart=25, iter.max=50)
  mean(silhouette(km$cluster, dist(Z))[,3])
})
kbest <- cfg$ks[which.max(sil)]
cat("Silhouette:", paste(sprintf("k%d=%.3f",cfg$ks,sil),collapse=" "), "→ best k=", kbest, "\n")
km <- kmeans(Z, centers=kbest, nstart=50, iter.max=100)

asg <- data.table(eid=rownames(Z), cluster=km$cluster)
ev  <- readRDS(file.path(MART,"intl_event_table.rds"))
ev[, eid := paste0(stn,"_",onset_min)]
asg <- merge(asg, ev[, .(eid, stn, pre_state, ts0, ta0, rh0, ice_dur, m_liq, year,
                         onset_tm)], by="eid")
asg[, region := substr(as.character(stn),1,2)]
asg[, month  := data.table::month(onset_tm)]

ari <- function(a,b){
  t<-table(a,b); n<-sum(t); a_<-rowSums(t); b_<-colSums(t)
  idx<-sum(choose(t,2)); e<-sum(choose(a_,2))*sum(choose(b_,2))/choose(n,2)
  m<-(sum(choose(a_,2))+sum(choose(b_,2)))/2; (idx-e)/(m-e)
}
cat("\nARI (cluster vs pre_state):", round(ari(asg$cluster, asg$pre_state),3), "\n")
cat("\nCluster sizes:\n"); print(table(asg$cluster))
cat("\nCluster × prior state:\n"); print(table(asg$cluster, asg$pre_state))
cat("\nOnset characteristics by cluster (median):\n")
print(asg[, .(n=.N, ts0=round(median(ts0,na.rm=T),1), ta0=round(median(ta0,na.rm=T),1),
              rh0=round(median(rh0,na.rm=T),1), dur=round(median(ice_dur),0)), by=cluster][order(cluster)])

prof <- d[eid %in% asg$eid]
prof <- merge(prof, asg[,.(eid,cluster)], by="eid")
pm <- prof[, lapply(.SD, mean, na.rm=TRUE), by=.(cluster,lag), .SDcols=cfg$vars]
png(file.path(FIG,"intl_pathway_clusters.png"), width=1100, height=380, res=110)
op<-par(mfrow=c(1,3), mar=c(4,4,2,1))
cols <- rainbow(kbest)
for (v in cfg$vars){
  rng<-range(pm[[v]],na.rm=TRUE)
  plot(NA, xlim=rev(range(lags)), ylim=rng, xlab="Minutes before onset (min)", ylab=v,
       main=v)
  if (v %in% c("ts","m_dep")) abline(h=0, lty=3, col="grey50")
  for (k in 1:kbest) { s<-pm[cluster==k]; lines(s$lag, s[[v]], col=cols[k], lwd=2) }
  if (v==cfg$vars[1]) legend("topright", legend=paste0("C",1:kbest),
                             col=cols, lwd=2, bty="n", cex=.8)
}
par(op); dev.off()

saveRDS(asg, file.path(MART,"intl_cluster_assign.rds"))
fwrite(pm,  file.path(MART,"intl_cluster_profiles.csv"))
cat("\nSaved: intl_cluster_assign.rds, intl_cluster_profiles.csv, 05_REPORT/FIGURES/intl_pathway_clusters.png\n")
