# Checks that the discrete (Euler) leaky-integrator approximation does not change the discrimination results, by comparison with the exact solution.
# Inputs: 02_MART/intl_cumstress.rds; 1-min RWIS files.
# Outputs: 02_MART/leaky_approx_check.csv, leaky_gap_distribution.csv.
# Reads the full 1-min archive (slow).

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
MART<-file.path(ROOT,"02_MART")
source(code_file("07_intl_cumstress.R"))

CS<-as.data.table(readRDS(file.path(MART,"intl_cumstress.rds")))[,.(stn,tm_min,type,cold07=cold_dose,dep07=dep_dose)]
files<-list.files(DDIR,pattern="\\.csv$",full.names=TRUE)
sof<-as.integer(sub(".*rwis_(\\d+)_.*","\\1",basename(files)))

accum<-function(tm,x,tau,gap,exact){
  n<-length(tm); dt<-c(0,diff(tm)); inwin<-!is.na(dt)&dt<=gap
  w<-if(exact) tau*(1-exp(-dt/tau)) else dt
  u<-ifelse(inwin & !is.na(x), x*w, 0)
  E<-cumsum(dt/tau)
  blk<-floor(E/500); S<-numeric(n); carry<-0; Eprev_end<-0
  for(b in unique(blk)){ ii<-which(blk==b); E0<-E[ii[1]]
    c0<-carry*exp(-(E0-Eprev_end))
    acc<-c0+cumsum(u[ii]*exp(E[ii]-E0))
    S[ii]<-acc*exp(-(E[ii]-E0))
    carry<-S[ii[length(ii)]]; Eprev_end<-E[ii[length(ii)]] }
  S }

R<-list(); GAP<-integer(16)
for(i in seq_along(files)){ sid<-sof[i]; q<-CS[stn==sid]; if(!nrow(q)) next
  d<-read_station(files[i]); if(is.null(d)||nrow(d)<10L) next
  dt<-diff(d$tm_min); GAP<-GAP+tabulate(pmin(dt,16L),16L)
  fp<-frost_point(d$ta,d$rh); xc<-pmax(0,-d$ts); xd<-pmax(0,-(d$ts-fp))
  out<-data.table(tm_min=d$tm_min)
  for(tau in c(60,1440)) for(ex in c(FALSE,TRUE)){
    tg<-sprintf("%s_%d",if(ex)"exact"else"euler",tau)
    out[[paste0("cold_",tg)]]<-accum(d$tm_min,xc,tau,cfg$gap_max,ex)
    out[[paste0("dep_",tg)]] <-accum(d$tm_min,xd,tau,cfg$gap_max,ex) }
  setkey(out,tm_min)
  R[[length(R)+1L]]<-cbind(q, out[.(q$tm_min),roll="nearest",on="tm_min"][,-"tm_min"])
  if(length(R)%%20L==0L) message("  ",length(R)," stations") }
R<-rbindlist(R)
cat(sprintf("Sample: onset %d · control %d · stations %d\n",sum(R$type=="onset"),sum(R$type=="control"),uniqueN(R$stn)))

cat("\n(1) Euler re-implementation vs values saved by script 07 (τ=1440)\n")
rel<-function(a,b){ m<-pmax(abs(a),abs(b)); ifelse(m>1e-9,abs(a-b)/m,0) }
e1<-max(rel(R$cold_euler_1440,R$cold07)); e2<-max(rel(R$dep_euler_1440,R$dep07))
cat(sprintf("  Max relative difference: cold dose %.2e · deposition dose %.2e  → %s\n",e1,e2,if(max(e1,e2)<1e-6)"reproduced" else "mismatch"))

auc<-function(s,p){r<-rank(s);np<-as.numeric(sum(p));nn<-as.numeric(sum(!p));(sum(r[p])-np*(np+1)/2)/(np*nn)}
pos<-R$type=="onset"
cat("\n(2) Euler vs exact solution\n")
T2<-rbindlist(lapply(c(60,1440),function(tau) rbindlist(lapply(c("cold","dep"),function(v){
  a<-R[[sprintf("%s_euler_%d",v,tau)]]; b<-R[[sprintf("%s_exact_%d",v,tau)]]
  data.table(indicator=v,tau_min=tau,auc_euler=auc(a,pos),auc_exact=auc(b,pos),
             auc_diff=auc(a,pos)-auc(b,pos),rank_corr=cor(a,b,method="spearman"),
             ratio_median=median((a/b)[b>1e-9]),ratio_max=max((a/b)[b>1e-9]),ratio_min=min((a/b)[b>1e-9])) }))))
print(T2)
cat(sprintf("\n(3) Observation intervals (consecutive rows): 1 min %.4f%% · 2–15 min %.4f%% · over 15 min %.4f%%\n",
            100*GAP[1]/sum(GAP),100*sum(GAP[2:15])/sum(GAP),100*GAP[16]/sum(GAP)))
fwrite(T2,file.path(MART,"leaky_approx_check.csv"),bom=TRUE)
fwrite(data.table(interval=c("1 min","2–15 min","over 15 min"),share_pct=100*c(GAP[1],sum(GAP[2:15]),GAP[16])/sum(GAP)),
       file.path(MART,"leaky_gap_distribution.csv"),bom=TRUE)
cat("Saved: 02_MART/leaky_approx_check.csv\n")
