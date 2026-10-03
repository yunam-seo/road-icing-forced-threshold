# Evaluates event-control AUC of RWIS-only and surrounding-field models (concurrent, 1 h and 3 h earlier) under station, date and winter holdout.
# Inputs: 02_MART/grid_full_raw.rds, grid_full_features.csv, grid_past_scale_raw.rds.
# Outputs: 02_MART/grid_holdout_auc.csv, grid_spatial_shap.csv.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages({library(data.table);library(lightgbm)})
MART<-file.path(ROOT, "02_MART")
source(code_file("_grid_features.R"))

FEAT_3KM<-c("ta_anom","td_anom","ta_grad","td_grad","ta_lap","td_lap","dd","dd_anom",
            "d_ta","d_td","d_dd","d_ta_anom","d_td_anom","adv_ta","adv_td")
B<-fread(file.path(MART,"grid_full_features.csv"),select=c("row","stn","label","rts","rta","rrh","rws",FEAT_3KM))
setnames(B,FEAT_3KM,paste0(FEAT_3KM,"_3km"))
raw<-readRDS(file.path(MART,"grid_full_raw.rds")); P<-as.data.table(raw$P)[,.(tm_min)]; P[,row:=.I]; rm(raw); gc(verbose=FALSE)
B<-merge(B,P,by="row")
X<-as.data.table(readRDS(file.path(MART,"grid_past_scale_raw.rds")))
stopifnot(all(c(0,30,60,90,180,210) %in% unique(X$lag)))

feat_at<-function(s){
  a<-X[lag==s]; b<-X[lag==s+30L,.(row,ta_vc_prev=ta_vc,td_vc_prev=td_vc)]
  G<-merge(a,b,by="row"); stopifnot(nrow(G)==nrow(B))
  F<-build_final_features(G, keys="row")
  cat(sprintf("  reference %3d min before: rows masked by wind gap >30 min %d (events %d)\n", s, sum(F$wind_masked),
              sum(F$wind_masked & G$label[match(F$row,G$row)]==1L)))
  F[,wind_masked:=NULL]
  setnames(F, FEAT_SPATIAL_FINAL, paste0(FEAT_SPATIAL_FINAL,"_s",s)); F }
cat("Assembling surrounding-field features\n")
M<-Reduce(function(a,b) merge(a,b,by="row"), list(B, feat_at(0L), feat_at(60L), feat_at(180L)))
frw<-c("rts","rta","rrh","rws")
M<-M[complete.cases(M[,..frw])]
M[,day:=as.integer(tm_min%/%1440L)]
source(code_file("00_analysis_window.R"))
M[,ptm:=as.POSIXct(tm_min*60,origin="1970-01-01",tz="UTC")]
M[,winter:=winter_label(as.integer(format(ptm,"%Y")),as.integer(format(ptm,"%m")))]; M[,ptm:=NULL]
print(M[,.(n=.N,n_events=sum(label),start_min=min(tm_min),end_min=max(tm_min)),by=winter][order(winter)])
WIN_OK<-2023:2025
cat(sprintf("Excluded outside the analysis window: %d rows (events %d)
", M[!winter %in% WIN_OK,.N], M[!winter %in% WIN_OK,sum(label)]))
M<-M[winter %in% WIN_OK]
cat(sprintf("Sample: %d rows (events %d · controls %d), stations %d, winters %s\n", nrow(M), sum(M$label), sum(M$label==0L),
            uniqueN(M$stn), paste(sort(unique(M$winter)),collapse="/")))

MODELS<-list(RWIS=frw,
             `concurrent surrounding fields`=c(frw,paste0(FEAT_SPATIAL_FINAL,"_s0")),
             `surrounding fields 1 h earlier`=c(frw,paste0(FEAT_SPATIAL_FINAL,"_s60")),
             `surrounding fields 3 h earlier`=c(frw,paste0(FEAT_SPATIAL_FINAL,"_s180")),
             initial_3km=c(frw,paste0(FEAT_3KM,"_3km")))
auc<-function(s,p){o<-is.finite(s);s<-s[o];p<-p[o];np<-as.numeric(sum(p));nn<-as.numeric(sum(!p))
  if(!np||!nn)return(NA_real_);r<-rank(s);(sum(r[p])-np*(np+1)/2)/(np*nn)}
pars<-list(objective="binary",learning_rate=0.05,num_leaves=15,min_data_in_leaf=20,
  feature_fraction=0.9,bagging_fraction=0.9,bagging_freq=1,verbosity=-1)

OUT<-list()
for(grp in c("stn","day","winter")){
  set.seed(2026)
  g<-M[[grp]]; ug<-sort(unique(g))
  folds<-if(grp=="winter") setNames(seq_along(ug),ug) else setNames(sample(rep(1:5,length.out=length(ug))),as.character(ug))
  kf<-folds[as.character(g)]
  for(mn in names(MODELS)){ fs<-MODELS[[mn]]; pr<-rep(NA_real_,nrow(M))
    for(k in sort(unique(kf))){ tr<-which(kf!=k); te<-which(kf==k)
      if(!sum(M$label[tr])||!sum(!M$label[tr])) next
      m<-lgb.train(pars,lgb.Dataset(as.matrix(M[tr,..fs]),label=M$label[tr]),nrounds=200)
      p<-predict(m,as.matrix(M[te,..fs])); if(is.matrix(p)) p<-p[,1]; pr[te]<-p }
    a<-auc(pr,M$label==1L)
    OUT[[length(OUT)+1L]]<-data.table(holdout=grp,model=mn,AUC=a,n=sum(is.finite(pr)),n_events=sum(M$label[is.finite(pr)]),n_folds=length(unique(kf)))
    cat(sprintf("[%-6s] %-8s AUC %.4f\n",grp,mn,a)) } }
OUT<-rbindlist(OUT)
OUT[,gain_vs_rwis_only:=AUC-AUC[model=="RWIS"],by=holdout]
print(OUT)
fwrite(OUT,file.path(MART,"grid_holdout_auc.csv"),bom=TRUE)
cat("Saved: 02_MART/grid_holdout_auc.csv\n")

fs<-MODELS$`concurrent surrounding fields`
m<-lgb.train(pars,lgb.Dataset(as.matrix(M[,..fs]),label=M$label),nrounds=200)
sh<-predict(m,as.matrix(M[,..fs]),type="contrib"); colnames(sh)<-c(fs,"bias")
SH<-data.table(feature=fs, mean_abs_shap=colMeans(abs(sh[,fs])))[order(-mean_abs_shap)]
SH[,feature:=sub("_s0$","",feature)]
print(SH)
fwrite(SH,file.path(MART,"grid_spatial_shap.csv"),bom=TRUE)
cat("Saved: 02_MART/grid_spatial_shap.csv\n")
