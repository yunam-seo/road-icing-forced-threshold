# Matches each RWIS station to its 500 m objective-analysis grid cell.
# Inputs: meta/STN_AWS3_MTU_OBS.csv; grid coordinate table (RICE_GRID_LATLON); 02_MART/intl_event_table.rds.
# Output: 02_MART/road_grid500_match.csv.

local({a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f)) else getwd()
  while (!file.exists(file.path(d, "_config.R")) && dirname(d) != d) d <- dirname(d)
  assign("CODE", d, envir = globalenv()); source(file.path(d, "_config.R"))})
suppressPackageStartupMessages(library(data.table))
OUT<-file.path(ROOT, "02_MART")

ev<-readRDS(file.path(ROOT, "02_MART/intl_event_table.rds"))
stns<-sort(unique(ev$stn))
meta<-fread(meta_file("STN_AWS3_MTU_OBS.csv"),
            select=c("STN_ID","LAT","LON"),encoding="Latin-1")
meta<-unique(meta[,.(STN_ID,LAT=as.numeric(LAT),LON=as.numeric(LON))])
meta<-meta[!is.na(LAT)&!is.na(LON)]
m<-meta[STN_ID %in% stns]
cat("Event stations:",length(stns)," coordinates matched:",nrow(m)," missing:",length(setdiff(stns,m$STN_ID)),"\n")

g<-fread(grid_latlon())
cat("Grid points:",nrow(g)," lon range",paste(round(range(g$lon),1),collapse="~"),
    " lat",paste(round(range(g$lat),1),collapse="~"),"\n")

haversine<-function(lon1,lat1,lon2,lat2){R<-6371;d2r<-pi/180
  dlon<-(lon2-lon1)*d2r;dlat<-(lat2-lat1)*d2r
  a<-sin(dlat/2)^2+cos(lat1*d2r)*cos(lat2*d2r)*sin(dlon/2)^2;2*R*asin(sqrt(a))}

res<-rbindlist(lapply(seq_len(nrow(m)),function(i){
  la<-m$LAT[i];lo<-m$LON[i]
  sub<-g[abs(lat-la)<0.05 & abs(lon-lo)<0.06]
  if(!nrow(sub)) sub<-g[abs(lat-la)<0.15 & abs(lon-lo)<0.18]
  d<-haversine(lo,la,sub$lon,sub$lat);j<-which.min(d)
  data.table(STN_ID=m$STN_ID[i],lat=la,lon=lo,
             xii=sub$xii[j],yii=sub$yii[j],grid_lon=sub$lon[j],grid_lat=sub$lat[j],
             dist_km=round(d[j],3))
}))
cat("Mapping done:",nrow(res)," max distance km:",round(max(res$dist_km),3),
    " median:",round(median(res$dist_km),3),"\n")
fwrite(res,file.path(OUT,"road_grid500_match.csv"))
cat("Saved: 02_MART/road_grid500_match.csv\n"); print(head(res,4))
